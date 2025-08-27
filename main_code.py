# Author: Chris Iliadis
# Host: Newcastle University

"""
Global SSD generator (Synthetic Storm Drains) — kerb-side + staggered + min-distance + CityCAT export

- Configure paths/CRS + CityCAT header in CONFIG (no CLI needed).
- Downloads roads with OSMnx (or uses your own roads).
- Places SSDs on both kerb lines, staggered, then enforces a global min spacing.
- Excludes SSDs inside buildings (optional).
- Exports:
    * Shapefiles: ssd_candidates.shp, ssds_final.shp (+ optional QA layers)
    * CityCAT text file: Inlets.txt (first row = header numbers, then id x y 1)
"""

import sys
from pathlib import Path
import json
import warnings

import geopandas as gpd
from shapely.geometry import LineString, MultiLineString
from shapely.ops import linemerge
import numpy as np
import pandas as pd

# ---------- Version-safe union helper (Shapely 1.8/2.x) ----------
try:
    from shapely import union_all as _union_all  # Shapely >= 2
    def union_all_geoms(geoms):
        return _union_all(list(geoms))
except Exception:
    from shapely.ops import unary_union as _unary_union  # Shapely < 2
    def union_all_geoms(geoms):
        return _unary_union(list(geoms))

# ========================= CONFIG =========================
# REQUIRED: study polygon
STUDY_AREA_PATH = r"copy the shape file of the study area"

# OPTIONAL: buildings polygons (for excluding SSDs within footprints)
BUILDINGS_PATH  = r"copy the shape file of the buildings"
# BUILDINGS_PATH = None

# OPTIONAL: use your own roads instead of OSM download
ROADS_PATH = None
# ROADS_PATH = r"D:\data\roads.gpkg"

# OUTPUT folder
OUTPUT_DIR = Path(r"copy the folder path for your outputs")

# CRS in metres: set to None for auto-UTM, or provide "EPSG:2100" / "EPSG:27700", etc.
TARGET_CRS = "add the projected coordinate system (e.g. EPSG:27700)"   # or None

# Export QA layers (buffered polygons + backbone lines)?
EXPORT_BACKBONE = True

# --- Placement parameters (metres) ---
# Per-class SSD spacing (this is the "global spacing" between consecutive inlets along the road).
SPACING_MAP = {
    "motorway": 60, "motorway_link": 60,
    "trunk": 60, "trunk_link": 60,
    "primary": 50, "primary_link": 50,
    "secondary": 50, "secondary_link": 50,
    "tertiary": 50, "tertiary_link": 50,
    "_default": 50
}

# Per-class kerb offsets from centerline to each side
KERB_OFFSET_MAP = {
    "motorway": 5.0, "motorway_link": 5.0,
    "trunk": 4.0, "trunk_link": 4.0,
    "primary": 3.5, "primary_link": 3.5,
    "secondary": 3.0, "secondary_link": 3.0,
    "tertiary": 3.0, "tertiary_link": 3.0,
    "_default": 3.0
}

# Global minimum distance between ANY two SSD points (metres)
MIN_SSD_SPACING = 50.0

# Road class → buffer (for QA polygons)
BUFFERS_M = {
    "motorway": 25, "motorway_link": 25,
    "trunk": 20, "trunk_link": 20,
    "primary": 20, "primary_link": 20,
    "secondary": 10, "secondary_link": 10,
    "tertiary": 10, "tertiary_link": 10,
    # Optionally include locals:
    # "residential": 8, "unclassified": 8, "service": 6, "living_street": 6,
}

# --- CityCAT export configuration ---
# First row numbers (e.g., inlet size(s) and friction coefficient) – adjust as needed
CITYCAT_HEADER = (0.30, 0.30, 0.50)   # example matches your format
CITYCAT_FILENAME = "Inlets.txt"
CITYCAT_COORD_DECIMALS = 4            # number of decimals for x/y in txt
CITYCAT_HEADER_DECIMALS = 2           # decimals for first-row numbers
CITYCAT_ACTIVE_FLAG = 1               # last column value
# =========================================================


def ensure_projected(gdf: gpd.GeoDataFrame, name: str):
    if gdf.crs is None:
        raise ValueError(f"{name} CRS is undefined.")
    if not gdf.crs.is_projected:
        raise ValueError(f"{name} CRS is geographic (degrees). Reproject to a projected CRS in metres.")

def auto_target_crs(study_gdf: gpd.GeoDataFrame) -> str:
    """Pick a local UTM based on centroid (pyproj>=3 suggested)."""
    try:
        from pyproj import CRS
        centroid = study_gdf.to_crs(4326).unary_union.centroid
        crs = CRS.estimate_utm_crs(latlon=(centroid.y, centroid.x))
        return crs.to_authority()[-1]  # 'EPSG:XXXX'
    except Exception:
        warnings.warn("Auto UTM failed; falling back to EPSG:3857 (metres, but distorted).")
        return "EPSG:3857"

def download_roads_osm(study_gdf_4326: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Download OSM roads inside the study polygon (expects 4326)."""
    try:
        import osmnx as ox
    except ImportError:
        raise RuntimeError("osmnx not installed. Install with: pip install osmnx")

    polygon = union_all_geoms(study_gdf_4326.geometry.values)
    G = ox.graph_from_polygon(polygon, network_type="drive")
    edges = ox.graph_to_gdfs(G, nodes=False, edges=True)
    roads = edges[["geometry", "highway"]].copy()
    # Clip to be safe
    roads = gpd.overlay(roads, study_gdf_4326, how="intersection")
    return roads  # still 4326

def normalize_highway(roads: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    roads = roads.copy()
    if "highway" in roads.columns:
        roads["__hw__"] = roads["highway"].apply(
            lambda v: ",".join(v).lower() if isinstance(v, (list, tuple, set)) else str(v).lower()
        )
    else:
        roads["__hw__"] = ""
    return roads

def buffer_roads(roads_proj: gpd.GeoDataFrame, buffers_map: dict) -> gpd.GeoDataFrame:
    """Buffered roads merged to a single polygon layer (QA)."""
    ensure_projected(roads_proj, "Roads")
    roads = normalize_highway(roads_proj)

    parts = []
    for klass, buf in buffers_map.items():
        sel = roads[roads["__hw__"].str.contains(klass, na=False)]
        cnt = len(sel)
        print(f"[buffer] class='{klass}' matched {cnt} features; buf={buf} m")
        if cnt > 0:
            parts.append(sel.buffer(buf, cap_style=2, join_style=2))

    if not parts:
        return gpd.GeoDataFrame(geometry=[], crs=roads.crs)

    merged = pd.concat(parts, ignore_index=True)
    dissolved = union_all_geoms(merged.geometry.values)
    if dissolved is None or dissolved.is_empty:
        return gpd.GeoDataFrame(geometry=[], crs=roads.crs)

    gser = gpd.GeoSeries(dissolved, crs=roads.crs).explode(index_parts=False).buffer(0)
    return gpd.GeoDataFrame(geometry=gser, crs=roads.crs)

# --- Kerb-side helpers ---

def _iter_lines(geom):
    """Yield LineStrings from LineString or MultiLineString."""
    if geom is None or geom.is_empty:
        return
    if isinstance(geom, LineString):
        yield geom
    elif isinstance(geom, MultiLineString):
        for g in geom.geoms:
            if isinstance(g, LineString):
                yield g
    else:
        try:
            m = linemerge(geom)
            yield from _iter_lines(m)
        except Exception:
            return

def _offset_lines(line: LineString, offset: float, side: str):
    """Create a parallel offset line; returns list (offset may produce MultiLineString)."""
    try:
        off = line.parallel_offset(distance=offset, side=side, join_style=2, mitre_limit=5.0)
        if off is None or off.is_empty:
            return []
        if isinstance(off, LineString):
            return [off]
        elif isinstance(off, MultiLineString):
            return [g for g in off.geoms if isinstance(g, LineString)]
        else:
            return []
    except Exception:
        return []

def points_on_kerb_sides_staggered(roads_proj: gpd.GeoDataFrame, spacing_map: dict, kerb_map: dict) -> gpd.GeoDataFrame:
    """
    Place SSD points on BOTH kerb sides with a staggered pattern:
    even indices on left, odd on right. Step along the road by spacing/2 so
    consecutive SSDs (alternating sides) are ~spacing apart globally.
    """
    ensure_projected(roads_proj, "Roads")
    roads = normalize_highway(roads_proj)

    pts, attrs = [], []
    covered = pd.Index([])

    # Explicit classes
    for klass, spacing in spacing_map.items():
        if klass.startswith("_"):
            continue
        kerb = float(kerb_map.get(klass, kerb_map.get("_default", 3.0)))
        sel = roads[roads["__hw__"].str.contains(klass, na=False)]
        if sel.empty:
            continue
        covered = covered.union(sel.index)

        for geom in sel.geometry:
            for line in _iter_lines(geom):
                left_lines  = _offset_lines(line, kerb, "left")
                right_lines = _offset_lines(line, kerb, "right")
                if not left_lines or not right_lines:
                    continue
                left, right = left_lines[0], right_lines[0]
                L = min(left.length, right.length)
                if L <= 0:
                    continue

                # stagger: positions every spacing/2; alternate sides
                step = spacing / 2.0
                dists = np.arange(step, L, step)
                for i, d in enumerate(dists):
                    side = "left" if i % 2 == 0 else "right"
                    line_use = left if side == "left" else right
                    if d < line_use.length:
                        pts.append(line_use.interpolate(d))
                        attrs.append({"class": klass, "side": side, "spacing_m": spacing, "kerb_m": kerb})

    # Remaining segments (default spacing/kerb)
    default_spacing = float(spacing_map.get("_default", 50))
    default_kerb    = float(kerb_map.get("_default", 3.0))
    rest = roads[~roads.index.isin(covered)]
    for geom in rest.geometry:
        for line in _iter_lines(geom):
            left_lines  = _offset_lines(line, default_kerb, "left")
            right_lines = _offset_lines(line, default_kerb, "right")
            if not left_lines or not right_lines:
                continue
            left, right = left_lines[0], right_lines[0]
            L = min(left.length, right.length)
            if L <= 0:
                continue

            step = default_spacing / 2.0
            dists = np.arange(step, L, step)
            for i, d in enumerate(dists):
                side = "left" if i % 2 == 0 else "right"
                line_use = left if side == "left" else right
                if d < line_use.length:
                    pts.append(line_use.interpolate(d))
                    attrs.append({"class": "other", "side": side, "spacing_m": default_spacing, "kerb_m": default_kerb})

    if not pts:
        return gpd.GeoDataFrame(columns=["class", "side", "spacing_m", "kerb_m", "geometry"], crs=roads.crs)

    return gpd.GeoDataFrame(attrs, geometry=pts, crs=roads.crs)

# --- Global min-distance thinning ---

def enforce_min_distance(points_gdf: gpd.GeoDataFrame, min_dist: float) -> gpd.GeoDataFrame:
    """
    Greedy thinning to enforce a global minimum spacing between SSD points.
    Sorts by (x, y) for determinism; uses STRtree when available.
    """
    if points_gdf.empty:
        return points_gdf

    pts = points_gdf.copy().reset_index(drop=True)
    pts["_x"] = pts.geometry.x
    pts["_y"] = pts.geometry.y
    pts = pts.sort_values(by=["_x", "_y"]).reset_index(drop=True)

    try:
        # Shapely STRtree fast path
        from shapely.strtree import STRtree
        geoms = list(pts.geometry.values)
        tree = STRtree(geoms)
        id_to_idx = {id(g): i for i, g in enumerate(geoms)}

        kept_mask = np.zeros(len(pts), dtype=bool)
        kept_indices = []

        for i, g in enumerate(geoms):
            if kept_mask[i]:
                continue
            # Check against already-kept neighbors within min_dist bbox
            neighbors = tree.query(g.buffer(min_dist))
            too_close = False
            for h in neighbors:
                j = id_to_idx.get(id(h), None)
                if j is not None and kept_mask[j]:
                    if g.distance(h) < min_dist:
                        too_close = True
                        break
            if not too_close:
                kept_mask[i] = True
                kept_indices.append(i)

        out = pts.iloc[kept_indices].drop(columns=["_x", "_y"])
        out.reset_index(drop=True, inplace=True)
        return out

    except Exception:
        # Fallback O(n^2)
        kept_idx = []
        kept_pts = []
        for i, g in enumerate(pts.geometry):
            ok = True
            for h in kept_pts:
                if g.distance(h) < min_dist:
                    ok = False
                    break
            if ok:
                kept_idx.append(i)
                kept_pts.append(g)
        out = pts.iloc[kept_idx].drop(columns=["_x", "_y"])
        out.reset_index(drop=True, inplace=True)
        return out

# --- CityCAT export ---

def write_citycat_inlets(points_gdf: gpd.GeoDataFrame,
                         out_path: Path,
                         header_triplet=(0.30, 0.30, 0.50),
                         coord_decimals: int = 4,
                         header_decimals: int = 2,
                         active_flag: int = 1,
                         id_start: int = 0) -> None:
    """
    Write CityCAT inlet file:
      - First row: three numbers (e.g., inlet size(s), friction coeff)
      - Next rows: id  x  y  active_flag
    """
    if points_gdf.empty:
        raise ValueError("No points to write to CityCAT file.")

    # Deterministic ordering (x then y)
    pts = points_gdf.copy()
    pts["_x"] = pts.geometry.x
    pts["_y"] = pts.geometry.y
    pts = pts.sort_values(by=["_x", "_y"]).reset_index(drop=True)

    fmt_h = f"{{:.{header_decimals}f}} {{:.{header_decimals}f}} {{:.{header_decimals}f}}\n"
    fmt_p = f"{{:d}} {{:.{coord_decimals}f}} {{:.{coord_decimals}f}} {{:d}}\n"

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(fmt_h.format(*header_triplet))
        for i, row in pts.iterrows():
            idx = id_start + i
            f.write(fmt_p.format(idx, float(row["_x"]), float(row["_y"]), active_flag))

# --- Building exclusion ---

def exclude_points_in_buildings(points: gpd.GeoDataFrame, buildings_path: str | None, target_crs: str) -> gpd.GeoDataFrame:
    """Remove points inside building footprints (sjoin with fallback)."""
    if not buildings_path or not Path(buildings_path).exists():
        return points
    bld = gpd.read_file(buildings_path)
    if bld.crs is None:
        raise ValueError("Buildings CRS is undefined.")
    bld = bld.to_crs(target_crs)
    try:
        inside = gpd.sjoin(points, bld, predicate="within", how="left")
        keep = inside[inside.index_right.isna()].drop(columns=["index_right"])
        return gpd.GeoDataFrame(keep, geometry="geometry", crs=points.crs)
    except Exception:
        # Fallback: tiny buffer + overlay if spatial index unavailable
        eps = 0.01
        pts_eps = points.copy()
        pts_eps["geometry"] = pts_eps.buffer(eps)
        dropped = gpd.overlay(pts_eps, bld, how="intersection")
        ids_to_drop = set(dropped.index)
        keep = points[~points.index.isin(ids_to_drop)]
        return gpd.GeoDataFrame(keep, geometry="geometry", crs=points.crs)

# --- Main pipeline ---

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1) Study area
    study = gpd.read_file(STUDY_AREA_PATH)
    if study.crs is None:
        raise ValueError("Study area CRS is undefined.")
    target_crs = TARGET_CRS or auto_target_crs(study)
    study = study.to_crs(target_crs)
    ensure_projected(study, "Study area")
    print(f"[info] Target CRS: {target_crs}")

    # 2) Roads (download or user-provided), clipped to study, in target CRS
    if ROADS_PATH:
        roads = gpd.read_file(ROADS_PATH)
        if roads.crs is None:
            roads = roads.set_crs(study.crs)
        roads = roads.to_crs(study.crs)
        roads = gpd.overlay(roads, study, how="intersection")
    else:
        roads_4326 = study.to_crs(4326)
        roads = download_roads_osm(roads_4326).to_crs(study.crs)

    if roads.empty:
        raise RuntimeError("No roads found after download/clip. Check study area geometry/extent and OSM coverage.")

    print(f"[info] Roads loaded: {len(roads)}")
    if "highway" in roads.columns:
        hw = roads["highway"].apply(lambda v: ",".join(v) if isinstance(v, (list, tuple, set)) else str(v))
        print("[info] Unique highway tags (sample):", sorted(hw.str.lower().unique())[:50])

    # 3) Kerb-side SSD placement (staggered, both sides)
    ssds_candidates = points_on_kerb_sides_staggered(roads, SPACING_MAP, KERB_OFFSET_MAP)
    if ssds_candidates.empty:
        raise RuntimeError("No SSD candidate points generated. Adjust spacing/kerb offsets or check roads.")

    # 3b) Enforce global minimum spacing between ANY two SSDs (e.g., 50 m)
    ssds_candidates = enforce_min_distance(ssds_candidates, MIN_SSD_SPACING)
    if ssds_candidates.empty:
        raise RuntimeError("All SSD candidates removed by min-distance filter. Lower MIN_SSD_SPACING or review inputs.")

    ssds_candidates.to_file(OUTPUT_DIR / "ssd_candidates.shp", driver="ESRI Shapefile")

    # 4) Optional: exclude buildings
    ssds_final = exclude_points_in_buildings(ssds_candidates, BUILDINGS_PATH, target_crs)

    # 4b) (Optional) enforce spacing again after exclusion (extra safety)
    ssds_final = enforce_min_distance(ssds_final, MIN_SSD_SPACING)

    ssds_final.to_file(OUTPUT_DIR / "ssds_final.shp", driver="ESRI Shapefile")

    # 5) CityCAT text export
    citycat_path = OUTPUT_DIR / CITYCAT_FILENAME
    write_citycat_inlets(
        ssds_final,
        citycat_path,
        header_triplet=CITYCAT_HEADER,
        coord_decimals=CITYCAT_COORD_DECIMALS,
        header_decimals=CITYCAT_HEADER_DECIMALS,
        active_flag=CITYCAT_ACTIVE_FLAG,
        id_start=0
    )
    print(f"[info] CityCAT inlet file written: {citycat_path}")

    # 6) Optional QA exports (buffered polygons + backbone lines)
    if EXPORT_BACKBONE:
        # buffered polygons for QA
        ensure_projected(roads, "Roads")
        road_poly_parts = []
        roads_norm = normalize_highway(roads)
        for klass, buf in BUFFERS_M.items():
            sel = roads_norm[roads_norm["__hw__"].str.contains(klass, na=False)]
            if not sel.empty:
                road_poly_parts.append(sel.buffer(buf, cap_style=2, join_style=2))
        if road_poly_parts:
            merged = pd.concat(road_poly_parts, ignore_index=True)
            dissolved = union_all_geoms(merged.geometry.values)
            if dissolved and not dissolved.is_empty:
                gser = gpd.GeoSeries(dissolved, crs=roads.crs).explode(index_parts=False).buffer(0)
                road_poly = gpd.GeoDataFrame(geometry=gser, crs=roads.crs)
                road_poly.to_file(OUTPUT_DIR / "roads_buffered_dissolved.shp", driver="ESRI Shapefile")
                # backbone lines from polygon rings
                lines = []
                for geom in road_poly.geometry:
                    if geom is None or geom.is_empty:
                        continue
                    lines.append(LineString(geom.exterior.coords))
                    for ring in getattr(geom, "interiors", []):
                        lines.append(LineString(ring.coords))
                if lines:
                    gpd.GeoDataFrame(geometry=lines, crs=road_poly.crs)\
                      .to_file(OUTPUT_DIR / "roads_backbone.shp", driver="ESRI Shapefile")

    # 7) Save run config
    (OUTPUT_DIR / "ssd_config.json").write_text(json.dumps({
        "target_crs": target_crs,
        "spacing_map": SPACING_MAP,
        "kerb_offset_map": KERB_OFFSET_MAP,
        "min_ssd_spacing": MIN_SSD_SPACING,
        "buffers_m_QA": BUFFERS_M,
        "citycat_header": CITYCAT_HEADER,
        "citycat_filename": CITYCAT_FILENAME,
        "study_area_path": str(STUDY_AREA_PATH),
        "buildings_path": str(BUILDINGS_PATH) if BUILDINGS_PATH else None,
        "roads_path": str(ROADS_PATH) if ROADS_PATH else None
    }, indent=2), encoding="utf-8")

    print("\n[done] Outputs:", OUTPUT_DIR.resolve())
    print(" - ssd_candidates.shp  (staggered kerb-side, min-distance enforced)")
    print(" - ssds_final.shp      (after building exclusion + min-distance)")
    print(f" - {CITYCAT_FILENAME}   (CityCAT inlet file)")
    if EXPORT_BACKBONE:
        print(" - roads_buffered_dissolved.shp (QA)")
        print(" - roads_backbone.shp          (QA)")

if __name__ == "__main__":
    sys.exit(main())
