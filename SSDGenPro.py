# SSDGen Pro — Synthetic Storm Drain Generator (kerb-side, staggered, QC)
# Author: Chris Iliadis (Newcastle University)

from __future__ import annotations

import sys
import json
import math
import time
import warnings
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import linemerge
from shapely.errors import ShapelyDeprecationWarning

warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)

# ---------- Version-safe union helper (Shapely 1.8/2.x) ----------
try:
    from shapely import union_all as _union_all  # Shapely >= 2
    def union_all_geoms(geoms):
        return _union_all(list(geoms))
except Exception:
    from shapely.ops import unary_union as _unary_union  # Shapely < 2
    def union_all_geoms(geoms):
        return _unary_union(list(geoms))


# ============================= CONFIG =============================

CONFIG = dict(
    # INPUTS
    STUDY_AREA_PATH=r"copy the shape file of the study area",
    BUILDINGS_PATH=r"copy the shape file of the buildings",
    ROADS_PATH=None,  # Optional: use your own roads; else download from OSM
    OSM_CACHE_PATH=r"copy the folder path for your outputs",  # Optional: path to cache downloaded OSM roads as GPKG (fast re-runs). "" disables caching.

    # CRS
    TARGET_CRS="add the projected coordinate system (e.g. EPSG:27700)",  # Set None to auto-UTM 
    VALIDATE_CRS=True,

    # ROAD FILTERS
    INCLUDE_HIGHWAYS=["motorway", "motorway_link", "trunk", "trunk_link",
                      "primary", "primary_link", "secondary", "secondary_link",
                      "tertiary", "tertiary_link"],
    INCLUDE_LOCALS=False,  # if True adds ["residential","unclassified","service","living_street"]

    # KERB OFFSETS (m)
    KERB_OFFSET_MAP={
        "motorway": 5.0, "motorway_link": 5.0,
        "trunk": 4.0, "trunk_link": 4.0,
        "primary": 3.5, "primary_link": 3.5,
        "secondary": 3.0, "secondary_link": 3.0,
        "tertiary": 3.0, "tertiary_link": 3.0,
        "_default": 3.0
    },
    # Optionally, derive kerb offset from OSM 'width' tag (metres) if present
    USE_OSM_WIDTH_WHEN_AVAILABLE=True,  # if True, kerb_offset = max(width/2 * 0.9, class_default)

    # SPACING (m) — this is the intended distance between consecutive inlets globally (staggered)
    SPACING_MAP={
        "motorway": 60, "motorway_link": 60,
        "trunk": 60, "trunk_link": 60,
        "primary": 50, "primary_link": 50,
        "secondary": 50, "secondary_link": 50,
        "tertiary": 50, "tertiary_link": 50,
        "_default": 50
    },
    # Stagger phase: for true staggering, use 0 for left, 0.5 for right (times spacing)
    STAGGER_PHASE_LEFT=0.0,
    STAGGER_PHASE_RIGHT=0.5,

    # MINIMUM DISTANCE ENFORCEMENT (m) — global
    MIN_SSD_SPACING=50.0,

    # QA BUFFERS (m) — optional visualisation
    QA_BUFFERS_M={
        "motorway": 25, "motorway_link": 25,
        "trunk": 20, "trunk_link": 20,
        "primary": 20, "primary_link": 20,
        "secondary": 10, "secondary_link": 10,
        "tertiary": 10, "tertiary_link": 10
    },

    # OUTPUTS
    OUTPUT_DIR=r"D:\PostDoc\Papers\Synthetic\Script\test\Results1",
    EXPORT_QA=True,         # buffered polygons + backbone
    EXPORT_GPKG=True,       # also export a GeoPackage alongside Shapefiles
    GPKG_NAME="ssdgen_outputs.gpkg",

    # CityCAT
    CITYCAT_HEADER=(0.30, 0.30, 0.50),
    CITYCAT_FILENAME="Inlets.txt",
    CITYCAT_COORD_DECIMALS=4,
    CITYCAT_HEADER_DECIMALS=2,
    CITYCAT_ACTIVE_FLAG=1,

    # LOGGING
    VERBOSE=True,
    RUN_TAG="",  # optional free text to store in QC report
)
# =================================================================


# ============================= CORE ==============================

@dataclass
class SSDGenConfig:
    study_area_path: str
    buildings_path: Optional[str]
    roads_path: Optional[str]
    osm_cache_path: Optional[str]
    target_crs: Optional[str]
    validate_crs: bool
    include_highways: List[str]
    include_locals: bool
    kerb_offset_map: Dict[str, float]
    use_osm_width_when_available: bool
    spacing_map: Dict[str, float]
    stagger_phase_left: float
    stagger_phase_right: float
    min_ssd_spacing: float
    qa_buffers_m: Dict[str, float]
    output_dir: str
    export_qa: bool
    export_gpkg: bool
    gpkg_name: str
    citycat_header: Tuple[float, float, float]
    citycat_filename: str
    citycat_coord_decimals: int
    citycat_header_decimals: int
    citycat_active_flag: int
    verbose: bool
    run_tag: str

    @staticmethod
    def from_dict(d: dict) -> "SSDGenConfig":
        # expand locals if requested
        include = list(d["INCLUDE_HIGHWAYS"])
        if d.get("INCLUDE_LOCALS", False):
            include.extend(["residential", "unclassified", "service", "living_street"])
        return SSDGenConfig(
            study_area_path=d["STUDY_AREA_PATH"],
            buildings_path=(d["BUILDINGS_PATH"] or None),
            roads_path=(d["ROADS_PATH"] or None),
            osm_cache_path=(d.get("OSM_CACHE_PATH") or None),
            target_crs=d.get("TARGET_CRS"),
            validate_crs=bool(d.get("VALIDATE_CRS", True)),
            include_highways=include,
            include_locals=bool(d.get("INCLUDE_LOCALS", False)),
            kerb_offset_map=d["KERB_OFFSET_MAP"],
            use_osm_width_when_available=bool(d.get("USE_OSM_WIDTH_WHEN_AVAILABLE", True)),
            spacing_map=d["SPACING_MAP"],
            stagger_phase_left=float(d.get("STAGGER_PHASE_LEFT", 0.0)),
            stagger_phase_right=float(d.get("STAGGER_PHASE_RIGHT", 0.5)),
            min_ssd_spacing=float(d["MIN_SSD_SPACING"]),
            qa_buffers_m=d["QA_BUFFERS_M"],
            output_dir=d["OUTPUT_DIR"],
            export_qa=bool(d["EXPORT_QA"]),
            export_gpkg=bool(d["EXPORT_GPKG"]),
            gpkg_name=d["GPKG_NAME"],
            citycat_header=tuple(d["CITYCAT_HEADER"]),
            citycat_filename=d["CITYCAT_FILENAME"],
            citycat_coord_decimals=int(d["CITYCAT_COORD_DECIMALS"]),
            citycat_header_decimals=int(d["CITYCAT_HEADER_DECIMALS"]),
            citycat_active_flag=int(d["CITYCAT_ACTIVE_FLAG"]),
            verbose=bool(d["VERBOSE"]),
            run_tag=str(d.get("RUN_TAG", "")),
        )


class SSDGen:
    def __init__(self, cfg: SSDGenConfig):
        self.cfg = cfg
        self.out_dir = Path(cfg.output_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self._log_header()

    # ------------------------- Logging -------------------------
    def log(self, msg: str):
        if self.cfg.verbose:
            print(msg)

    def _log_header(self):
        self.log("=== SSDGen Pro — Synthetic Storm Drain Generator ===")
        self.log(f"[config] Output dir: {self.out_dir}")
        if self.cfg.run_tag:
            self.log(f"[config] Run tag: {self.cfg.run_tag}")

    # -------------------- CRS helpers -------------------------
    @staticmethod
    def ensure_projected(gdf: gpd.GeoDataFrame, name: str):
        if gdf.crs is None:
            raise ValueError(f"{name} CRS is undefined.")
        if not gdf.crs.is_projected:
            raise ValueError(f"{name} CRS is geographic (degrees). Reproject to a projected CRS in metres.")

    @staticmethod
    def auto_target_crs(study_gdf: gpd.GeoDataFrame) -> str:
        try:
            from pyproj import CRS
            centroid = study_gdf.to_crs(4326).unary_union.centroid
            crs = CRS.estimate_utm_crs(latlon=(centroid.y, centroid.x))
            return crs.to_authority()[-1]
        except Exception:
            warnings.warn("Auto UTM failed; falling back to EPSG:3857.")
            return "EPSG:3857"

    # -------------------- Data loading ------------------------
    def load_study(self) -> gpd.GeoDataFrame:
        gdf = gpd.read_file(self.cfg.study_area_path)
        if gdf.crs is None:
            raise ValueError("Study area CRS is undefined.")
        target = self.cfg.target_crs or self.auto_target_crs(gdf)
        gdf = gdf.to_crs(target)
        if self.cfg.validate_crs:
            self.ensure_projected(gdf, "Study area")
        self.log(f"[CRS] Target CRS: {target}")
        return gdf

    def _download_roads_osm(self, poly_4326: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        try:
            import osmnx as ox
        except ImportError:
            raise RuntimeError("osmnx not installed. Install with: pip install osmnx")

        polygon = union_all_geoms(poly_4326.geometry.values)
        G = ox.graph_from_polygon(polygon, network_type="drive")
        edges = ox.graph_to_gdfs(G, nodes=False, edges=True)
        roads = edges[["geometry", "highway", "width"]].copy() if "width" in edges else edges[["geometry", "highway"]].copy()
        roads = gpd.overlay(roads, poly_4326, how="intersection")
        return roads

    def load_roads(self, study: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        if self.cfg.roads_path:
            roads = gpd.read_file(self.cfg.roads_path).to_crs(study.crs)
            roads = gpd.overlay(roads, study, how="intersection")
            self.log(f"[roads] Loaded user roads: {len(roads)}")
            return roads

        # Optional cache for OSM download
        if self.cfg.osm_cache_path:
            cache_path = Path(self.cfg.osm_cache_path)
            if cache_path.exists():
                try:
                    roads = gpd.read_file(cache_path)
                    roads = roads.to_crs(study.crs)
                    roads = gpd.overlay(roads, study, how="intersection")
                    self.log(f"[roads] Loaded OSM from cache: {cache_path} ({len(roads)})")
                    return roads
                except Exception:
                    self.log("[roads] Failed to read OSM cache; will download fresh.")

        # Fresh download
        roads_4326 = study.to_crs(4326)
        roads = self._download_roads_osm(roads_4326).to_crs(study.crs)
        self.log(f"[roads] Downloaded OSM roads: {len(roads)}")

        # Save cache if requested
        if self.cfg.osm_cache_path:
            cache_path = Path(self.cfg.osm_cache_path)
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            # write a single-layer GeoPackage
            roads.to_file(cache_path, layer="roads", driver="GPKG")
            self.log(f"[roads] Cached OSM roads to: {cache_path}")

        return roads

    # -------------------- Road processing ---------------------
    @staticmethod
    def normalize_highway(roads: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        roads = roads.copy()
        if "highway" in roads.columns:
            roads["__hw__"] = roads["highway"].apply(
                lambda v: ",".join(v).lower() if isinstance(v, (list, tuple, set)) else str(v).lower()
            )
        else:
            roads["__hw__"] = ""
        return roads

    def filter_by_highway(self, roads: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        roads = self.normalize_highway(roads)
        inc = self.cfg.include_highways
        mask = pd.Series(False, index=roads.index)
        for k in inc:
            mask = mask | roads["__hw__"].str.contains(k, na=False)
        filtered = roads[mask]
        self.log(f"[roads] Filtered by class → {len(filtered)} / {len(roads)} remain")
        return filtered

    # ---------------- Kerb + spacing helpers ------------------
    @staticmethod
    def _iter_lines(geom):
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
                yield from SSDGen._iter_lines(m)
            except Exception:
                return

    @staticmethod
    def _offset_lines(line: LineString, offset: float, side: str):
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

    def _kerb_offset_for_feature(self, row, klass: str) -> float:
        base = float(self.cfg.kerb_offset_map.get(klass, self.cfg.kerb_offset_map.get("_default", 3.0)))
        if self.cfg.use_osm_width_when_available and ("width" in row and pd.notna(row["width"])):
            try:
                w = float(str(row["width"]).split(";")[0].strip())  # handle "7; 8"
                derived = max(0.5 * w * 0.9, base)  # 90% half-width as kerb offset, min base
                return float(derived)
            except Exception:
                return base
        return base

    def points_on_kerbs_staggered(self, roads: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        self.ensure_projected(roads, "Roads")
        norm = self.normalize_highway(roads)
        pts, attrs = [], []
        covered = pd.Index([])

        # Explicit classes first
        for klass, spacing in self.cfg.spacing_map.items():
            if klass.startswith("_"):
                continue
            sel = norm[norm["__hw__"].str.contains(klass, na=False)]
            if sel.empty:
                continue
            covered = covered.union(sel.index)

            for _, row in sel.iterrows():
                for line in self._iter_lines(row.geometry):
                    kerb = self._kerb_offset_for_feature(row, klass)
                    left_lines = self._offset_lines(line, kerb, "left")
                    right_lines = self._offset_lines(line, kerb, "right")
                    if not left_lines or not right_lines:
                        continue
                    left, right = left_lines[0], right_lines[0]
                    L = min(left.length, right.length)
                    if L <= 0:
                        continue

                    # stagger with phases
                    step = spacing / 2.0
                    # left starts at phase_left * spacing ; right at phase_right * spacing
                    phase_L = (self.cfg.stagger_phase_left % 1.0) * spacing
                    phase_R = (self.cfg.stagger_phase_right % 1.0) * spacing

                    dists_L = np.arange(max(step, phase_L), L, step*2)  # every 2 steps
                    dists_R = np.arange(max(step, phase_R), L, step*2)

                    for d in dists_L:
                        if d < left.length:
                            pts.append(left.interpolate(d))
                            attrs.append({"class": klass, "side": "left",
                                          "spacing_m": spacing, "kerb_m": kerb})
                    for d in dists_R:
                        if d < right.length:
                            pts.append(right.interpolate(d))
                            attrs.append({"class": klass, "side": "right",
                                          "spacing_m": spacing, "kerb_m": kerb})

        # Default for remaining
        default_spacing = float(self.cfg.spacing_map.get("_default", 50))
        rest = norm[~norm.index.isin(covered)]
        for _, row in rest.iterrows():
            for line in self._iter_lines(row.geometry):
                kerb = self._kerb_offset_for_feature(row, "other")
                left_lines = self._offset_lines(line, kerb, "left")
                right_lines = self._offset_lines(line, kerb, "right")
                if not left_lines or not right_lines:
                    continue
                left, right = left_lines[0], right_lines[0]
                L = min(left.length, right.length)
                if L <= 0:
                    continue

                step = default_spacing / 2.0
                phase_L = (self.cfg.stagger_phase_left % 1.0) * default_spacing
                phase_R = (self.cfg.stagger_phase_right % 1.0) * default_spacing

                dists_L = np.arange(max(step, phase_L), L, step*2)
                dists_R = np.arange(max(step, phase_R), L, step*2)

                for d in dists_L:
                    if d < left.length:
                        pts.append(left.interpolate(d))
                        attrs.append({"class": "other", "side": "left",
                                      "spacing_m": default_spacing, "kerb_m": kerb})
                for d in dists_R:
                    if d < right.length:
                        pts.append(right.interpolate(d))
                        attrs.append({"class": "other", "side": "right",
                                      "spacing_m": default_spacing, "kerb_m": kerb})

        if not pts:
            return gpd.GeoDataFrame(columns=["class", "side", "spacing_m", "kerb_m", "geometry"], crs=roads.crs)
        return gpd.GeoDataFrame(attrs, geometry=pts, crs=roads.crs)

    # --------------- Min distance enforcement ----------------
    @staticmethod
    def _thinning_strtree(points: gpd.GeoDataFrame, min_dist: float) -> gpd.GeoDataFrame:
        pts = points.copy().reset_index(drop=True)
        pts["_x"] = pts.geometry.x
        pts["_y"] = pts.geometry.y
        pts = pts.sort_values(by=["_x", "_y"]).reset_index(drop=True)

        try:
            from shapely.strtree import STRtree
            geoms = list(pts.geometry.values)
            tree = STRtree(geoms)
            id_to_idx = {id(g): i for i, g in enumerate(geoms)}
            kept_mask = np.zeros(len(pts), dtype=bool)
            kept = []

            for i, g in enumerate(geoms):
                if kept_mask[i]:
                    continue
                neighbors = tree.query(g.buffer(min_dist))
                too_close = False
                for h in neighbors:
                    j = id_to_idx.get(id(h))
                    if j is not None and kept_mask[j]:
                        if g.distance(h) < min_dist:
                            too_close = True
                            break
                if not too_close:
                    kept_mask[i] = True
                    kept.append(i)
            out = pts.iloc[kept].drop(columns=["_x", "_y"])
            out.reset_index(drop=True, inplace=True)
            return out
        except Exception:
            # O(n^2) fallback
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

    # ---------------- Building exclusion ----------------------
    def exclude_points_in_buildings(self, points: gpd.GeoDataFrame, buildings_path: Optional[str], target_crs: str) -> gpd.GeoDataFrame:
        if not buildings_path or not Path(buildings_path).exists():
            return points
        bld = gpd.read_file(buildings_path).to_crs(target_crs)
        try:
            inside = gpd.sjoin(points, bld, predicate="within", how="left")
            keep = inside[inside.index_right.isna()].drop(columns=["index_right"])
            return gpd.GeoDataFrame(keep, geometry="geometry", crs=points.crs)
        except Exception:
            eps = 0.01
            pts_eps = points.copy()
            pts_eps["geometry"] = pts_eps.buffer(eps)
            dropped = gpd.overlay(pts_eps, bld, how="intersection")
            ids_to_drop = set(dropped.index)
            keep = points[~points.index.isin(ids_to_drop)]
            return gpd.GeoDataFrame(keep, geometry="geometry", crs=points.crs)

    # ---------------------- QA layers -------------------------
    def buffer_roads_for_qa(self, roads: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        roads = self.normalize_highway(roads)
        parts = []
        for klass, buf in self.cfg.qa_buffers_m.items():
            sel = roads[roads["__hw__"].str.contains(klass, na=False)]
            if not sel.empty:
                parts.append(sel.buffer(buf, cap_style=2, join_style=2))
        if not parts:
            return gpd.GeoDataFrame(geometry=[], crs=roads.crs)
        merged = pd.concat(parts, ignore_index=True)
        dissolved = union_all_geoms(merged.geometry.values)
        if dissolved is None or dissolved.is_empty:
            return gpd.GeoDataFrame(geometry=[], crs=roads.crs)
        gser = gpd.GeoSeries(dissolved, crs=roads.crs).explode(index_parts=False).buffer(0)
        return gpd.GeoDataFrame(geometry=gser, crs=roads.crs)

    # ---------------------- Exports ---------------------------
    def _sorted_points(self, gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        out = gdf.copy()
        out["_x"] = out.geometry.x
        out["_y"] = out.geometry.y
        out = out.sort_values(by=["_x", "_y"]).drop(columns=["_x", "_y"]).reset_index(drop=True)
        return out

    def export_layers(self, roads: gpd.GeoDataFrame, ssds_raw: gpd.GeoDataFrame, ssds_final: gpd.GeoDataFrame):
        # Shapefile
        ssds_raw.to_file(self.out_dir / "ssd_candidates.shp", driver="ESRI Shapefile")
        ssds_final.to_file(self.out_dir / "ssds_final.shp", driver="ESRI Shapefile")
        # GPKG
        if self.cfg.export_gpkg:
            gpkg = self.out_dir / self.cfg.gpkg_name
            if gpkg.exists():
                gpkg.unlink()
            ssds_raw.to_file(gpkg, layer="ssd_candidates", driver="GPKG")
            ssds_final.to_file(gpkg, layer="ssds_final", driver="GPKG")

        # QA
        if self.cfg.export_qa:
            road_poly = self.buffer_roads_for_qa(roads)
            if not road_poly.empty:
                road_poly.to_file(self.out_dir / "roads_buffered_dissolved.shp", driver="ESRI Shapefile")
                # backbone lines
                lines = []
                for geom in road_poly.geometry:
                    if geom is None or geom.is_empty:
                        continue
                    lines.append(LineString(geom.exterior.coords))
                    for ring in getattr(geom, "interiors", []):
                        lines.append(LineString(ring.coords))
                if lines:
                    gpd.GeoDataFrame(geometry=lines, crs=road_poly.crs)\
                        .to_file(self.out_dir / "roads_backbone.shp", driver="ESRI Shapefile")
            if self.cfg.export_gpkg:
                gpkg = self.out_dir / self.cfg.gpkg_name
                if not road_poly.empty:
                    road_poly.to_file(gpkg, layer="roads_buffered_dissolved", driver="GPKG")
                    lines_gdf = gpd.GeoDataFrame(geometry=lines, crs=road_poly.crs) if lines else gpd.GeoDataFrame(geometry=[], crs=road_poly.crs)
                    if not lines_gdf.empty:
                        lines_gdf.to_file(gpkg, layer="roads_backbone", driver="GPKG")

    def export_citycat(self, points: gpd.GeoDataFrame):
        path = self.out_dir / self.cfg.citycat_filename
        pts = self._sorted_points(points)
        fmt_h = f"{{:.{self.cfg.citycat_header_decimals}f}} {{:.{self.cfg.citycat_header_decimals}f}} {{:.{self.cfg.citycat_header_decimals}f}}\n"
        fmt_p = f"{{:d}} {{:.{self.cfg.citycat_coord_decimals}f}} {{:.{self.cfg.citycat_coord_decimals}f}} {{:d}}\n"
        with open(path, "w", encoding="utf-8") as f:
            f.write(fmt_h.format(*self.cfg.citycat_header))
            for i, row in pts.iterrows():
                x, y = row.geometry.x, row.geometry.y
                f.write(fmt_p.format(i, x, y, self.cfg.citycat_active_flag))
        self.log(f"[export] CityCAT file written: {path}")

    def export_config_and_qc(self, study: gpd.GeoDataFrame, roads: gpd.GeoDataFrame,
                             ssds_raw: gpd.GeoDataFrame, ssds_final: gpd.GeoDataFrame, duration_s: float):
        # Snapshot config
        snapshot = asdict(self.cfg)
        snapshot["CRS"] = study.crs.to_string()
        (self.out_dir / "ssd_config.json").write_text(json.dumps(snapshot, indent=2), encoding="utf-8")

        # QC report
        def total_length(gdf: gpd.GeoDataFrame) -> float:
            try:
                return float(gdf.length.sum())
            except Exception:
                return float("nan")

        qc = {
            "run_tag": self.cfg.run_tag,
            "crs": study.crs.to_string(),
            "roads_total": int(len(roads)),
            "roads_total_length_m": total_length(roads),
            "ssds_raw_count": int(len(ssds_raw)),
            "ssds_final_count": int(len(ssds_final)),
            "min_ssd_spacing_m": self.cfg.min_ssd_spacing,
            "export_gpkg": self.cfg.export_gpkg,
            "export_qa": self.cfg.export_qa,
            "duration_s": round(duration_s, 2)
        }
        pd.DataFrame([qc]).to_csv(self.out_dir / "ssd_qc_report.csv", index=False)
        self.log(f"[QC] Report written: {self.out_dir / 'ssd_qc_report.csv'}")

    # ------------------------- Run ----------------------------
    def run(self):
        t0 = time.time()

        # Study → CRS
        study = self.load_study()

        # Roads → filter classes
        roads = self.load_roads(study)
        roads = self.filter_by_highway(roads)
        if roads.empty:
            raise RuntimeError("No roads remain after highway filtering. Consider enabling INCLUDE_LOCALS or adjusting INCLUDE_HIGHWAYS.")

        self.log(f"[roads] Unique highway tags (sample): {sorted(pd.Series(roads['highway']).astype(str).str.lower().unique())[:30]}")

        # Generate kerb-side staggered SSDs
        ssds_raw = self.points_on_kerbs_staggered(roads)
        if ssds_raw.empty:
            raise RuntimeError("No SSD candidate points generated.")
        self.log(f"[place] SSD candidates: {len(ssds_raw)}")

        # Global min distance (pass 1)
        ssds_raw = self._thinning_strtree(ssds_raw, self.cfg.min_ssd_spacing)
        self.log(f"[thin] After pass 1 (global min dist): {len(ssds_raw)}")

        # Building exclusion
        ssds_final = self.exclude_points_in_buildings(ssds_raw, self.cfg.buildings_path, study.crs.to_string())
        self.log(f"[buildings] After exclusion: {len(ssds_final)}")

        # Global min distance (pass 2)
        ssds_final = self._thinning_strtree(ssds_final, self.cfg.min_ssd_spacing)
        self.log(f"[thin] After pass 2 (post-buildings): {len(ssds_final)}")

        if ssds_final.empty:
            raise RuntimeError("All SSDs removed by spacing/building filters. Review spacing, kerb offsets, or inputs.")

        # Exports
        self.export_layers(roads, ssds_raw, ssds_final)
        self.export_citycat(ssds_final)
        self.export_config_and_qc(study, roads, ssds_raw, ssds_final, time.time() - t0)

        self.log("\n[done] SSDGen Pro completed successfully.")
        self.log(f"[done] Outputs at: {self.out_dir.resolve()}")


# ============================ MAIN ==============================

def main():
    cfg = SSDGenConfig.from_dict(CONFIG)
    SSDGen(cfg).run()

if __name__ == "__main__":
    sys.exit(main())
