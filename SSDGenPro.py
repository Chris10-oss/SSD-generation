# SSDGen — Synthetic Storm Drain Generator
# Authors: Christos Iliadis
# Host: Newcastle University, University of Thessaly
# Final version, November 2026
# Run in Spyder, an IDE, or Python. Edit CONFIG below for your own datasets.
# Required packages: numpy, pandas, geopandas, shapely, pyproj and a vector
# file engine (pyogrio or fiona). OSM downloads also require osmnx.
# The synthetic demonstration also requires matplotlib.
# Geometry computations and published-study reproducibility should be checked
# with the same Shapely/GEOS versions and spatial datasets as the original run.

from __future__ import annotations

import sys
import json
import math
import time
import os
import tempfile
import warnings
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString, Point, box
from shapely.prepared import prep
from shapely.ops import linemerge
from shapely.errors import ShapelyDeprecationWarning

warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)

# ---------- Version-safe union helper (Shapely 1.8/2.x) ----------
try:
    from shapely import union_all as _union_all  # Shapely >= 2
    def union_all_geoms(geoms):
        return _union_all(list(geoms))
except ImportError:
    from shapely.ops import unary_union as _unary_union  # Shapely < 2
    def union_all_geoms(geoms):
        return _unary_union(list(geoms))


try:
    from shapely import line_interpolate_point as _interpolate_points
except ImportError:
    _interpolate_points = None

# Set True to generate local synthetic inputs, run the complete workflow, and
# create a demonstration plot. Set False to run your own CONFIG below.
SYNTHETIC_DEMO = False

# Performance changes preserve the original class matching, candidate order,
# offset selection and phase formula. The original phases 0.0 and 0.5 both
# start at spacing/2 under that formula. This version retains that rule.
# Corrected minimum-spacing enforcement can change results from Shapely 2 runs
# of the original script, whose STRtree query handling did not enforce spacing.

# ============================= Config =============================
# Edit below. No CLI needed.

CONFIG = dict(
    # INPUTS
    STUDY_AREA_PATH=r"copy the shape file of the study area",
    BUILDINGS_PATH=r"copy the shape file of the buildings",
    ROADS_PATH=None,  # Optional: use your own roads; else download from OSM
    OSM_CACHE_PATH="",  # Optional: path to cache downloaded OSM roads as GPKG (fast re-runs). "" disables caching. r"copy the folder path for your outputs"

    # CRS
    TARGET_CRS="add the projected coordinate system (e.g. EPSG:27700)",  # Set None to auto-UTM "add the projected coordinate system (e.g. EPSG:27700)"
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
    # Historical phase parameters (times spacing); see the placement note above CONFIG.
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
        self._validate_config()
        self.out_dir = Path(cfg.output_dir)
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self._log_header()

    def _validate_config(self):
        for label, values in (
            ("SPACING_MAP", self.cfg.spacing_map),
            ("KERB_OFFSET_MAP", self.cfg.kerb_offset_map),
            ("QA_BUFFERS_M", self.cfg.qa_buffers_m),
        ):
            for klass, value in values.items():
                if not math.isfinite(float(value)) or float(value) <= 0:
                    raise ValueError(f"{label}[{klass!r}] must be finite and positive.")
        if not math.isfinite(self.cfg.min_ssd_spacing) or self.cfg.min_ssd_spacing < 0:
            raise ValueError("MIN_SSD_SPACING must be finite and nonnegative.")
        phases = (self.cfg.stagger_phase_left, self.cfg.stagger_phase_right)
        if not all(math.isfinite(value) for value in phases):
            raise ValueError("Stagger phases must be finite.")
        if len(self.cfg.citycat_header) != 3 or not all(
            math.isfinite(float(value)) for value in self.cfg.citycat_header
        ):
            raise ValueError("CITYCAT_HEADER must contain three finite values.")
        if min(self.cfg.citycat_coord_decimals, self.cfg.citycat_header_decimals) < 0:
            raise ValueError("CityCAT decimal counts must be nonnegative.")

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
            raise ValueError(f"{name} CRS is geographic (degrees). Use a projected CRS in metres.")
        if any(not math.isclose(axis.unit_conversion_factor, 1.0)
               for axis in gdf.crs.axis_info[:2]):
            raise ValueError(f"{name} CRS must use metres for its horizontal axes.")

    @staticmethod
    def auto_target_crs(study_gdf: gpd.GeoDataFrame) -> str:
        crs = study_gdf.estimate_utm_crs()
        if crs is None:
            raise ValueError("Cannot estimate UTM CRS. Set TARGET_CRS explicitly.")
        return crs.to_string()

    # -------------------- Data loading ------------------------
    def load_study(self) -> gpd.GeoDataFrame:
        gdf = gpd.read_file(self.cfg.study_area_path)
        if gdf.crs is None:
            raise ValueError("Study area CRS is undefined.")
        if gdf.empty or gdf.geometry.is_empty.all():
            raise ValueError("Study area contains no usable polygons.")
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
        roads = gpd.overlay(roads, poly_4326[[poly_4326.geometry.name]], how="intersection")
        return roads

    def load_roads(self, study: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        if self.cfg.roads_path:
            roads = gpd.read_file(self.cfg.roads_path).to_crs(study.crs)
            roads = gpd.overlay(roads, study[[study.geometry.name]], how="intersection")
            self.log(f"[roads] Loaded user roads: {len(roads)}")
            return roads

        # Optional cache for OSM download
        if self.cfg.osm_cache_path:
            cache_path = Path(self.cfg.osm_cache_path)
            if cache_path.exists():
                try:
                    roads = gpd.read_file(cache_path)
                    roads = roads.to_crs(study.crs)
                    roads = gpd.overlay(roads, study[[study.geometry.name]], how="intersection")
                    self.log(f"[roads] Loaded OSM from cache: {cache_path} ({len(roads)})")
                    return roads
                except Exception:
                    self.log("[roads] Failed to read OSM cache; will download fresh.")

        # Fresh download
        roads_4326 = study.to_crs(4326)
        roads = self._download_roads_osm(roads_4326).to_crs(study.crs)
        # Use the same projected boundary on a fresh download and a cache reload.
        roads = gpd.overlay(roads, study[[study.geometry.name]], how="intersection")
        self.log(f"[roads] Downloaded OSM roads: {len(roads)}")

        # Save cache if requested
        if self.cfg.osm_cache_path:
            cache_path = Path(self.cfg.osm_cache_path)
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            # write a single-layer GeoPackage
            cache_roads = roads.copy()
            for column in cache_roads.select_dtypes(include="object").columns:
                if column == "width":
                    # Preserve the scalar-width fallback for multi-valued tags.
                    cache_roads[column] = cache_roads[column].map(
                        lambda value: str(value)
                        if isinstance(value, (list, tuple, set)) else value
                    )
                else:
                    cache_roads[column] = cache_roads[column].map(
                        lambda value: ",".join(map(str, value))
                        if isinstance(value, (list, tuple, set)) else value
                    )
            with tempfile.TemporaryDirectory(
                prefix=".ssd-cache-", dir=cache_path.parent
            ) as directory:
                temporary = Path(directory) / cache_path.name
                cache_roads.to_file(temporary, layer="roads", driver="GPKG", index=False)
                os.replace(temporary, cache_path)
            self.log(f"[roads] Cached OSM roads to: {cache_path}")

        return roads

    # -------------------- Road processing ---------------------
    @staticmethod
    def normalize_highway(roads: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        if "__hw__" in roads.columns:
            return roads
        roads = roads.copy()
        if "highway" in roads.columns:
            roads["__hw__"] = roads["highway"].apply(
                lambda v: ",".join(map(str, v)).lower() if isinstance(v, (list, tuple, set)) else str(v).lower()
            )
        else:
            roads["__hw__"] = ""
        return roads

    def filter_by_highway(self, roads: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        roads = self.normalize_highway(roads)
        inc = self.cfg.include_highways
        mask = pd.Series(False, index=roads.index)
        for k in inc:
            mask = mask | roads["__hw__"].str.contains(k, na=False, regex=False)
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
        base = float(self.cfg.kerb_offset_map.get(
            klass, self.cfg.kerb_offset_map.get("_default", 3.0)
        ))
        value = row.get("width")
        if not self.cfg.use_osm_width_when_available or value is None:
            return base
        # OSM simplified edges may hold several width tags. Preserve the old
        # scalar interpretation and use the class offset for multi-valued tags.
        try:
            width = float(str(value).split(";")[0].strip())
        except (TypeError, ValueError):
            return base
        if not math.isfinite(width):
            return base
        return max(0.5 * width * 0.9, base)

    def points_on_kerbs_staggered(self, roads: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        self.ensure_projected(roads, "Roads")
        norm = self.normalize_highway(roads)
        geometries = norm.geometry.to_numpy()
        widths = norm["width"].to_numpy() if "width" in norm else [None] * len(norm)
        covered = np.zeros(len(norm), dtype=bool)
        groups = []
        for klass, spacing in self.cfg.spacing_map.items():
            if klass.startswith("_"):
                continue
            selected = norm["__hw__"].str.contains(
                klass, na=False, regex=False
            ).to_numpy()
            covered |= selected
            groups.append((klass, spacing, np.flatnonzero(selected)))
        groups.append(("other", float(self.cfg.spacing_map.get("_default", 50)),
                       np.flatnonzero(~covered)))

        points, classes, sides, spacings, offsets = [], [], [], [], []
        for klass, spacing, indices in groups:
            # These station formulas deliberately match the supplied script.
            start_left = max(spacing / 2, (self.cfg.stagger_phase_left % 1) * spacing)
            start_right = max(spacing / 2, (self.cfg.stagger_phase_right % 1) * spacing)
            for index in indices:
                kerb = self._kerb_offset_for_feature({"width": widths[index]}, klass)
                for line in self._iter_lines(geometries[index]):
                    left_lines = self._offset_lines(line, kerb, "left")
                    right_lines = self._offset_lines(line, kerb, "right")
                    if not left_lines or not right_lines:
                        continue
                    left, right = left_lines[0], right_lines[0]
                    length = min(left.length, right.length)
                    for side, offset_line, start in (
                        ("left", left, start_left), ("right", right, start_right)
                    ):
                        distances = np.arange(start, length, spacing)
                        if not len(distances):
                            continue
                        if _interpolate_points is None:
                            generated = [offset_line.interpolate(d) for d in distances]
                        else:
                            generated = list(_interpolate_points(offset_line, distances))
                        count = len(generated)
                        points.extend(generated)
                        classes.extend([klass] * count)
                        sides.extend([side] * count)
                        spacings.extend([spacing] * count)
                        offsets.extend([kerb] * count)
        return gpd.GeoDataFrame(
            {"class": classes, "side": sides, "spacing_m": spacings, "kerb_m": offsets},
            geometry=points, crs=roads.crs,
        )

    # --------------- Min distance enforcement ----------------
    @staticmethod
    def _thinning_strtree(points: gpd.GeoDataFrame, min_dist: float) -> gpd.GeoDataFrame:
        """Greedy XY thinning using a spatial hash; retains the public method name.

        Candidates are considered in the original lexicographic XY order.
        Only retained points in neighbouring cells need a distance check.
        Equality with min_dist is allowed, matching the original strict rule.
        """
        if not math.isfinite(min_dist) or min_dist < 0:
            raise ValueError("Minimum spacing must be finite and nonnegative.")
        if points.empty:
            return points.copy().reset_index(drop=True)
        if points.geometry.isna().any() or not points.geom_type.eq("Point").all():
            raise ValueError("Spacing enforcement requires nonempty Point geometries.")
        pts = points.copy().reset_index(drop=True)
        pts["_x"] = pts.geometry.x
        pts["_y"] = pts.geometry.y
        pts = pts.sort_values(["_x", "_y"]).reset_index(drop=True)
        coordinates = pts[["_x", "_y"]].to_numpy()
        if not np.isfinite(coordinates).all():
            raise ValueError("Inlet coordinates must be finite.")
        if min_dist == 0:
            return pts.drop(columns=["_x", "_y"])

        cells = {}
        kept = []
        for index, (x, y) in enumerate(coordinates):
            cell_x, cell_y = math.floor(x / min_dist), math.floor(y / min_dist)
            neighbours = (
                neighbour
                for dx in (-1, 0, 1)
                for dy in (-1, 0, 1)
                for neighbour in cells.get((cell_x + dx, cell_y + dy), ())
            )
            if any(math.hypot(x - hx, y - hy) < min_dist for hx, hy in neighbours):
                continue
            cells.setdefault((cell_x, cell_y), []).append((x, y))
            kept.append(index)
        return pts.iloc[kept].drop(columns=["_x", "_y"]).reset_index(drop=True)

    # ---------------- Building exclusion ----------------------
    def exclude_points_in_buildings(self, points: gpd.GeoDataFrame,
                                    buildings_path: Optional[str],
                                    target_crs: str) -> gpd.GeoDataFrame:
        if points.empty or not buildings_path:
            return points.copy()
        if not Path(buildings_path).exists():
            warnings.warn(f"Buildings file missing; exclusion skipped: {buildings_path}")
            return points.copy()
        # Geometry-only reads avoid allocating building attribute tables.
        buildings = gpd.read_file(buildings_path, columns=[])
        if buildings.crs is None:
            raise ValueError("Buildings CRS is undefined.")
        buildings = buildings.to_crs(target_crs)
        buildings = buildings.loc[
            buildings.geometry.notna() & ~buildings.geometry.is_empty
        ]
        if buildings.empty:
            return points.copy()
        # Envelope preselection reduces the spatial index to nearby footprints.
        min_x, min_y, max_x, max_y = points.total_bounds
        # Expand by one floating-point step so a single point or collinear
        # points produce a valid envelope on older GEOS versions.
        envelope = box(
            np.nextafter(min_x, -np.inf), np.nextafter(min_y, -np.inf),
            np.nextafter(max_x, np.inf), np.nextafter(max_y, np.inf),
        )
        buildings = buildings.loc[buildings.geometry.intersects(envelope)]
        if buildings.empty:
            return points.copy()
        try:
            spatial_index = buildings.sindex
        except ImportError:
            footprints = [prep(geometry) for geometry in buildings.geometry]
            inside = np.fromiter(
                (any(footprint.contains(point) for footprint in footprints)
                 for point in points.geometry),
                dtype=bool, count=len(points),
            )
            return points.iloc[np.flatnonzero(~inside)].copy().reset_index(drop=True)

        # query_bulk supports GeoPandas 0.10-0.13; query supports newer releases.
        query = getattr(spatial_index, "query_bulk", spatial_index.query)
        keep = np.ones(len(points), dtype=bool)
        batch_size = 50_000
        for start in range(0, len(points), batch_size):
            chunk = points.geometry.iloc[start:start + batch_size]
            pairs = query(chunk, predicate="within")
            keep[start + np.unique(pairs[0])] = False
        return points.iloc[np.flatnonzero(keep)].copy().reset_index(drop=True)

    # ---------------------- QA layers -------------------------
    def buffer_roads_for_qa(self, roads: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        roads = self.normalize_highway(roads)
        parts = []
        seen = {}
        for klass, buf in self.cfg.qa_buffers_m.items():
            selected = roads["__hw__"].str.contains(
                klass, na=False, regex=False
            ).to_numpy()
            previous = seen.setdefault(float(buf), np.zeros(len(roads), dtype=bool))
            new_indices = np.flatnonzero(selected & ~previous)
            previous |= selected
            if len(new_indices):
                parts.append(roads.geometry.iloc[new_indices].buffer(
                    float(buf), cap_style=2, join_style=2
                ))
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

    def export_layers(self, roads: gpd.GeoDataFrame,
                      ssds_raw: gpd.GeoDataFrame, ssds_final: gpd.GeoDataFrame):
        layers = [
            ("ssd_candidates", ssds_raw),
            ("ssds_final", ssds_final),
        ]
        if self.cfg.export_qa:
            road_poly = self.buffer_roads_for_qa(roads)
            if not road_poly.empty:
                layers.append(("roads_buffered_dissolved", road_poly))
                lines = []
                for geometry in road_poly.geometry:
                    if geometry is None or geometry.is_empty:
                        continue
                    lines.append(LineString(geometry.exterior.coords))
                    lines.extend(LineString(ring.coords) for ring in geometry.interiors)
                if lines:
                    layers.append(("roads_backbone", gpd.GeoDataFrame(
                        geometry=lines, crs=road_poly.crs
                    )))

        for name, frame in layers:
            frame.to_file(self.out_dir / f"{name}.shp",
                          driver="ESRI Shapefile", index=False)
        if self.cfg.export_gpkg:
            path = self.out_dir / self.cfg.gpkg_name
            # Finish every layer before replacing an existing GeoPackage.
            with tempfile.TemporaryDirectory(prefix=".ssd-gpkg-", dir=path.parent) as directory:
                temporary = Path(directory) / path.name
                for name, frame in layers:
                    frame.to_file(temporary, layer=name, driver="GPKG", index=False)
                os.replace(temporary, path)

    def export_citycat(self, points: gpd.GeoDataFrame):
        path = self.out_dir / self.cfg.citycat_filename
        pts = self._sorted_points(points)
        header = " ".join(f"{float(value):.{self.cfg.citycat_header_decimals}f}"
                          for value in self.cfg.citycat_header) + "\n"
        coordinates = zip(pts.geometry.x.to_numpy(), pts.geometry.y.to_numpy())
        decimals = self.cfg.citycat_coord_decimals
        # Publish only a completely written CityCAT file.
        handle, temporary = tempfile.mkstemp(prefix=".inlets-", dir=path.parent)
        try:
            with os.fdopen(handle, "w", encoding="utf-8") as stream:
                stream.write(header)
                stream.writelines(
                    f"{index:d} {x:.{decimals}f} {y:.{decimals}f} "
                    f"{self.cfg.citycat_active_flag:d}\n"
                    for index, (x, y) in enumerate(coordinates)
                )
            os.replace(temporary, path)
        finally:
            Path(temporary).unlink(missing_ok=True)
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
        t0 = time.perf_counter()

        # Study → CRS
        study = self.load_study()

        # Roads → filter classes
        roads = self.load_roads(study)
        roads = self.filter_by_highway(roads)
        if roads.empty:
            raise RuntimeError("No roads remain after highway filtering. Consider enabling INCLUDE_LOCALS or adjusting INCLUDE_HIGHWAYS.")

        self.log(f"[roads] Unique highway tags (sample): {sorted(pd.Series(roads['__hw__']).astype(str).str.lower().unique())[:30]}")

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

        # Removing points cannot introduce a spacing violation, so a second
        # greedy thinning pass is unnecessary.
        ssds_final = self._sorted_points(ssds_final)

        if ssds_final.empty:
            raise RuntimeError("All SSDs removed by spacing/building filters. Review spacing, kerb offsets, or inputs.")

        # Exports
        self.export_layers(roads, ssds_raw, ssds_final)
        self.export_citycat(ssds_final)
        self.export_config_and_qc(study, roads, ssds_raw, ssds_final, time.perf_counter() - t0)

        self.log("\n[done] SSDGen Pro completed successfully.")
        self.log(f"[done] Outputs at: {self.out_dir.resolve()}")
        return ssds_raw, ssds_final


def run_synthetic_demo(base_dir: Optional[Path] = None):
    """Create deterministic synthetic inputs and exercise the normal run entry point."""
    from copy import deepcopy
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch, Rectangle
    from matplotlib.ticker import ScalarFormatter

    if base_dir is None:
        base_dir = Path(__file__).resolve().parent / "SSDGen_Synthetic_Demo"
    base_dir = Path(base_dir)
    inputs = base_dir / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    crs = "EPSG:2100"
    origin_x, origin_y = 410000.0, 4498000.0
    area = box(origin_x, origin_y, origin_x + 1200, origin_y + 900)
    study = gpd.GeoDataFrame({"name": ["Synthetic study area"]},
                             geometry=[area], crs=crs)

    roads, highway, width = [], [], []
    # A connected street grid with varying widths and road classes.
    for i, y in enumerate((140, 350, 560, 770)):
        roads.append(LineString([(origin_x - 40, origin_y + y),
                                 (origin_x + 1240, origin_y + y)]))
        highway.append("primary" if i == 1 else "secondary")
        width.append("10" if i == 1 else "7")
    for x in (150, 390, 630, 870, 1110):
        roads.append(LineString([(origin_x + x, origin_y - 40),
                                 (origin_x + x, origin_y + 940)]))
        highway.append("tertiary")
        width.append("6")
    roads.append(LineString([
        (origin_x - 30, origin_y + 15), (origin_x + 300, origin_y + 150),
        (origin_x + 680, origin_y + 330), (origin_x + 1240, origin_y + 600),
    ]))
    highway.append("primary_link")
    width.append("8; 9")
    road_frame = gpd.GeoDataFrame({"highway": highway, "width": width},
                                  geometry=roads, crs=crs)

    footprints = []
    for x in (150, 390, 630, 870):
        for y in (140, 350, 560):
            for offset_x in (28, 108):
                footprints.append(box(
                    origin_x + x + offset_x, origin_y + y + 30,
                    origin_x + x + offset_x + 66, origin_y + y + 138,
                ))
    # Two deliberate road/footprint conflicts exercise building exclusion.
    footprints.extend([
        box(origin_x + 455, origin_y + 344, origin_x + 495, origin_y + 360),
        box(origin_x + 1080, origin_y + 740, origin_x + 1140, origin_y + 805),
    ])
    buildings = gpd.GeoDataFrame({"building_id": np.arange(1, len(footprints) + 1)},
                                 geometry=footprints, crs=crs)
    study_path = inputs / "study_area.gpkg"
    roads_path = inputs / "roads.gpkg"
    buildings_path = inputs / "buildings.gpkg"
    study.to_file(study_path, layer="study_area", driver="GPKG", index=False)
    road_frame.to_file(roads_path, layer="roads", driver="GPKG", index=False)
    buildings.to_file(buildings_path, layer="buildings", driver="GPKG", index=False)

    settings = deepcopy(CONFIG)
    settings.update(STUDY_AREA_PATH=str(study_path), ROADS_PATH=str(roads_path),
                    BUILDINGS_PATH=str(buildings_path), OSM_CACHE_PATH=None,
                    OUTPUT_DIR=str(base_dir / "results"), TARGET_CRS=crs,
                    EXPORT_QA=True, EXPORT_GPKG=True,
                    RUN_TAG="Reproducible synthetic demonstration")
    generator = SSDGen(SSDGenConfig.from_dict(settings))
    raw, final = generator.run()
    clipped_roads = generator.load_roads(generator.load_study())
    candidates = generator.points_on_kerbs_staggered(
        generator.filter_by_highway(clipped_roads)
    )
    coordinates = np.column_stack([final.geometry.x, final.geometry.y])
    # Independent, exhaustive validation is inexpensive for this small demo.
    distances = np.hypot(coordinates[:, 0, None] - coordinates[None, :, 0],
                         coordinates[:, 1, None] - coordinates[None, :, 1])
    np.fill_diagonal(distances, np.inf)
    minimum_distance = float(distances.min())
    assert minimum_distance >= settings["MIN_SSD_SPACING"] - 1e-9
    assert not any(point.within(building) for point in final.geometry
                   for building in buildings.geometry)
    citycat_lines = (generator.out_dir / settings["CITYCAT_FILENAME"]).read_text().splitlines()
    assert len(citycat_lines) == len(final) + 1
    assert len(gpd.read_file(generator.out_dir / settings["GPKG_NAME"],
                             layer="ssds_final")) == len(final)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.titlesize": 12, "axes.labelsize": 10})
    fig = plt.figure(figsize=(14, 8.5), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=(1.55, 1))
    overview = fig.add_subplot(grid[:, 0])
    zoom = fig.add_subplot(grid[0, 1])
    checks = fig.add_subplot(grid[1, 1])
    for ax in (overview, zoom):
        buildings.plot(ax=ax, color="#d7d9dd", edgecolor="#abb1ba", linewidth=0.4)
        road_frame.plot(ax=ax, color="#747f8f", linewidth=1.2)
        ax.set_aspect("equal")
        ax.set_xlabel("Easting (m)")
        ax.set_ylabel("Northing (m)")
        ax.tick_params(axis="x", rotation=30)
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(False)
        ax.xaxis.set_major_formatter(formatter)
        formatter_y = ScalarFormatter(useOffset=False)
        formatter_y.set_scientific(False)
        ax.yaxis.set_major_formatter(formatter_y)
        ax.spines[["top", "right"]].set_visible(False)
    candidates.plot(ax=overview, color="#dfe7f0", markersize=8, zorder=3)
    final.plot(ax=overview, color="#086d8c", markersize=15, zorder=4)
    study.boundary.plot(ax=overview, color="#384453", linewidth=0.8)
    overview.set_xlim(origin_x - 30, origin_x + 1230)
    overview.set_ylim(origin_y - 30, origin_y + 930)
    overview.set_title("Generated synthetic inlets", loc="left", fontweight="bold")
    zoom_bounds = (origin_x + 335, origin_y + 285, 235, 160)
    overview.add_patch(Rectangle(
        zoom_bounds[:2], zoom_bounds[2], zoom_bounds[3], fill=False,
        edgecolor="#c66c12", linewidth=1.5, zorder=5,
    ))
    removed = raw.loc[~raw.geometry.isin(final.geometry)]
    if not removed.empty:
        removed.plot(ax=overview, color="#d4483c", marker="x",
                     markersize=28, linewidth=1.2, zorder=5)
    raw.plot(ax=zoom, color="#b5c6d8", markersize=35, zorder=3)
    final.plot(ax=zoom, color="#086d8c", markersize=40, zorder=4)
    if not removed.empty:
        removed.plot(ax=zoom, color="#d4483c", marker="x",
                     markersize=75, linewidth=1.8, zorder=5)
    zoom.set_xlim(zoom_bounds[0], zoom_bounds[0] + zoom_bounds[2])
    zoom.set_ylim(zoom_bounds[1], zoom_bounds[1] + zoom_bounds[3])
    zoom.set_title("Building exclusion at a road intersection", loc="left")
    legend = [
        Line2D([], [], color="#747f8f", label="Synthetic roads"),
        Patch(facecolor="#d7d9dd", edgecolor="#abb1ba", label="Synthetic buildings"),
        Line2D([], [], linestyle="", marker="o", color="#086d8c",
               label="Final synthetic inlets"),
        Line2D([], [], linestyle="", marker="x", color="#d4483c",
               label="Removed inside buildings"),
    ]
    overview.legend(handles=legend, loc="upper left", frameon=True, fontsize=9)
    labels = ["Generated", "After spacing", "Final"]
    counts = [len(candidates), len(raw), len(final)]
    bars = checks.bar(labels, counts, color=["#b5c6d8", "#4895ad", "#086d8c"], width=0.55)
    checks.bar_label(bars, padding=4)
    checks.set_ylim(0, max(counts) * 1.23)
    checks.set_ylabel("Number of inlets")
    checks.set_title("Workflow checks", loc="left")
    checks.spines[["top", "right"]].set_visible(False)
    checks.text(
        0, -0.29,
        f"Minimum final spacing   {minimum_distance:.2f} m\n"
        f"Inlets inside buildings   0\n"
        f"CityCAT and GeoPackage counts   verified",
        transform=checks.transAxes, fontsize=10, linespacing=1.7,
    )
    fig.suptitle("SSDGen synthetic demonstration", fontsize=17, fontweight="bold")
    plot_path = base_dir / "SSDGen_Synthetic_Demo.png"
    fig.savefig(plot_path, dpi=600, facecolor="white")
    plt.close(fig)
    generator.log(f"[demo] Synthetic illustration: {plot_path.resolve()}")
    generator.log(f"[demo] Minimum final spacing: {minimum_distance:.6f} m")
    return {"plot": plot_path, "candidates": len(candidates),
            "after_spacing": len(raw), "final": len(final),
            "minimum_spacing_m": minimum_distance}


# ============================ MAIN ==============================

def main():
    if SYNTHETIC_DEMO:
        run_synthetic_demo()
    else:
        cfg = SSDGenConfig.from_dict(CONFIG)
        SSDGen(cfg).run()

if __name__ == "__main__":
    sys.exit(main())
