##SSDGen Pro — Synthetic Storm Drain Generator

SSDGen Pro is an automated tool for generating synthetic storm drains (SSDs) in data-scarce cities. It places inlets along urban road kerb lines using spatial analysis and exports GIS layers (Shapefile/GeoPackage) plus a CityCAT-ready inlet file. The method assumes infinite downstream capacity and is designed to couple surface and sub-surface processes in hydrodynamic models.

License: Apache 2.0

Repository: https://github.com/Chris10-oss/SSD-generation

Latest release: v1.0 (2025)

Features

Kerb-side inlet placement (left/right offsets), staggered pattern (~50 m)

Global min-distance filter (two-pass) to avoid clustering

Optional building exclusion (removes SSDs inside footprints)

Road data from OSM (auto-download) or user-provided roads

Auto-UTM CRS or user-defined (e.g., EPSG:2100 / EPSG:27700)

Outputs: Shapefiles + optional GeoPackage, CityCAT inlet TXT, QC report, config snapshot

Optional OSM caching to speed up re-runs

Quick start

Clone and enter the repo

git clone https://github.com/Chris10-oss/SSD-generation.git
cd SSD-generation


Create a Python environment (recommended)

# using conda (recommended on Windows)
conda create -n ssdgen python=3.11 -y
conda activate ssdgen


Install dependencies

# core geospatial stack
conda install -c conda-forge geopandas shapely pyproj rtree gdal fiona -y
# network + analysis
conda install -c conda-forge osmnx pandas numpy -y


Edit configuration in the script

Open ssdgen_pro.py and update the CONFIG block at the top:

CONFIG = dict(
    STUDY_AREA_PATH=r"PATH/TO/study_area.shp",   # polygon
    BUILDINGS_PATH=r"PATH/TO/buildings.shp",     # optional
    ROADS_PATH=None,                             # optional (use your own road lines)
    OSM_CACHE_PATH=r"",                          # optional GPKG path to cache OSM roads

    TARGET_CRS="EPSG:2100",                      # or "EPSG:27700" or None for auto-UTM

    # spacing and offsets (edit if needed)
    SPACING_MAP={"primary":50,"secondary":50,"tertiary":50,"_default":50},
    KERB_OFFSET_MAP={"primary":3.5,"secondary":3.0,"tertiary":3.0,"_default":3.0},

    MIN_SSD_SPACING=50.0,                        # global minimum distance between inlets (m)

    OUTPUT_DIR=r"outputs",                       # where results are written
    EXPORT_GPKG=True,
    EXPORT_QA=True,

    # CityCAT export
    CITYCAT_HEADER=(0.30, 0.30, 0.50),           # size_x size_y friction
    CITYCAT_FILENAME="citycat_inlets.txt",
)

Run

python ssdgen_pro.py

Inputs

Study area (required): polygon boundary (SHP/GPKG/GeoJSON). Must have a CRS.

Buildings (optional): polygons to exclude SSDs inside buildings.

Roads (optional): your own line dataset. If not provided, roads are downloaded from OSM within the study area.

CRS: set TARGET_CRS="EPSG:xxxx" or leave None for auto-UTM.

Tip: For the UK use EPSG:27700; for Greece use EPSG:2100.

Outputs (in OUTPUT_DIR)

ssd_candidates.shp — pre-filter SSD points

ssds_final.shp — final SSD points after spacing + building exclusion

ssdgen_outputs.gpkg (if EXPORT_GPKG=True) — same layers in GPKG

roads_buffered_dissolved.shp, roads_backbone.shp (if EXPORT_QA=True)

citycat_inlets.txt — CityCAT-ready inlet file (see format below)

ssd_config.json — run configuration snapshot

ssd_qc_report.csv — summary (counts, spacing, timing)

CityCAT file format
0.30 0.30 0.50          # header: size_x size_y friction (from CONFIG.CITYCAT_HEADER)
0   425194.7739 564040.3189  1
1   425498.7496 564098.4754  1
...


Columns after header: id x y active_flag

active_flag = 1 for active inlets

Coordinates are written in the TARGET_CRS units (metres)

How it works (in brief)

Road acquisition: Use OSM (via OSMnx) or user roads; filter by OSM highway classes.

Kerb generation: Offset centreline left/right by class-based kerb distances (optionally using OSM width tags where available).

Staggered placement: Alternate left/right with ~50 m design spacing.

Global spacing filter: Enforce a minimum distance (default 50 m) across all inlets (two-pass: before and after building exclusion).

Building exclusion: Remove SSDs inside building footprints.

Export: GIS layers, CityCAT file, config snapshot, QC report.

Configuration tips

Spacing per class: edit SPACING_MAP (e.g., higher spacing for motorways).

Kerb offsets: edit KERB_OFFSET_MAP or enable width-based refinement by keeping USE_OSM_WIDTH_WHEN_AVAILABLE=True (in script).

Include locals: set INCLUDE_LOCALS=True to include residential/service/unclassified classes.

Cache OSM: set OSM_CACHE_PATH="path/to/roads_cache.gpkg" to avoid re-downloading.

Troubleshooting

“Study area CRS is undefined”
→ Define the CRS of your study polygon (e.g., in QGIS) before running.

No roads found / “Filtered by class → 0 remain”
→ Enable INCLUDE_LOCALS=True or broaden INCLUDE_HIGHWAYS. Ensure the study area is correct.

ImportError: osmnx not installed
→ Install with conda-forge: conda install -c conda-forge osmnx.

GDAL/Fiona install issues (Windows)
→ Prefer conda-forge and install geopandas shapely pyproj rtree gdal fiona via conda, not pip.

PROJ / CRS errors
→ Ensure pyproj and proj data are installed (conda install -c conda-forge pyproj proj).

Very close SSDs remain
→ Increase MIN_SSD_SPACING (e.g., 60–75 m). Remember spacing is enforced globally after staggering.

Reproducibility

Each run writes:

ssd_config.json — exact configuration used

ssd_qc_report.csv — counts and summary metrics

Commit these with your results to track provenance.

Contributing

Issues and pull requests are welcome. Please:

Open an issue describing the bug/feature first.

Use feature branches and include small, focused PRs.

Add a short description to ssd_qc_report.csv if you change core logic.

License

Apache License 2.0 (January 2004). See LICENSE file for details.
