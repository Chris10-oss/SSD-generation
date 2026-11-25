# SSDGen Pro — Synthetic Storm Drain Generator

SSDGen Pro generates synthetic storm drain inlet locations for urban flood modelling in areas with limited drainage records. It derives left and right kerb lines from road centrelines, places candidate inlets using class-specific offsets and spacing, and applies a minimum distance across the study area. An optional building footprint layer allows the tool to remove candidates that fall inside buildings.

The script exports inlet locations as Shapefiles, an optional GeoPackage, and a CityCAT-ready text file. In the associated study, CityCAT represents the inlets as local abstractions from the 2D surface with unlimited downstream capacity. Captured water leaves the surface model, and the simulations omit pipe routing and surcharge return.

**Repository** · [Chris10-oss/SSD-generation](https://github.com/Chris10-oss/SSD-generation)  
**License** · [Apache License 2.0](LICENSE)

## Quick start

Clone the repository:

```bash
git clone https://github.com/Chris10-oss/SSD-generation.git
cd SSD-generation
```

Create a Python environment. On Windows, you can use Miniforge or Anaconda:

```bash
conda create -n ssdgen -c conda-forge python=3.11 geopandas shapely pyproj pyogrio matplotlib -y
conda activate ssdgen
```

If you want to download roads from OpenStreetMap, install OSMnx:

```bash
conda install -n ssdgen -c conda-forge osmnx -y
```

Open `SSDGenPro.py` and edit the entries below in the existing `CONFIG` block. Keep the other entries in that block so the script has all the settings it expects.

```python
STUDY_AREA_PATH=r"C:\data\study_area.shp",
BUILDINGS_PATH=r"C:\data\buildings.shp",  # Set to None to skip building exclusion
ROADS_PATH=None,                        # Set a local path or use OSM with None
OSM_CACHE_PATH=r"cache\osm_roads.gpkg",  # Use "" to disable caching

TARGET_CRS="EPSG:2100",                 # For this Greek example
MIN_SSD_SPACING=50.0,                   # Metres

OUTPUT_DIR=r"outputs",
EXPORT_GPKG=True,
EXPORT_QA=True,

CITYCAT_HEADER=(0.30, 0.30, 0.50),
CITYCAT_FILENAME="Inlets.txt",
```

Run the script in Spyder using the `ssdgen` environment, or run it from a terminal:

```bash
python SSDGenPro.py
```

## Inputs and configuration

| Input or setting | Purpose |
|---|---|
| `STUDY_AREA_PATH` | Required polygon boundary with a defined coordinate reference system (CRS). Common formats include Shapefile, GeoPackage and GeoJSON. |
| `BUILDINGS_PATH` | Optional polygon layer. Points strictly inside building footprints are removed. |
| `ROADS_PATH` | Optional local road line layer. Include a `highway` attribute with classes such as `primary`, `secondary` or `tertiary` so the configured filters can select roads. |
| `TARGET_CRS` | Projected CRS in metres. Use a CRS appropriate for the study area, such as `EPSG:27700` for Great Britain or `EPSG:2100` for Greece. Set `None` for automatic UTM selection. |
| `SPACING_MAP` | Candidate spacing along each kerb, in metres by road class. |
| `KERB_OFFSET_MAP` | Lateral offsets from the road centreline, in metres by road class. Set `USE_OSM_WIDTH_WHEN_AVAILABLE=True` to use suitable road width tags where available. |
| `MIN_SSD_SPACING` | Minimum straight-line distance between any two retained inlets, including inlets on opposite kerbs or different roads. |
| `INCLUDE_LOCALS` | Set `True` to include residential, service, unclassified and living streets. |
| `OSM_CACHE_PATH` | Path to an optional GeoPackage cache. Use a separate cache for each study area. |

## How it works

1. **Select roads.** Read a local road layer or download roads from OpenStreetMap using OSMnx. Clip the roads to the study area and filter by highway class.
2. **Generate kerb lines.** Offset each road centreline to the left and right using the configured kerb distances.
3. **Place candidate inlets.** Sample points along both kerbs using the road class spacing and phase parameters. With the default phase settings, sampling starts at half the class spacing on both kerbs.
4. **Enforce global spacing.** Sort candidates by easting and northing, then retain points that meet the minimum distance from previously retained points. This step also controls clustering near junctions.
5. **Exclude buildings.** Remove retained points that lie inside the supplied building footprints. Since removing points cannot introduce a spacing conflict, the final workflow uses one spacing pass.
6. **Export results.** Write the inlet layers, CityCAT input, run configuration and quality-control report.

## Synthetic demonstration

Set `SYNTHETIC_DEMO = True` near the top of the script and run it. The demonstration generates its own study boundary, road network and building footprints. It runs locally without downloading OpenStreetMap data and writes its inputs, outputs and a 600 dpi plot to `SSDGen_Synthetic_Demo/` beside the script.

The demonstration generates **476 candidates**, retains **160** after spacing enforcement, and removes **7** inside building footprints. The **153 final inlets** have a minimum separation of **50 m**, and none lies inside a building footprint.

![SSDGen synthetic demonstration showing road centrelines, buildings, retained inlets and an example of building exclusion](SSDGen_Synthetic_Demo.png)

## Output files

Files are written to `OUTPUT_DIR`.

| File | Contents |
|---|---|
| `ssd_candidates.shp` | Candidate inlet locations after the global minimum spacing filter and before building exclusion. |
| `ssds_final.shp` | Final inlet locations after building exclusion. |
| `ssdgen_outputs.gpkg` | Optional GeoPackage containing both inlet layers and any enabled road layers for visual checks. |
| `roads_buffered_dissolved.shp` | Optional dissolved road buffers for visual checks. |
| `roads_backbone.shp` | Optional linework derived from the dissolved road buffer boundaries. |
| `Inlets.txt` | CityCAT inlet input file. Change `CITYCAT_FILENAME` to use a different name. |
| `ssd_config.json` | Saved run settings and the CRS used. |
| `ssd_qc_report.csv` | Road and inlet counts, configured minimum spacing and elapsed time. |

### CityCAT file format

```text
0.30 0.30 0.50
0 425194.7739 564040.3189 1
1 425498.7496 564098.4754 1
```

The first line contains the three values from `CITYCAT_HEADER`, in the order `size_x size_y friction`. Each subsequent line contains `id x y active_flag`. Coordinates use the projected CRS in metres. An `active_flag` of `1` marks an active inlet.

## Troubleshooting

- **Study area CRS is undefined.** Assign the study polygon its correct source CRS before running SSDGen.
- **No roads remain after filtering.** Check that the road layer overlaps the study area, inspect the `highway` values, and adjust `INCLUDE_HIGHWAYS` or `INCLUDE_LOCALS` as needed.
- **OSMnx is missing.** Install `osmnx` in the environment running SSDGen, or supply a local road line layer.
- **Too many or too few inlets.** Adjust `SPACING_MAP` to change candidate density along each kerb. Adjust `MIN_SSD_SPACING` to change the minimum separation across the full study area.

## Reproducibility

Keep `ssd_config.json` and `ssd_qc_report.csv` with the generated layers. Record the script version, input dataset versions and Python package versions. For OSM-based runs, retain the road cache with the results so that later changes to OpenStreetMap data do not affect a repeat run.

The final script corrects a Shapely 2 compatibility issue that could leave minimum spacing unenforced in earlier versions. Final inlet counts may therefore differ from earlier Shapely 2 runs with the same inputs.

## Associated paper

“Development of Synthetic Storm Drain Networks for Urban Flood Risk Modelling in Data-Limited Areas.” Accepted for publication in the *Journal of Flood Risk Management*, special issue “Physics-based Urban Flood Modelling”.

## Contributing

Issues and pull requests are welcome. For a bug report, include the Python and package versions, relevant configuration settings and, where possible, a small dataset that reproduces the issue. For changes to inlet placement or filtering, describe the expected effect on the output.

## License

SSDGen Pro is distributed under the [Apache License 2.0](LICENSE).
