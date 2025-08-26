# SSD-generation
1) Study area selection → STUDY_AREA_PATH (or polygon layer in QGIS).

2) Buffers by road class → buffer_by_class with {primary:20, secondary:10}.

3) Dissolve → done inside buffer_by_class (union + dissolve).

4) Polygon → line → polygon_to_line().

5) SSD placement every 50 m → points_along_lines() with SPACING_M=50.

6) Remove SSDs inside buildings → spatial join within and filter.
