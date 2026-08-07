# ---
# title: "xenium subset tool"
# date: "22-July-2026"
# author: "Jingjing Zhao @10x Genomics"
# ---

import argparse
import os
import json
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import spatialdata
import spatialdata_io
import sopa.io.explorer
import zarr

# Set rectilinear chunking in zarr
zarr.config.set({"array.rectilinear_chunks": True})
warnings.filterwarnings("ignore")


def main():
    parser = argparse.ArgumentParser(description="Subset SpatialData using a polygon.")
    parser.add_argument("-i", "--input", required=True, help="Path to input XOA bundle")
    parser.add_argument("-p", "--polygon", required=True, help="Path to selection GEOJSON")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-z", "--zarr_out", help="Optional: Path to save subsetted .zarr")

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # 1. Load polygon geometry
    # -------------------------------------------------------------------------
    print(f"Extracting geometry from: {args.polygon}")
    polygon_geom = gpd.read_file(args.polygon)
    target_polygon = polygon_geom.geometry[0]
    
    # -------------------------------------------------------------------------
    # 2. Load SpatialData bundle
    # -------------------------------------------------------------------------
    print(f"Reading data: {args.input}")
    sdata = spatialdata_io.xenium(args.input, cells_as_circles=False)

    # -------------------------------------------------------------------------
    # 3. Two-Stage Spatial Query (Memory Optimized)
    # -------------------------------------------------------------------------
    minx, miny, maxx, maxy = target_polygon.bounds
    
    print("Stage 1/2: Pre-filtering via bounding box query...")
    bb_sdata = spatialdata.bounding_box_query(
        sdata,
        axes=("x", "y"),
        min_coordinate=[minx, miny],
        max_coordinate=[maxx, maxy],
        target_coordinate_system="global"
    )

    print("Stage 2/2: Fine-grained spatial polygon query...")
    query_sdata = spatialdata.polygon_query(
        bb_sdata,
        polygon=target_polygon,
        target_coordinate_system="global"
    )

    # -------------------------------------------------------------------------
    # 4. Rechunk Images to Prevent OOM during Export
    # -------------------------------------------------------------------------
    if query_sdata.images:
        print("Rechunking image layers to save memory...")
        for img_key, data_tree in query_sdata.images.items():
            for scale_key in data_tree.keys():
                da = data_tree[scale_key]["image"]
                # Chunk spatial dimensions y and x to 1024x1024
                data_tree[scale_key]["image"] = da.chunk({"c": -1, "y": 1024, "x": 1024})

    # -------------------------------------------------------------------------
    # 5. Synchronize AnnData Table & Cell Boundaries
    # -------------------------------------------------------------------------
    print("Synchronizing AnnData table with spatial cell boundaries...")
    query_sdata["table"].obs["region"] = "cell_boundaries"
    query_sdata.set_table_annotates_spatialelement("table", region="cell_boundaries", instance_key="cell_id")

    # Match cell IDs between boundary GeoDataFrame and AnnData observations
    boundary_ids = set(query_sdata["cell_boundaries"].index.astype(str))
    obs_cell_ids = query_sdata["table"].obs["cell_id"].astype(str)
    
    mask = obs_cell_ids.isin(boundary_ids)
    query_sdata["table"] = query_sdata["table"][mask].copy()

    common_ids = list(query_sdata["table"].obs["cell_id"].astype(str))
    query_sdata["cell_boundaries"] = query_sdata["cell_boundaries"].loc[
        query_sdata["cell_boundaries"].index.astype(str).isin(common_ids)
    ]

    # -------------------------------------------------------------------------
    # 6. Write to Xenium Explorer format
    # -------------------------------------------------------------------------
    print(f"Writing Explorer files to: {args.output}")
    if not os.path.exists(args.output):
        os.makedirs(args.output, exist_ok=True)

    sopa.io.explorer.write(
        args.output, 
        query_sdata, 
        shapes_key="cell_boundaries",
        lazy=True
    )

    # -------------------------------------------------------------------------
    # 7. Update experiment.xenium Metadata File
    # -------------------------------------------------------------------------
    print(f"Updating experiment.xenium file in {args.output}")
    source_xenium_path = os.path.join(args.input, "experiment.xenium")
    target_xenium_path = os.path.join(args.output, "experiment.xenium")
    
    keys_to_extract = [
        "major_version", "minor_version", "patch_version", "run_name",
        "run_start_time", "region_name", "preservation_method", "cassette_name",
        "slide_id", "panel_type", "panel_design_id", "panel_predesigned_id",
        "panel_designer_version", "panel_name", "panel_organism", "panel_tissue_type",
        "panel_num_targets_predesigned", "panel_num_targets_custom", "chemistry_version",
        "pixel_size", "instrument_sn", "instrument_sw_version", "z_step_size",
        "well_uuid", "calibration_uuid"
    ]
    
    source_data = {}
    if os.path.exists(source_xenium_path):
        try:
            with open(source_xenium_path, 'r') as f:
                source_data = json.load(f)
        except json.JSONDecodeError:
            print(f"Error: {source_xenium_path} is not valid JSON.")

    extracted_data = {k: source_data[k] for k in keys_to_extract if k in source_data}

    target_data = {}
    if os.path.exists(target_xenium_path):
        try:
            with open(target_xenium_path, 'r') as f:
                target_data = json.load(f)
        except json.JSONDecodeError:
            target_data = {}

    if 'analysis_sw_version' in target_data and isinstance(target_data['analysis_sw_version'], str):
        target_data['analysis_sw_version'] += "-subset2zarr"
        
    num_cells = query_sdata["table"].shape[0]
    transcripts_per_cell = float(np.mean(query_sdata["table"].X.sum(axis=1)))
    num_transcripts = float(query_sdata["table"].X.sum())
    
    transcripts_df = query_sdata.points["transcripts"].compute()
    num_transcripts_high_qv = int((transcripts_df["qv"] >= 20).sum())
    
    target_data.update({
        'num_cells': num_cells,
        'transcripts_per_cell': transcripts_per_cell,
        'num_transcripts': num_transcripts,
        'num_transcripts_high_qv': num_transcripts_high_qv
    })
    target_data.update(extracted_data)
    
    with open(target_xenium_path, 'w') as f:
        json.dump(target_data, f, indent=4)
    
    # -------------------------------------------------------------------------
    # 8. Save Subsetted SpatialData Zarr Storage (Optional)
    # -------------------------------------------------------------------------
    if args.zarr_out:
        print(f"Saving subsetted Zarr to: {args.zarr_out}")
        query_sdata.write(args.zarr_out, overwrite=True)

    print("Done!")

if __name__ == "__main__":
    main()
