import pandas as pd
import numpy as np
import os
from pathlib import Path

# ------------------------
# Setup paths
# ------------------------
epi = "/home/aliya/Liver/1111/Gene_epigenomics_v0.csv"
pseudo_dir = "/home/tooba/Skin_Cell_Model/pseudobulk/"               
output_dir = "/home/tooba/Skin_Cell_Model/processed_IbG"         

os.makedirs(output_dir, exist_ok=True)

ibg_files = sorted([f for f in os.listdir(pseudo_dir) if f.endswith('.tsv')])
print(f"Found {len(ibg_files)} IbG files")

# Track statistics
successful_files = 0
skipped_existing = 0
skipped_errors = 0

# ------------------------
# Load & prepare epigenomics file
# ------------------------
pred_df = pd.read_csv(epi)
print("Original columns:", pred_df.columns.tolist()[:10])

# Rename numeric columns to feature_N format
rename_dict = {}
for col in pred_df.columns:
    if col not in ['gene_name', 'chromo']:
        try:
            num = int(col)
            rename_dict[col] = f'feature_{num}'
        except (ValueError, TypeError):
            pass

pred_df = pred_df.rename(columns=rename_dict)
print("Renamed columns:", pred_df.columns.tolist()[:10])

# ------------------------
# Process each IbG file
# ------------------------
for ibg_file in ibg_files:
    ibg_path = os.path.join(pseudo_dir, ibg_file)

    # Extract cell type from filename
    # Handles formats: IbG_{cell_type}.tsv OR processed_IbG_{cell_type}.tsv OR {cell_type}.tsv
    stem = Path(ibg_file).stem
    if stem.startswith("processed_IbG_"):
        cell_type = stem[len("processed_IbG_"):]
    elif stem.startswith("IbG_"):
        cell_type = stem[len("IbG_"):]
    else:
        cell_type = stem

    out_filename = f"processed_IbG_{cell_type}.csv"   # Changed
    out_path = os.path.join(output_dir, out_filename)

    if os.path.exists(out_path):
        print(f"Skipping (already exists): {out_filename}")
        skipped_existing += 1
        continue

    print(f"Processing: {ibg_file} → {out_filename}")

    try:
        # ------------------------
        # Load expression data
        # ------------------------
        train_expr = pd.read_csv(ibg_path, sep='\t')
        train_expr["gene_name"] = train_expr["gene_name"].astype(str).str.strip('"')

        # Average across all donor columns
        donor_cols = [c for c in train_expr.columns if c != "gene_name"]
        if len(donor_cols) == 0:
            print(f"  ⚠ WARNING: No donor columns found in {ibg_file}, skipping")
            skipped_errors += 1
            continue

        train_expr[donor_cols] = train_expr[donor_cols].apply(pd.to_numeric, errors="coerce")
        train_expr["mean_expression"] = train_expr[donor_cols].mean(axis=1)
        train_expr = train_expr[["gene_name", "mean_expression"]]

        # Sanity check
        nonzero = (train_expr["mean_expression"] > 0).sum()
        print(f"  Genes with non-zero ctrl mean expression: {nonzero} / {len(train_expr)}")

        # ------------------------
        # Merge with epigenomics
        # ------------------------
        merged_df = pred_df.merge(
            train_expr[['gene_name', 'mean_expression']],
            on='gene_name',
            how='left'
        )

        merged_df['percentile'] = merged_df['mean_expression'].rank(method='average', pct=True)
        merged_df['percentile'] = merged_df['percentile'].fillna(0)

        # Sanity check
        print(f"  Percentile stats: min={merged_df['percentile'].min():.3f}, "
              f"max={merged_df['percentile'].max():.3f}, "
              f"zeros={(merged_df['percentile'] == 0).sum()}")

        # ------------------------
        # Reorder columns
        # ------------------------
        feature_cols = sorted(
            [c for c in merged_df.columns if c.startswith("feature_")],
            key=lambda x: int(x.split("_")[1])
        )

        merged_df = merged_df[['gene_name', 'chromo'] + feature_cols + ['percentile']]

        # ------------------------
        # Save
        # ------------------------
        merged_df.to_csv(out_path, index=False)
        print(f"  ✓ Saved: {out_filename} (Shape: {merged_df.shape})")
        successful_files += 1

    except Exception as e:
        print(f"  ⚠ ERROR processing {ibg_file}: {str(e)}")
        skipped_errors += 1
        continue

print(f"\n{'='*60}")
print(f"PROCESSING COMPLETE")
print(f"{'='*60}")
print(f"Successfully generated: {successful_files} new CSV files")
print(f"Skipped (already exist): {skipped_existing} files")
print(f"Skipped (errors):        {skipped_errors} files")
print(f"\nOutput directory: {output_dir}")
print(f"{'='*60}")
