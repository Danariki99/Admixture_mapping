#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate LocusZoom plots for the exact lead SNPs already saved in probs_folder.

Workflow
- Read lead SNPs from: probs/{hit}/{ancestry}/*.tsv  (your pipeline's shortlist: 1 per ancestry)
- Build LocusZoom input (CHR BP SNP P) from the matching PLINK2 results
  (filtered to TEST == ADD by default)
- Call locuszoom-cli to produce one HTML per lead SNP

Requirements
- pandas, numpy
- locuszoom-cli in PATH (conda: conda install -c conda-forge locuszoom-cli)
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import shutil
import subprocess

# ---- Default paths: keep aligned with your existing project structure ----
HIT_FOLDER_NAME_DEFAULT = "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_ancestries_PCA_verbose"
PROBS_FOLDER_DEFAULT    = "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/probs"
PLOTS_OUT_DEFAULT       = "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/plots/locuszoom"

# Map NAT -> AMR (for LD panel)
LD_PANEL_BY_ANCESTRY = {
    "AFR": "1000G_AFR",
    "EAS": "1000G_EAS",
    "EUR": "1000G_EUR",
    "SAS": "1000G_SAS",
    "WAS": "1000G_SAS",  # closest available reference
    "NAT": "1000G_AMR"   # NAT -> AMR
}

def parse_args():
    ap = argparse.ArgumentParser(description="Generate LocusZoom plots for lead SNPs found in probs_folder.")
    ap.add_argument("--probs-folder", default=PROBS_FOLDER_DEFAULT,
                    help="Root folder with probs/{hit}/{ancestry}/*.tsv (lead SNPs).")
    ap.add_argument("--hit-folder",   default=HIT_FOLDER_NAME_DEFAULT,
                    help="Root folder with hit subfolders containing PLINK2 outputs.")
    ap.add_argument("--out",          default=PLOTS_OUT_DEFAULT,
                    help="Output root for LocusZoom TSV/HTML files.")
    ap.add_argument("--build",        default="GRCh38", choices=["GRCh37", "GRCh38"],
                    help="Genome build to use for LocusZoom LD server.")
    ap.add_argument("--test-filter",  default="ADD",
                    help="Which TEST row to keep from PLINK results (default: ADD).")
    ap.add_argument("--ld-panel",     default=None,
                    help="Override LD panel (e.g., 1000G_EUR). If not set, it is inferred from ancestry.")
    ap.add_argument("--dry-run",      action="store_true",
                    help="Do not call locuszoom-cli; only write TSV input files.")
    ap.add_argument("--prefer-rsid",  action="store_true",
                    help="Kept for compatibility; IDs are already normalized (rsID or CHR:POS).")
    ap.add_argument("--only-hit",     default=None,
                    help="Process a single hit (exact folder name under probs/).")
    return ap.parse_args()

def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def pick_ld_panel(ancestry: str, override=None) -> str:
    """Return LD panel name for LocusZoom --source, possibly overridden."""
    if override:
        return override
    return LD_PANEL_BY_ANCESTRY.get(ancestry, "1000G_EUR")

def read_plink_results(plink_path: str, test_filter: str = "ADD"):
    """
    Read a PLINK2 .glm.* file with known columns, filter TEST, and guarantee an ID.

    Expected columns (subset):
    - '#CHROM', 'POS', 'ID', 'TEST', 'P', ...
    We keep only TEST == test_filter (default: ADD), i.e., the additive genetic effect.

    Returns
    -------
    df : pd.DataFrame (filtered to TEST==test_filter; with valid 'ID')
    col_chr, col_pos, col_id, col_p : str
        Canonical column names for downstream ('#CHROM', 'POS', 'ID', 'P')
    """
    # Tab-separated PLINK output
    df = pd.read_csv(plink_path, sep="\t")

    # Keep only the genetic-effect rows
    if "TEST" in df.columns:
        df = df[df["TEST"].astype(str).str.upper() == str(test_filter).upper()].copy()
        if df.empty:
            raise ValueError(f"No rows with TEST={test_filter} in {plink_path}")
    # If TEST is missing, assume file already filtered (rare case)

    # Check required columns
    required = ["#CHROM", "POS", "P"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns {missing} in {plink_path}")

    # Guarantee an ID (prefer rsID if present and valid; fallback to CHR:POS)
    if "ID" not in df.columns:
        df["ID"] = df["#CHROM"].astype(str) + ":" + df["POS"].astype(str)
    else:
        invalid = df["ID"].isna() | (df["ID"].astype(str).isin([".", ""]))
        if invalid.any():
            df.loc[invalid, "ID"] = (
                df.loc[invalid, "#CHROM"].astype(str) + ":" +
                df.loc[invalid, "POS"].astype(str)
            )

    return df, "#CHROM", "POS", "ID", "P"

def build_locuszoom_input(df: pd.DataFrame, col_chr: str, col_pos: str, col_id: str, col_p: str) -> pd.DataFrame:
    """
    Create a LocusZoom 'metal' table with canonical columns: CHR BP SNP P.
    """
    out = pd.DataFrame({
        "CHR": df[col_chr].astype(str),
        "BP":  pd.to_numeric(df[col_pos], errors="coerce"),
        "SNP": df[col_id].astype(str),
        "P":   pd.to_numeric(df[col_p], errors="coerce")
    }).dropna(subset=["BP", "P"])

    # Remove non-positive or non-finite p-values
    out = out[(out["P"] > 0) & np.isfinite(out["P"])]
    return out[["CHR", "BP", "SNP", "P"]]

def call_locuszoom(lz_tsv: str, refsnp: str, build: str, ld_source: str, prefix: str) -> None:
    """
    Run locuszoom-cli for a single refsnp.
    """
    if shutil.which("locuszoom") is None:
        print("WARNING: locuszoom-cli not found in PATH. Skipping CLI call.\n"
              "Install with: conda install -c conda-forge locuszoom-cli")
        return

    cmd = [
        "locuszoom", "--metal", lz_tsv,
        "--markercol", "SNP", "--pvalcol", "P", "--chrcol", "CHR", "--poscol", "BP",
        "--refsnp", refsnp,
        "--build", build,
        "--source", ld_source,
        "--plotonly",
        "--prefix", prefix
    ]
    print("Running:", " ".join(cmd))
    subprocess.check_call(cmd)

def try_find_plink_file(hit_folder: str, hit: str, ancestry: str, pheno: str):
    """
    Try common PLINK output filenames for robustness.
    """
    candidates = [
        f"{hit}_output.{ancestry}.{pheno}.glm.logistic.hybrid",
        f"{hit}_output.{ancestry}.{pheno}.glm.logistic",
        f"{hit}_output.{ancestry}.{pheno}.glm.linear",
        f"{hit}_output.{ancestry}.{pheno}.glm"
    ]
    for name in candidates:
        path = os.path.join(hit_folder, hit, name)
        if os.path.exists(path):
            return path
    return None

def main():
    args = parse_args()
    ensure_dir(args.out)

    # Iterate hits listed in probs folder
    for hit in sorted(os.listdir(args.probs_folder)):
        if args.only_hit and hit != args.only_hit:
            continue

        hit_path = os.path.join(args.probs_folder, hit)
        if not os.path.isdir(hit_path):
            continue

        # Parse hit name: IMPANCESTRY_<pheno_with_underscores>_<chrom>
        hit_parts = hit.split("_")
        if len(hit_parts) < 3:
            print(f"WARNING: Unexpected hit name format: {hit}  (skipping)")
            continue

        imp_ancestry = hit_parts[0]
        chrom = hit_parts[-1]
        pheno = "_".join(hit_parts[1:-1])  # middle part(s) can contain underscores

        for ancestry in sorted(os.listdir(hit_path)):
            ancestry_dir = os.path.join(hit_path, ancestry)
            if not os.path.isdir(ancestry_dir):
                continue

            # Locate the matching PLINK2 results for this (hit, ancestry)
            plink_file = try_find_plink_file(args.hit_folder, hit, ancestry, pheno)
            if plink_file is None:
                print(f"WARNING: Missing PLINK file for {hit}/{ancestry} (pheno={pheno}) → skipping")
                continue

            # Read and prepare the LocusZoom input once per (hit, ancestry)
            try:
                df_plink, col_chr, col_pos, col_id, col_p = read_plink_results(
                    plink_file, test_filter=args.test_filter
                )
                lz_input_df = build_locuszoom_input(df_plink, col_chr, col_pos, col_id, col_p)
            except Exception as e:
                print(f"ERROR reading/preparing {plink_file}: {e}")
                continue

            # Write the metal TSV for this (hit, ancestry)
            out_dir = os.path.join(args.out, hit, ancestry)
            ensure_dir(out_dir)
            lz_tsv = os.path.join(out_dir, f"{hit}.{ancestry}.{pheno}.locuszoom_input.tsv")
            lz_input_df.to_csv(lz_tsv, sep="\t", index=False)
            print(f"Wrote LocusZoom input: {lz_tsv}  (rows={len(lz_input_df)})")

            # Decide LD source
            ld_source = pick_ld_panel(ancestry, args.ld_panel)

            # For each lead SNP found in probs/{hit}/{ancestry}/*.tsv → one plot
            for fn in sorted(os.listdir(ancestry_dir)):
                if not fn.endswith(".tsv"):
                    continue
                refsnp = os.path.splitext(fn)[0]  # filename without .tsv is the SNP key you saved
                prefix = os.path.join(out_dir, f"{refsnp}")

                print(f"[{hit}/{ancestry}] refsnp={refsnp}  →  {prefix}.html")
                if args.dry_run:
                    continue

                try:
                    call_locuszoom(
                        lz_tsv=lz_tsv,
                        refsnp=refsnp,
                        build=args.build,
                        ld_source=ld_source,
                        prefix=prefix
                    )
                except subprocess.CalledProcessError as e:
                    print(f"locuszoom-cli failed for {refsnp}: {e}")
                    # continue to next SNP

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted by user.")
        sys.exit(130)
