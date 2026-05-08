import argparse
import sys
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Build UKBB reference panel description table.")
    parser.add_argument(
        "--panel-file",
        required=True,
        help="TSV file with columns: Sample, Panel (collapsed ancestry)",
    )
    parser.add_argument(
        "--metadata-file",
        required=True,
        help="TSV with columns: Sample, Population code, Population, Superpopulation code, Superpopulation",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV/CSV/XLSX file path for the detailed reference panel table.",
    )

    args = parser.parse_args()

    # Load files (use pd.read_table as in daniel_code.py to guarantee correct parsing)
    required_meta_cols = {"Sample", "Population code", "Population", "Superpopulation code", "Superpopulation"}
    try:
        df_panel = pd.read_table(args.panel_file)
    except:
        print("ERROR: Cannot read panel file.")
        sys.exit(1)

    # rfmix_prelim_sample_map.tsv uses '#Sample'; align to 'Sample' expected by merge
    if "#Sample" in df_panel.columns and "Sample" not in df_panel.columns:
        df_panel = df_panel.rename(columns={"#Sample": "Sample"})

    required_panel_cols = {"Sample", "Panel"}
    missing_panel = required_panel_cols - set(df_panel.columns)
    if missing_panel:
        print(f"ERROR: Panel file missing columns: {', '.join(sorted(missing_panel))}")
        sys.exit(1)

    try:
        df_meta = pd.read_table(args.metadata_file)
        missing_meta = required_meta_cols - set(df_meta.columns)
        if missing_meta:
            print(f"ERROR: Metadata file missing columns: {', '.join(sorted(missing_meta))}")
            sys.exit(1)
    except:
        print("ERROR: Cannot read metadata file.")
        sys.exit(1)

    # rfmix_prelim_sample_map.tsv uses '#Sample'; align to 'Sample' expected by merge
    if "#Sample" in df_panel.columns and "Sample" not in df_panel.columns:
        df_panel = df_panel.rename(columns={"#Sample": "Sample"})

    # Merge on Sample ID
    merged = df_meta.merge(df_panel, on="Sample", how="left")

    # Check for missing ancestry labels
    missing = merged[merged["Panel"].isna()]
    if not missing.empty:
        print("\nWARNING: Some samples do not have a collapsed ancestry assignment:\n")
        print(missing[["Sample", "Population", "Superpopulation"]].head())
        print("\n → Fill missing ancestry in your input file.\n")

    # Count N individuals per population
    pop_counts = merged.groupby("Population")["Sample"].count().reset_index()
    pop_counts.columns = ["Population", "N_individuals"]

    merged = merged.merge(pop_counts, on="Population", how="left")

    # Add empty column for source dataset (to fill manually)
    merged["Source_dataset"] = ""

    # Reorder columns
    final = merged[
        [
            "Sample",
            "Population code",
            "Population",
            "Superpopulation",
            "Panel",               # Collapsed ancestry
            "N_individuals",
            "Source_dataset"      # to fill manually (1000G / HGDP / SGDP / GenomeAsia / Ancient)
        ]
    ]

    # Sort by ancestry and population
    final = final.sort_values(by=["Panel", "Population"])

    # Output
    if args.output.endswith((".csv", ".tsv")):
        sep = "\t" if args.output.endswith(".tsv") else ","
        final.to_csv(args.output, sep=sep, index=False)
    else:
        final.to_excel(args.output, index=False)

    print(f"\n✓ Reference panel table written to: {args.output}\n")


if __name__ == "__main__":
    main()
