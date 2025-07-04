import os
import pandas as pd
import sys

# Define paths
result_folder = sys.argv[1]
counts_file = os.path.join(result_folder, 'ancestry_counts', 'final_counts_aggregated.csv')
output_folder = os.path.join(result_folder, 'keep_files')
os.makedirs(output_folder, exist_ok=True)

# Load the aggregated ancestry count file
print(f"🔍 Loading: {counts_file}")
df_counts = pd.read_csv(counts_file)

# Set '#IID' as index and clean data
df_counts = df_counts.set_index('#IID')
df_counts = df_counts.dropna()         # Remove individuals with any missing value
df_counts = df_counts.astype(int)      # Ensure counts are integers
df_counts = df_counts.reset_index()    # Move '#IID' back to regular column

# Generate one keep file for each ancestry
for ancestry in df_counts.columns[1:]:
    print(f"▶ Generating keep file for: {ancestry}")
    
    # Select individuals for which this ancestry has the maximum count
    max_ancestry_samples = df_counts.loc[
        df_counts[ancestry] == df_counts.iloc[:, 1:].max(axis=1), '#IID'
    ]

    # Build FID/IID format DataFrame
    keep_df = pd.DataFrame()
    keep_df['#IID'] = max_ancestry_samples

    # Write keep file
    keep_file = os.path.join(output_folder, f"{ancestry}_keep.txt")
    keep_df.to_csv(keep_file, index=False, header=True, sep='\t')
    print(f"✅ Saved: {keep_file} ({len(keep_df)} individuals)")

print("🎉 Done. All keep files written to:", output_folder)
