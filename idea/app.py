import sys
import pandas as pd

# Get the file paths from the command line arguments
files = sys.argv[1:]

# Split the file paths into conditions and replicates
conditions = []
replicates = []
for file_path in files:
    condition, replicate = file_path.split("/")[-1].split("_")
    conditions.append(condition)
    replicates.append(replicate[:-13])  # Remove ".genes.results" from the replicate

# Create a DataFrame with the conditions and replicates
samples = pd.DataFrame({"run": replicates, "condition": conditions}, index=replicates)

# Initialize an empty DataFrame
df_final = pd.DataFrame()

# Loop over each file
for i, file in enumerate(files):
    # Read the file into a DataFrame
    df = pd.read_csv(file, sep="\t", usecols=[0, 4], names=["gene_id", samples.iloc[i, 0]], header=None)

    # If this is the first file, copy it to df_final
    if df_final.empty:
        df_final = df
    else:
        # Otherwise, merge the new DataFrame with df_final
        df_final = pd.merge(df_final, df, on="gene_id", how="outer")

# Filter rows based on the sum of values, excluding the 'gene_id' column
numeric_cols = df_final.columns[1:]  # All columns except 'gene_id'
df_final[numeric_cols] = df_final[numeric_cols].apply(pd.to_numeric, errors='coerce')
min_sum = 10  # Replace with your desired minimum sum
filtered_df = df_final.loc[df_final[numeric_cols].sum(axis=1) >= min_sum, :]

# Write the filtered DataFrame to a new file
filtered_df.to_csv("filtered_merged_file.csv", index=False)

# Final part is actually implementing deseq, need to get out "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj" from the dds
