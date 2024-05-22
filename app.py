import pandas as pd


files = ["../public/lab4/Chow_Rep1.genes.results",
         "../public/lab4/Chow_Rep2.genes.results",
         "../public/lab4/Chow_Rep3.genes.results",
         "../public/lab4/HFD_Rep1.genes.results",
         "../public/lab4/HFD_Rep2.genes.results",
         "../public/lab4/HFD_Rep3.genes.results"]

conditions = ["Chow"] * 3 + ["HFD"] * 3

samples = pd.DataFrame({"run": [f.split("/")[-1] for f in files],
                        "condition": conditions},
                       index=[f.split("/")[-1] for f in files])


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

#Final part is actually implementing deseq, need to get out "gene_id""baseMean","log2FoldChange","lfcSE","stat","pvalue","padj" from the csv
