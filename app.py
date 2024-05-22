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

# Write the final DataFrame to a CSV file
df_final.to_csv("merged_file.csv", index=False)

# I believe this has the same data as the dds from the lab, just as a csv, not sure what txi is necessary for, if someone could figure that part out.

# Have to include filtering out low counts

#Final part is actually implementing deseq, need to get out "gene_id""baseMean","log2FoldChange","lfcSE","stat","pvalue","padj" from the dds
