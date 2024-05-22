import pandas as pd


files = ["../public/lab4/Chow_Rep1.genes.results",
         "../public/lab4/Chow_Rep2.genes.results",
         "../public/lab4/Chow_Rep3.genes.results",
         "../public/lab4/HFD_Rep1.genes.results",
         "../public/lab4/HFD_Rep2.genes.results",
         "../public/lab4/HFD_Rep3.genes.results"]

# This might not be necessary, had it in my lab though...
conditions = ["Chow"] * 3 + ["HFD"] * 3

samples = pd.DataFrame({"run": [f.split("/")[-1] for f in files],
                        "condition": conditions},
                       index=[f.split("/")[-1] for f in files])


# Initialize an empty DataFrame
df_final = pd.DataFrame()

# Loop over each file
for i, file in enumerate(files):
    # Read the file into a DataFrame
    df = pd.read_csv(file, sep="\t", usecols=[0, 5], names=["gene_id", f"tpm_{i+1}"], header=None)
    
    # If this is the first file, copy it to df_final
    if df_final.empty:
        df_final = df
    else:
        # Otherwise, merge the new DataFrame with df_final
        df_final = pd.merge(df_final, df, on="gene_id")

# Write the final DataFrame to a CSV file
df_final.to_csv("merged_file.csv", index=False)

#Do I need to set lengths to 0, this was in read me:
# 2. Note, we also used txi$length[txi$length == 0] <- 1 


#This gets a file in the same format as the RSEM file needed to run Deseq
