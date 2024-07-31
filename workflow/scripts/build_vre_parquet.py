import pandas as pd

#pd.read_csv(snakemake.input[0], index_col=["time", "technology", "spatial"]).to_parquet(
#    snakemake.output[0], compression="zstd"
#)

df = pd.read_csv(
    snakemake.input[0], 
    index_col=["time", "technology", "spatial"],
    dtype={
    0: float
    },
    header=0,  # Ensure the first row is treated as a header
)

# Convert to Parquet with zstd compression
df.to_parquet(
    snakemake.output[0], 
    compression="zstd"
)
