import pandas as pd

pd.read_parquet(snakemake.input[0]).to_csv(snakemake.output["csvgdx"])
