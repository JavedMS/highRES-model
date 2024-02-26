import pathlib
import sqlite3

import pandas as pd

path2 = pathlib.Path(snakemake.input[0])
p = path2.parent / "results"
p.mkdir(parents=True, exist_ok=True)

conn = sqlite3.connect(path2)
c = conn.cursor()

for table in c.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall():
    t = table[0]
    # print(t)
    pd.read_sql("SELECT * from " + t, conn).to_parquet(
        p / (t + ".parquet"), compression="zstd"
    )
