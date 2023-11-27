import pathlib
import sqlite3
from multiprocessing import Pool

import pandas as pd

paths = list(
    pathlib.Path("/cluster/work/projects/ec85/joint-wind/model-aggregated/").rglob(
        "results.db"
    )
)

paths2 = []
for path in paths:
    paths2.append(str(path))


def convertdbparq(paths):
    path2 = pathlib.Path(paths)
    p = path2.parent / "results"
    # print(p)
    p.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(path2)
    c = conn.cursor()

    for table in c.execute(
        "SELECT name FROM sqlite_master WHERE type='table';"
    ).fetchall():
        t = table[0]
        # print(t)
        pd.read_sql("SELECT * from " + t, conn).to_parquet(
            p / (t + ".parquet"), compression="zstd"
        )


processes = 10

with Pool(processes) as pool:
    processed = pool.map(convertdbparq, paths2)
