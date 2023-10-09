"""
Extract cap2area parameter from ods file and only update output file if input differs.
Intermediate cheap step so snakemake does not run resource intensive rule
build_weather again if not necessary.
"""
import pathlib

import pandas as pd

gen = (
    pd.read_excel(snakemake.input["odsdatabase"], sheet_name="gen", header=1)
    .set_index("Technology Name (highRES)")
    .loc["Windonshore", "cap2area"]
)

genfilepath = pathlib.Path(snakemake.output["windcap2area"])

GEN = str(gen)

if genfilepath.is_file():
    if genfilepath.read_text(encoding="utf8") == GEN:
        pass
    else:
        genfilepath.write_text(GEN, encoding="utf8")
else:
    genfilepath.write_text(GEN, encoding="utf8")
