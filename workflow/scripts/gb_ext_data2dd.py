import calendar
import datetime
import pathlib

import pandas as pd

from data2dd_funcs import euro_demand2dd, scen2dd, temporal2dd, trans_links

root = pathlib.Path(snakemake.output[0]).parent     #TO extract the root/parent directory of output
data_root = root
out = pathlib.Path(".")

pscens = ["BASE"]


for psys in pscens:
    psys_scen = psys
    esys_scen = "BASE"

    scen_db = snakemake.input[1]
    # pathlib.Path(pathlib.PureWindowsPath(
    # "C:\\science\\studies\\innopaths\\highres\\scenarios"
    # "\\gb_ext_scenarios.xls"))
    f_techno = snakemake.input[2]

    params_to_write = {}

    params_to_write["gen"] = {}
    params_to_write["gen"]["set"] = [
        "g",
        "non_vre",
        "vre",
        "hydro_res",
        "uc_int",
        "uc_lin",
    ]
    params_to_write["gen"]["parameter"] = {}
    params_to_write["gen"]["parameter"]["all"] = [
        "emis fac",
        "max ramp",
        "min down",
        "min up",
        "startup cost",
        "inertia",
        "unit size",
        "capex2050",
        "varom",
        "fom",
        "fuelcost2050",
        "cap2area",
        "af",
        "min gen",
        "peak af",
    ]

    params_to_write["store"] = {}
    params_to_write["store"]["set"] = ["s", "uc_lin"]
    params_to_write["store"]["parameter"] = {}
    params_to_write["store"]["parameter"]["all"] = [
        "max_freq",
        "max_res",
        "eff_in",
        "eff_out",
        "loss_per_hr",
        "p_to_e",
        "varom",
        "e capex2050",
        "p capex2050",
        "fom",
        "af",
        "min down",
        "min up",
        "startup cost",
        "inertia",
        "unit size",
        "min gen",
        "max ramp",
    ]

    zones = pd.read_csv(
        # data_root/"zonal_def"/"zones.csv"
        snakemake.input[0]
    ).loc[:, "zone"]  # .values

    scen2dd(
        snakemake.output[1],
        root,
        f_techno,
        params_to_write,
        scen_db,
        psys_scen,
        esys_scen,
        zones,
        out=out,
        esys_cap=False,
        exist_cap=True,
    )

    trans_links(root, f_techno, out=out)


years = [snakemake.wildcards.year]
date_range=snakemake.params.date_range

for yr in years:
    
    date_range=[yr+"-"+date for date in date_range]
    yr = int(yr)

    dstart = datetime.datetime.fromisoformat(date_range[0])
    dstop = (datetime.datetime.fromisoformat(date_range[1])
            +datetime.timedelta(hours=23))

    #dstart = datetime.datetime(yr, 1, 1, 0)
    #dstop = datetime.datetime(yr, 12, 31, 23)

    # 2012 is close to central over the 2010-2016 period

    if calendar.isleap(yr):
        rleap = False
    else:
        rleap = True

    print(rleap)
    temporal2dd(dstart, dstop, root / out, snakemake.output[0])

    euro_demand2dd(
        snakemake.input["europedemandcsvlocation"],
        snakemake.input["europecountriescsvlocation"],
        root,
        root / out,
        dstart,
        dstop,
        scen_db,
        esys_scen,
        yr,
        "norescale",
    )
