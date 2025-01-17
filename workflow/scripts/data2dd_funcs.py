import calendar
import datetime
import itertools
import os

# import salem
import pathlib
import random
import re
import warnings

import numpy as np
import pandas as pd

# import matplotlib.pyplot as plt
# import geopandas as gpd
import xarray as xr


def wrapdd(d, out_par, otype, outfile=""):
    d = np.atleast_2d(d)

    if otype == "set":
        top = np.array([["set"], [out_par + " /"]])
        bottom = np.array([["/"], [""]])

    if otype == "scalar":
        top = np.array([["scalar"], [out_par + " /"]])
        bottom = np.array([["/"], [""]])

    if otype == "parameter":
        top = np.array([["parameter", ""], [out_par + " /", ""]])
        bottom = np.array([["/", ""], ["", ""]])

    #    if d.shape[1]<2:
    #
    #        top=np.array([["scalar"],[out_par+" /"]])
    #        bottom=np.array([["/"],[""]])
    #    else:
    #        top=np.array([["parameter",""],[out_par+" /",""]])
    #        bottom=np.array([["/",""],["",""]])

    d = np.concatenate((top, d, bottom), axis=0)

    if outfile != "":
        np.savetxt(outfile, d, delimiter=" ", fmt="%s")
    else:
        return d


def data2dd(data, sets, outfile="", all_combin=False, rounddp=8):
    # rounddp=5

    # each set in sets needs to be 1D at the moment

    # sets and data must be in correct order -> last set must be column
    # headers, previous sets are rows

    if type(data) != "numpy.ndarray":
        data = np.array(data)
    sets = np.array(sets, dtype="object")

    if len(sets) > 0:
        data = np.round(data.astype(float), rounddp)
        if len(sets) == 1:
            sets_out = sets[0]
            sets_new = np.array(sets_out).astype(str)

            data_out = np.hstack((sets_new.reshape(-1, 1), data.reshape(-1, 1)))
        else:
            if all_combin:
                sets_out = [item for item in list(itertools.product(*sets))]

                sets_out = np.char.array(sets_out)

                sets_new = sets_out[:, 0]

                for i in range(sets_out.shape[1] - 1):
                    sets_new = sets_new + "." + sets_out[:, i + 1]

                data_out = np.hstack((sets_new.reshape(-1, 1), data.reshape(-1, 1)))

            else:
                lens = np.array([item.shape[0] for item in sets])

                sets_new = []
                for i, s in enumerate(sets):
                    if lens[i] != lens.max():
                        sets_new.append(np.repeat(s, lens.max()))
                    else:
                        sets_new.append(s)

                sets_new = np.array(
                    [".".join(item) for item in np.char.array(list(zip(*sets_new)))]
                )

                data_out = np.hstack((sets_new.reshape(-1, 1), data.reshape(-1, 1)))

    else:
        data_out = data.reshape(-1, 1)
    if outfile != "":
        np.savetxt(outfile, data_out, fmt="%s", delimiter=" ")
    else:
        return data_out


def temporal2dd(dstart, dend, opath, temporaloutputpath):
    # nyears=(dend.year-dstart.year)+1
    years = np.arange(dstart.year, dend.year + 1)
    ntime = ((dend - dstart).total_seconds() / 3600) + 1

    hr2yr = []
    for nyr, yr in enumerate(years):
        shour = ((datetime.datetime(yr, 1, 1, 0) - dstart).total_seconds() / 3600) + 1
        if dend.year == yr:
            ehour = ((dend - dstart).total_seconds() / 3600) + 1
        else:
            ehour = (
                (datetime.datetime(yr, 12, 31, 23) - dstart).total_seconds() / 3600
            ) + 1

        hrs = np.arange(shour - 1, ehour).astype(int)
        hr2yr.append(list(zip(np.repeat(nyr, hrs.shape[0]).astype(int), hrs)))

    hr2yr = np.char.array(np.vstack(hr2yr).astype(str))

    hr2yr = (hr2yr[:, 0] + "." + hr2yr[:, 1]).reshape(-1, 1)

    if dstart.year != dend.year:
        yr = str(dstart.year) + "-" + str(dend.year)
    else:
        yr = str(dstart.year)

    out = []

    out.append(wrapdd(np.arange(ntime).reshape(-1, 1).astype(int), "h", "set"))
    out.append(
        wrapdd(np.arange(years.shape[0]).reshape(-1, 1).astype(int), "yr", "set")
    )
    out.append(wrapdd(hr2yr, "hr2yr_map", "set"))

    np.savetxt(temporaloutputpath, np.concatenate(out, axis=0), fmt="%s")


def year2dem(demyrs, yrs, fout):
    leaps = []
    rest = []
    for yr in np.arange(demyrs[0], demyrs[1] + 1):
        dstart = datetime.datetime(yr, 1, 1, 0)
        dend = datetime.datetime(yr, 12, 31, 23)
        if ((dend - dstart).total_seconds() / 3600.0) + 1 > 8760.0:
            leaps.append(yr)
        else:
            rest.append(yr)

    out = []
    for yr in np.arange(yrs[0], yrs[1] + 1):
        if yr >= 2005:
            demyr = yr
        else:
            dstart = datetime.datetime(yr, 1, 1, 0)
            dend = datetime.datetime(yr, 12, 31, 23)
            if ((dend - dstart).total_seconds() / 3600.0) + 1 > 8760.0:
                demyr = random.choice(leaps)
            else:
                demyr = random.choice(rest)

        out.append([yr, demyr])

    out = pd.DataFrame(out, columns=["year", "demand_year"])

    out.to_csv(fout, index=False)


def sensitivity(f, out):
    sense_runs = pd.read_excel(f, sheet_name="sensitivity_scenario_names", skiprows=0)

    if np.isnan(sense_runs.loc[0, "All"]):
        s_all = [""]
    else:
        s_all = list(map(str.strip, sense_runs.loc[0, "All"].split(",")))

    for j, r in sense_runs.iterrows():
        name = r["Name"]
        s = list(map(str.strip, r["Sensitivities"].split(","))) + s_all

        outdd = []

        sense1d = pd.read_excel(f, sheet_name="sensitivity_param_1d", skiprows=0)

        sense1d = sense1d.loc[sense1d["Sensitivity"].isin(s), :]

        pars = np.unique(sense1d["Parameter"].values)

        for p in pars:
            pdd = []
            sense = sense1d.loc[sense1d["Parameter"] == p, :]
            for i, row in sense.iterrows():
                tech = np.atleast_1d(row["Technology"])
                lim = np.atleast_1d(row["Value"])

                pdd.append(data2dd(lim.T, [tech]))
            outdd.append(wrapdd(np.concatenate(pdd, axis=0), p, "parameter"))

        sense2d = pd.read_excel(f, sheet_name="sensitivity_constraint_2d", skiprows=0)

        sense2d = sense2d.loc[sense2d["Sensitivity"].isin(s), :]

        pars = np.unique(sense2d["Parameter"].values)

        for p in pars:
            pdd = []
            sense = sense2d.loc[sense2d["Parameter"] == p, :]
            for i, row in sense.iterrows():
                zones = row["Z1_1":"Z17"][~row.loc["Z1_1":"Z17"].isnull()].index.values

                # divide by 1000 if in GW

                lim = row["Z1_1":"Z17"][~row.loc["Z1_1":"Z17"].isnull()].values / 1e3
                tech = np.atleast_1d(row["Technology"])
                lt = np.atleast_1d(row["LimType"])

                pdd.append(data2dd(lim.T, [zones, tech, lt]))

            outdd.append(wrapdd(np.concatenate(pdd, axis=0), p, "parameter"))

        np.savetxt(
            out / ("sensitivity_" + name + ".dd"), np.concatenate(outdd), fmt="%s"
        )


def trans_links(root, f, out="work"):
    tech_type = "trans"

    links_allowed = pd.read_excel(f, sheet_name="transmission_allowed", skiprows=1)

    links_out = np.array(
        links_allowed["Zone1"]
        + "."
        + links_allowed["Zone2"]
        + "."
        + links_allowed["Tech"]
    )

    links_tech = pd.unique(links_allowed["Tech"])

    set_outdd = []
    param_outdd = []

    set_outdd.append(wrapdd(data2dd(links_tech, []), "trans", "set"))
    set_outdd.append(wrapdd(data2dd(links_out, []), "trans_links", "set"))

    # data2dd(links_out,[],outfile=root/out/(tech_type+"_links.dd"))

    out_data = ["links_dist", "links_cap"]

    for od in out_data:
        out_par = tech_type + "_" + od

        param_outdd.append(
            wrapdd(data2dd(links_allowed[od].values, [links_out]), out_par, "parameter")
        )

    params = pd.read_excel(f, sheet_name="transmission", skiprows=1)

    nans = params.isnull()
    params = params.where(~nans, other=0.0)

    out_data = ["loss", "line_capex", "sub_capex", "varom"]

    #    trans_set=params["Technology Name (highRES)"]

    params = params[params["Technology Name (highRES)"].isin(links_tech)]

    for od in out_data:
        out_par = tech_type + "_" + od

        param_outdd.append(
            wrapdd(
                data2dd(params[od], [params["Technology Name (highRES)"]]),
                out_par,
                "parameter",
            )
        )

    param_outdd = np.concatenate(param_outdd, axis=0)
    set_outdd = np.concatenate(set_outdd, axis=0)

    pad = np.repeat(np.array(""), set_outdd.shape[0]).reshape(set_outdd.shape[0], 1)

    outdd = np.concatenate((np.hstack((set_outdd, pad)), param_outdd), axis=0)

    #    outdd=np.concatenate(outdd,axis=0)

    np.savetxt(root / out / (tech_type + ".dd"), outdd, delimiter=" ", fmt="%s")


def water(root, d1, d2, q, out="work"):
    fw_abs = pd.read_csv(
        pathlib.Path(
            pathlib.PureWindowsPath(
                "C:\\science\\highRES\\data\\nexus\\water\\fw_scenarios_" + q + ".csv"
            )
        )
    )
    sw_abs = pd.read_csv(
        pathlib.Path(
            pathlib.PureWindowsPath(
                "C:\\science\\highRES\\data\\nexus\\water\\sw_scenarios.csv"
            )
        )
    )

    _, delta = datelist(d1, d2, "hours")

    # fw_abs.loc[:,"Water_abs (ML/year)"]=fw_abs.loc[
    #     :,"Water_abs (ML/year)"]*(delta/8760.)
    # sw_abs.loc[:,"Water_abs (ML/year)"]=sw_abs.loc[
    #     :,"Water_abs (ML/year)"]*(delta/8760.)

    # abs units changed from Ml to Kl on output

    for lvl in ["low", "med", "high"]:
        data = np.tile(
            (
                fw_abs.loc[fw_abs["scenario"] == lvl, "Water_abs (ML/year)"].values
                * 1e3
                / float(delta)
            ).reshape(
                fw_abs.loc[fw_abs["scenario"] == lvl, "sys_zone"].values.shape[0], -1
            ),
            (1, delta),
        )
        t = np.arange(delta)
        data2dd(
            data,
            [fw_abs["sys_zone"].drop_duplicates().values, t],
            all_combin=True,
            outfile=root / out / ("fw_abstract_lim_" + lvl + "_" + q + ".dd"),
            rounddp=0,
        )

    for lvl in ["low", "med", "high"]:
        data = np.tile(
            (
                sw_abs.loc[sw_abs["scenario"] == lvl, "Water_abs (ML/year)"].values
                * 1e3
                / float(delta)
            ).reshape(
                sw_abs.loc[sw_abs["scenario"] == lvl, "zone"].values.shape[0], -1
            ),
            (1, delta),
        )
        t = np.arange(delta)
        data2dd(
            data,
            [sw_abs["zone"].drop_duplicates().values, t],
            all_combin=True,
            outfile=root / out / ("sw_abstract_lim_" + lvl + "_" + q + ".dd"),
            rounddp=0,
        )


#    abstract=pd.read_excel(fin,sheet_name="fw_abs_potentials",skiprows=2)
#    consum=pd.read_excel(fin,sheet_name="sw_abs_potentials",skiprows=2)
#
#    _,delta=datelist(d1,d2,"hours")
#
#
#
# for lvl in ["low", "med", "high"]:
#     data = np.tile((abstract[lvl].values/float(delta)
#                     ).reshape(
#                         abstract[lvl].values.shape[0], -1), (1, delta))
#     t = np.arange(delta)
#     data2dd(data, [abstract["zone"].values, t], all_combin=True,
#             outfile=root/out/("fw_abstract_lim_"+lvl+".dd"), rounddp=1)

# for lvl in ["low", "med", "high"]:
#     data = np.tile((consum[lvl].values/float(delta)
#                     ).reshape(
#                         consum[lvl].values.shape[0], -1), (1, delta))
#     t = np.arange(delta)
#     data2dd(data, [consum["zone"].values, t], all_combin=True,
#             outfile=root/out/("sw_abstract_lim_"+lvl+".dd"), rounddp=1)


def co2lim2dd(co2budgetddlocation, root, run, esys, scen_db, out=""):
    scen = pd.read_excel(scen_db, sheet_name="scenario_co2_cap", skiprows=1)

    dout = scen.loc[scen["Esys Scenario"] == esys, 2050].values

    wrapdd(dout, "co2_budget", "scalar", outfile=co2budgetddlocation)


def getzlims(lim, techs, zones):
    lim = lim.loc[(lim["Year"] == 2050) & (lim["Technology"].isin(techs)), :]

    if lim.empty:
        return np.array([])

    have_lim = ~lim.loc[:, "limtype"].isnull()

    if (~have_lim).any():
        warnings.warn(
            "Warning, missing zonal new capacity limits for:"
            + ", ".join(lim.loc[~have_lim, "Technology"])
        )

    lim = lim.loc[have_lim, :]
    para_lim = lim["parameter"].drop_duplicates()

    outlims = []

    for p_lim in para_lim:
        nl = []
        for _, row in lim.loc[lim["parameter"] == p_lim, :].iterrows():
            zones = row[zones].index.values
            limval = row[zones].values.astype(float) / 1e3
            limval[np.isnan(limval)] = 0.0

            tech = np.atleast_1d(row["Technology"])
            limtype = np.atleast_1d(row["limtype"])

            nl.append(data2dd(limval.T, [zones, tech, limtype]))

        outlims.append(wrapdd(np.concatenate(nl, axis=0), p_lim, "parameter"))

    return np.concatenate(outlims, axis=0)


def scen2dd(
    co2budgetddlocation,
    root,
    fin,
    params2write,
    scen_db,
    run,
    esys,
    zones,
    out="work",
    esys_cap=False,
    exist_cap=False,
):
    co2lim2dd(co2budgetddlocation, root, run, esys, scen_db, out=out)

    scen = pd.read_excel(scen_db, sheet_name="scenario_tech_definition", skiprows=0)

    scen = scen.loc[(scen["Psys Scenario"] == run), :]

    techs = scen["Technology Name (highRES)"]
    tech_class = ["gen", "store"]

    for tech_type in tech_class:
        set_outdd = []
        param_outdd = []

        params = pd.read_excel(fin, sheet_name=None, skiprows=1)[tech_type]

        params = params[params["Technology Name (highRES)"].isin(techs)]

        sets = params2write[tech_type]["set"]
        data_out = params2write[tech_type]["parameter"]

        params = params.merge(scen, on="Technology Name (highRES)")

        # new_lim=pd.read_excel(fin,sheet_name=tech_type+"_lim_z",skiprows=0)

        param_outdd.append(
            getzlims(
                pd.read_excel(fin, sheet_name=tech_type + "_lim_z", skiprows=0),
                techs,
                zones,
            )
        )

        if exist_cap:
            lims = getzlims(
                pd.read_excel(fin, sheet_name=tech_type + "_exist_z", skiprows=0),
                techs,
                zones,
            )

            if lims.size != 0:
                param_outdd.append(
                    getzlims(
                        pd.read_excel(
                            fin, sheet_name=tech_type + "_exist_z", skiprows=0
                        ),
                        techs,
                        zones,
                    )
                )

        if esys_cap:
            dout = params.loc[
                (~params["Esys_capacity~2050"].isnull()), "Esys_capacity~2050"
            ]

            if dout.empty:
                continue
            else:
                wrapdd(
                    data2dd(
                        dout.values,
                        [
                            params[
                                (params["Technology Name (highRES)"] != "pgen")
                                & (~params["Esys_capacity~2050"].isnull())
                            ]["Technology Name (highRES)"].values
                        ],
                        rounddp=2,
                    ),
                    tech_type + "_fx_natcap",
                    root / out / (esys + "_" + tech_type + "_fx_natcap.dd"),
                )

        for s in sets:
            if s == "g" or s == "s":
                vals = params["Technology Name (highRES)"]

                set_outdd.append(wrapdd(data2dd(vals, []), s, "set"))
                continue

                # outsets.append(wrapdd(data2dd(vals,[])))
                # set_outdd.append(wrapdd(data2dd(vals,[],s,"set"))
            elif "uc" in s:
                vals = params.loc[(~params[s].isnull()), "Technology Name (highRES)"]

                # if vals.empty:
                #    continue

                set_outdd.append(
                    wrapdd(
                        data2dd(vals, []),
                        tech_type + "_" + s + "(" + sets[0] + ")",
                        "set",
                    )
                )
                # else:
                #   set_outdd.append(wrapdd(data2dd(
                #       vals,
                #       [],
                #       root/out/(run+"_"+tech_type+"_set_"+s+".dd"))))

            else:
                vals = params[params["set"] == s]["Technology Name (highRES)"]

                set_outdd.append(
                    wrapdd(data2dd(vals, []), s + "(" + sets[0] + ")", "set")
                )

                # data2dd(vals,[],root/out/(run+"_"+tech_type+"_set_"+s+".dd"))

        for s in data_out:
            for p in data_out[s]:
                if "varom" in p or "fuel" in p:
                    # rounddp=5
                    rounddp = 8
                elif "emis" in p:
                    # rounddp=3
                    rounddp = 8
                else:
                    # rounddp=2
                    rounddp = 8

                out_par = tech_type + "_" + str.replace(p, " ", "")

                if "capex" in p or "fuelcost" in p:
                    out_par = re.sub(r"\d+", "", out_par)

                vals = params.loc[~(params[p].isnull()), p].values
                t_out = params.loc[
                    ~(params[p].isnull()), "Technology Name (highRES)"
                ].values

                param_outdd.append(
                    wrapdd(
                        data2dd(vals, [t_out], rounddp=rounddp), out_par, "parameter"
                    )
                )

        param_outdd = np.concatenate(param_outdd, axis=0)
        set_outdd = np.concatenate(set_outdd, axis=0)

        pad = np.repeat(np.array(""), set_outdd.shape[0]).reshape(set_outdd.shape[0], 1)

        outdd = np.concatenate((np.hstack((set_outdd, pad)), param_outdd), axis=0)

        np.savetxt(
            root / out / (run + "_" + tech_type + ".dd"), outdd, delimiter=" ", fmt="%s"
        )

        # param_outdd=np.concatenate(param_outdd,axis=0)
        # set_outdd=np.concatenate(set_outdd,axis=0)

        # np.savetxt(
        #     root/out/(run+"_"+tech_type+"_parameters.dd"), param_outdd,
        #     delimiter=" ",
        #     fmt="%s")
        # np.savetxt(
        #     root/out/(run+"_"+tech_type+"_sets.dd"),
        #     set_outdd, delimiter=" ",
        #     fmt="%s")


def euro_demand2dd(
    europedemandcsvlocation,
    europecountriescsvlocation,
    dpath,
    opath,
    dstart,
    dstop,
    scen_db,
    esys_scen,
    yr,
    rescale="annual",
    zones="yes",
):
    # 2010 demand is ~ 3216 TWh over the 31 countries in ETM
    #

    countries = pd.read_csv(europecountriescsvlocation)
    # pd.read_csv(dpath/"zonal_def"/"europe_countries.csv")

    c_sel = countries.loc[countries["ETM"] == 1, "ISO2"]

    d = pd.read_csv(europedemandcsvlocation)

    d["datetime"] = pd.to_datetime(d["datetime"])
    d = d.set_index("datetime")

    d = d.loc[:, d.columns.isin(c_sel)]

    d = d[dstart:dstop]

    if d.shape[1] != c_sel.shape[0]:
        print("Countries missing...")

    d[(d == 0)] = np.nan

    d.interpolate(limit=2, inplace=True)

    if pd.isnull(d).any().any():
        for c in d.columns[pd.isnull(d).any(axis=0)]:
            print(d.loc[pd.isnull(d[c]), c])

        print("Zero demands for: ", d.columns[pd.isnull(d).any(axis=0)])

    if calendar.isleap(yr):
        d = pd.concat((d, d.iloc[0:24, :]))

    if rescale == "annual":
        out_flg = "annual"
        if d.shape[0] >= 8760.0:
            scen = pd.read_excel(scen_db, sheet_name="scenario_annual_dem", skiprows=0)

            euro31_dem = scen[scen["Esys Scenario"] == esys_scen][
                "Annual demand (2050)"
            ].iloc[0]

            d = d * (euro31_dem * 1e6 / d.sum().sum())

    #        else:
    #            dem_scaling=1.511
    #            demand=demand*dem_scaling

    else:
        out_flg = "norescale"

    t = np.arange(d.shape[0])
    z = d.columns.values

    wrapdd(
        data2dd(d.values.T, [z, t], all_combin=True),
        "demand",
        "parameter",
        outfile=opath / (esys_scen + "_" + out_flg + "_demand_" + str(yr) + ".dd"),
    )


def demand2dd(
    root,
    d1,
    d2,
    scen_db,
    uktm_scen,
    out="work",
    rescale="annual",
    zones="yes",
    yr_flg="",
    add_heat=False,
    add_ev=False,
):
    ndays = ((d2 - d1).days + 1) * 24
    print(d1)
    print(ndays)

    yr = d1.year

    demand_file = root / "demand" / "GB_demand_2002-2017.csv"

    dem = pd.read_csv(demand_file)

    dstart = datetime.datetime.strptime(
        dem["SETTLEMENT_DATE"].iloc[0], "%d/%m/%Y 00:00"
    )
    dend = datetime.datetime.strptime(
        dem["SETTLEMENT_DATE"].iloc[-1], "%d/%m/%Y 00:00"
    ) + datetime.timedelta(minutes=30.0 * (dem["SETTLEMENT_PERIOD"].iloc[-1] - 1))

    dem.index = pd.date_range(start=dstart, end=dend, freq="0.5H")

    delta = ((dend - dstart).days + 1) * 24

    # dates_hourly = [
    #     dstart+datetime.timedelta(minutes=60.*x) for x in range(delta)]

    if dem.shape[0] != delta * 2:
        print(dem.shape[0])
        print(delta * 2)
        print("Actual and expected demand time steps don't match, quitting")
        stop

    # indo + estimate of embedded wind/solar in MW

    dem["demand"] = (
        dem["ND"] + dem["EMBEDDED_WIND_GENERATION"] + dem["EMBEDDED_SOLAR_GENERATION"]
    )

    dem = dem["demand"]

    # Demand is average of settlement period from, e.g., 00:00:00 to 00:30:00.
    # Setup so resampled data is mean(23:30:00,00:00:00) => 00:00:00 demand
    # first settlement period is just duplicated
    #
    # This should ensure that the instantenous weather data aligns fairly
    # sensibly with the demand data.

    dem = pd.concat((dem.tshift(periods=-1).iloc[:1], dem.iloc[:-1]), axis=0)

    # half hours to hours

    dem = dem.resample("1H", closed="right", label="right").mean()

    # select hours from requested year

    dem = dem.loc[d1:d2]

    if int(dem.shape[0]) != ((d2 - d1).days + 1) * 24:
        print(dem.shape[0])
        print(((d2 - d1).days + 1) * 24)
        print(yr)
        print("Requested demand time series length doesn't match returned")
        return

    if rescale == "annual":
        # used to remove the current share of electricity demand from Northern
        # Ireland

        ireland = 0.974

        out_flg = "annual"

        # if dem.shape[0]>=8760.:

        scen = pd.read_excel(scen_db, sheet_name="scenario_annual_dem", skiprows=0)

        scen_dem = (
            scen[scen["Esys Scenario"] == uktm_scen]["Annual demand (2050)"].iloc[0]
            * ireland
        )

        # total demand over the selected time series -> could be longer than
        # 1 year

        dem_tot = np.sum(dem) / 1.0e6

        # scale up demand time series -> accounts for series longer than 1 year

        dem_scaling = (scen_dem * (dem.shape[0] / 8760)) / dem_tot
        dem = dem * dem_scaling

    #        else:
    #            dem_scaling=1.511
    #            #dem_scaling=(scen_dem*(dem.shape[0]/8760))/dem_tot
    #            dem=dem*dem_scaling

    elif rescale == "timesliced" and dem.shape[0] >= 8760:
        out_flg = "tsliced"

        uktm_dem = pd.read_excel(
            scen_db, sheet_name="uktm_elec_demand_timesliced", skiprows=0
        )
        uktm_dem = uktm_dem[uktm_dem["uktm_scenario"] == uktm_scen]
        # uktm_dem = pd.read_csv(
        #     root/"data"/"demand"/"lowghg_elc_demand_timesliced.csv")
        uktm_dem["demand_PJ"] = uktm_dem["demand_PJ"] * 277778.0
        uktm_dem.rename(columns={"demand_PJ": "demand_MWh"}, inplace=True)

        hour = demand.index.hour

        day = (hour > 7) & (hour <= 17)
        pek = (hour > 17) & (hour <= 20)
        eve = (hour > 20) | (hour <= 0)
        ngt = (hour > 0) & (hour <= 7)

        times = {}
        times["D"] = day
        times["P"] = pek
        times["E"] = eve
        times["N"] = ngt

        spr = (demand.index > datetime.datetime(yr, 3, 1, 0)) & (
            demand.index < datetime.datetime(yr, 6, 1, 1)
        )
        smr = (demand.index > datetime.datetime(yr, 6, 1, 0)) & (
            demand.index < datetime.datetime(yr, 9, 1, 1)
        )
        aut = (demand.index > datetime.datetime(yr, 9, 1, 0)) & (
            demand.index < datetime.datetime(yr, 12, 1, 1)
        )
        win = (demand.index > datetime.datetime(yr, 12, 1, 0)) | (
            demand.index < datetime.datetime(yr, 3, 1, 1)
        )

        seasons = {}
        seasons["P"] = spr
        seasons["S"] = smr
        seasons["A"] = aut
        seasons["W"] = win

        # plt.scatter(demand.index,demand,color="r")

        ireland = 0.974

        for i, ts in enumerate(uktm_dem["timeslice"]):
            select = seasons[ts[0]] & times[ts[1]]
            hist = dem[select].sum()
            uktm = uktm_dem["demand_MWh"].iloc[i]
            demand[select] = demand[select] * (uktm * ireland / hist)

        # print(demand.sum())
        # plt.scatter(demand.index,demand,color="g")
        # plt.show()

    elif rescale == "no":
        out_flg = "norescale"
        print("No demand rescaling done")

    if zones == "yes":
        # grid=root/"data"/"shapefiles"/"zones_wind_masked_intersect.shp"

        # shp=gpd.read_file(grid)

        # zones=np.unique(shp["sys_zones"])

        # zones = np.unique(arcpy.da.TableToNumPyArray(
        #     grid, "sys_zones").astype(str))

        dem_share = pd.read_csv(
            root / "demand" / "demand_share_sys_zones.csv", index_col=0
        )

        zones = dem_share.index

        # ind=[np.where(dem_share["zone"] == z)[0][0] for z in zones]

        # dem_share=dem_share.iloc[ind,1].values/100.

        demand_out = pd.DataFrame(
            dem.values.reshape(-1, 1) * dem_share.values.T / 100.0,
            index=dem.index,
            columns=zones,
        )

        # demand_out=np.zeros((zones.shape[0],dem.shape[0]))

        # demand_out[:,:]=dem.values.reshape(1,-1)*dem_share.reshape(-1,1)

        if add_heat:
            if d1.year != d2.year:
                heat = []
                for yy in np.arange(d1.year, d2.year + 1):
                    heat.append(
                        pd.read_csv(
                            root
                            / "demand"
                            / "heat"
                            / (uktm_scen + "_heat_demand_" + str(yy) + ".csv"),
                            index_col=0,
                            parse_dates=True,
                        )
                    )
                heat = pd.concat(heat)

            else:
                heat = pd.read_csv(
                    root
                    / "demand"
                    / "heat"
                    / (uktm_scen + "_heat_demand_" + str(d1.year) + ".csv"),
                    index_col=0,
                    parse_dates=True,
                )

            demand_out = demand_out + heat[d1:d2].loc[:, zones]

        if add_ev:
            if d1.year != d2.year:
                ev = []
                for yy in np.arange(d1.year, d2.year + 1):
                    ev.append(
                        pd.read_csv(
                            root
                            / "demand"
                            / "ev"
                            / (uktm_scen + "_ev_demand_" + str(yy) + ".csv"),
                            index_col=0,
                            parse_dates=True,
                        )
                    )
                ev = pd.concat(ev)
            else:
                ev = pd.read_csv(
                    root
                    / "demand"
                    / "ev"
                    / (uktm_scen + "_ev_demand_" + str(d1.year) + ".csv"),
                    index_col=0,
                    parse_dates=True,
                )

            demand_out = demand_out + ev[d1:d2].loc[:, zones]

        print(demand_out.sum().sum() / 1e6)
        print(demand_out.sum(axis=1).max())

        #        demand_out.loc[:,"Z14"].plot()
        #        plt.show()

        dem_csv = pd.DataFrame(demand_out, columns=zones)

        dem_csv.to_csv(
            pathlib.Path(
                pathlib.PureWindowsPath(
                    "C:\\science\\highRES\\gb\\data\\demand\\out\\"
                    + (uktm_scen + "_" + str(yr) + ".csv")
                )
            )
        )

        demand_out = data2dd(
            demand_out.T, [zones, np.arange(dem.shape[0])], all_combin=True, rounddp=2
        )

        if d1.year != d2.year:
            yr = str(d1.year) + "-" + str(d2.year)
        elif yr_flg != "":
            yr = yr_flg
        else:
            yr = str(d1.year)

        wrapdd(
            demand_out,
            "demand",
            "parameter",
            outfile=out / (uktm_scen + "_" + out_flg + "_demand_" + yr + ".dd"),
        )

        print(zones)

        zones.to_frame().to_csv(out / "zones.dd", header=False, index=False)

    elif zones == "no":
        data2dd(
            dem.values, [np.arange(dem.shape[0])], out / ("demand_" + str(yr) + ".dd")
        )


def gb_ext_demand2dd(f, agg, uk_shares, opath, dstart, dstop, scen_db, esys_scen, yr):
    d = pd.read_csv(f)

    agg = pd.read_csv(agg)

    c_sel = agg.loc[agg["gb_ext"] == 1, ["ISO2", "EU_simpl"]]

    d.index = pd.to_datetime(d["datetime"])
    d.drop("datetime", axis=1, inplace=True)

    d = d[dstart:dstop]

    d = (
        d.T.reset_index()
        .merge(c_sel, left_on="index", right_on="ISO2")
        .groupby("EU_simpl")
        .sum()
    )

    uk_shares = (
        pd.read_excel(uk_shares, sheet_name="uk_zonal_shares", index_col=0) / 100
    ).dot(d.loc[d.index == "UK", :].values)
    uk_shares.columns = d.columns

    d = d.append(uk_shares)
    d = d.loc[d.index != "UK", :]

    # temp = (uk_shares.loc[:, "share"].values.reshape(
    #     uk_shares.shape[0], 1)/100)*d.loc[d.index == "UK", :].values

    # d = d.loc[d.index.str.contains("UK") | d.index.str.contains("BNL"), :]

    t = np.arange(d.shape[1])
    z = d.index

    wrapdd(
        data2dd(d.values, [z, t], all_combin=True),
        "demand",
        "parameter",
        outfile=opath / (esys_scen + "_demand_" + str(yr) + ".dd"),
    )


def datelist(d1, d2, delta):
    if delta == "days":
        delta = int((d2 - d1).days + 1)
        return [d1 + datetime.timedelta(days=x) for x in range(delta)], delta
    if delta == "hours":
        delta = int(((d2 - d1).total_seconds() / 3600) + 1)
        return [d1 + datetime.timedelta(hours=x) for x in range(delta)], delta


def norestrict():
    vre_params = {}
    vre_params["Solar"] = [
        pathlib.Path(
            pathlib.PureWindowsPath(
                "C:\\science\\highRES\\data\\studies\\nexus\\areas\\base\\"
                "onshore_base_areas_zone_wgrid.csv"
            )
        ),
        "data/weather/solar/cmsaf/csv",
        "solarpv_cmsaf",
        0.09,
        40,
    ]
    vre_params["Windoffshore"] = [
        pathlib.Path(
            pathlib.PureWindowsPath(
                "C:\\science\\highRES\\data\\studies\\nexus\\areas\\offshore"
                "\\buildable\\offshore_tec_low_buildable_areas.csv"
            )
        ),
        "data/weather/wind/cfsr/csv",
        "offshore",
        0.2,
        5.0,
    ]
    vre_params["Windonshore"] = [
        pathlib.Path(
            pathlib.PureWindowsPath(
                "C:\\science\\highRES\\data\\studies\\nexus\\areas\\base"
                "\\onshore_base_areas_zone_wgrid.csv"
            )
        ),
        "data/weather/wind/cfsr/csv",
        "onshore",
        0.15,
        3.0,
    ]
    vre2dd_aggregated(
        root, d1, d2, vre_params, out=out, writexl=True, fmod="norestriction"
    )


def hydro_inflow(path, out, dstart, dend, sel_zones=None):
    tech = "HydroRes"

    ds_file = "hydro_res_inflow_" + str(dstart.year)

    ds = pd.read_csv(path / (ds_file + ".csv"))

    delta_h = str(int((dend - dstart).total_seconds() / 3600))

    ds = ds.loc[:, :delta_h]

    # cells=ds["zones"].values
    # zones=ds["zones"].values

    ds = ds.set_index("zones").stack().reset_index()
    ds = ds.rename(columns={"level_1": "h", "zones": "z", 0: "Value"})
    ds.insert(0, "Tech", np.repeat(tech, ds.shape[0]))
    ds = ds[["h", "z", "Tech", "Value"]]

    ds.to_csv(out / (ds_file + ".csv"), index=False)

    out_file = "hydro_inflow_" + str(dstart.year)

    os.system(
        "csv2gdx "
        + str(out / (ds_file + ".csv"))
        + " output="
        + str(out / (out_file + ".gdx"))
        + " ID=hydro_inflow Index=(1,2,3)"
        "Value=(4) UseHeader=True StoreZero=True"
    )
    # (out/(ds_file+'.csv')).unlink()


def hydro_res_inflow(path, out, dstart, dstop, agg="", rleap=False):
    tech = "HydroRes"

    # ds_file="hydro_res_inflow_"+str(dstart.year)
    ds_file = "hydro_res_inflow_2013"
    ds = pd.read_csv(path / (ds_file + ".csv"))

    print("Hydro inflow currently fixed to 2013!!")

    if rleap:
        ds = pd.pivot_table(ds, columns="zones")
        ds.index = ds.index.astype(int)

        ds = ds.sort_index()
        ds.index = pd.date_range(
            start=pd.Timestamp(2013, 1, 1, 0),
            end=pd.Timestamp(2013, 12, 31, 23),
            freq="H",
        )

        noleap = (ds.index >= pd.Timestamp(2012, 2, 29, 0)) & (
            ds.index <= pd.Timestamp(2012, 2, 29, 23)
        )

        ds = ds[~noleap]

        ds.index = np.arange(len(ds)).astype(str)

        ds = ds.T.reset_index()

    delta_h = str(int((dstop - dstart).total_seconds() / 3600))
    ds = ds.loc[:, :delta_h]

    if agg:
        agg = pd.read_csv(agg)
        c_sel = agg.loc[agg["gb_ext"] == 1, ["ISO2", "EU_simpl"]]

        ds = ds.merge(c_sel, left_on="zones", right_on="ISO2").groupby("EU_simpl").sum()

    ds = ds.stack().reset_index()
    ds = ds.rename(columns={"level_1": "h", "EU_simpl": "z", 0: "Value"})
    ds.insert(0, "Tech", np.repeat(tech, ds.shape[0]))
    ds = ds[["h", "z", "Tech", "Value"]]

    ds.to_csv(out / (ds_file + ".csv"), index=False)

    out_file = "hydro_res_inflow_" + str(dstart.year)

    os.system(
        "csv2gdx "
        + str(out / (ds_file + ".csv"))
        + " output="
        + str(out / (out_file + ".gdx"))
        + " ID=hydro_inflow Index=(1,2,3)"
        "Value=(4) UseHeader=True StoreZero=True"
    )
    (out / (ds_file + ".csv")).unlink()


def aggregate_zones(ds, agg, func):
    agg = pd.read_csv(agg)
    c_sel = agg.loc[agg["gb_ext"] == 1, ["ISO2", "EU_simpl"]]

    ds = ds.merge(c_sel, left_on="zones", right_on="ISO2").groupby("EU_simpl").agg(func)

    return ds


def vre2dd(
    out,
    dstart,
    dend,
    vre_params,
    sel_zones=None,
    zonal=True,
    fmod="",
    write_csv=True,
    csv_file="",
    agg="",
    rleap=False,
):
    min_area = 1

    cf_out = []
    area_out = []
    cells_out = []
    zones_out = []

    if write_csv:
        csv_out = []

    for tech in vre_params.keys():
        print(tech)

        cell_areas = vre_params[tech][0]
        data_path = vre_params[tech][1]
        cf_cutoff = vre_params[tech][2]

        years = "|".join(np.arange(dstart.year, dend.year + 1).astype(str))
        flist = pd.Series(os.listdir(data_path))
        flist = flist[flist.str.contains(years)]

        if len(flist) == 0 and tech != "HydroRoR":
            raise ValueError("No generation files found for " + tech)

        if tech == "HydroRoR":
            flist = pd.Series(os.listdir(data_path))
            flist = flist[flist.str.contains("2012")]

            print("Hydro inflow currently fixed to 2012!!")

            ds = pd.read_csv(data_path / flist.iloc[0])
            # ds = pd.read_csv(
            #     pathlib.Path(pathlib.PureWindowsPath(
            #         "C:\\science\\highRES\\europe\\data\\hydro\\ror"
            #         "\\ror_2013.csv")))

            if rleap:
                ds = pd.pivot_table(ds, columns="zones")
                ds.index = ds.index.astype(int)

                ds = ds.sort_index()
                ds.index = pd.date_range(
                    start=pd.Timestamp(2012, 1, 1),
                    end=pd.Timestamp(2012, 12, 31) + datetime.timedelta(hours=23),
                    freq="H",
                )

                noleap = (ds.index >= pd.Timestamp(2012, 2, 29, 0)) & (
                    ds.index <= pd.Timestamp(2012, 2, 29, 23)
                )

                ds = ds[~noleap]

                ds.index = np.arange(len(ds)).astype(str)

                ds = ds.T.reset_index()

            delta_h = str(int((dend - dstart).total_seconds() / 3600))

            ds = ds.loc[:, :delta_h]

            if agg:
                uk = ds.loc[ds["zones"].str.contains("UK"), :]
                ds = aggregate_zones(ds, agg, "mean")
                # ds=ds.loc[ds.index!="UK",:]
                ds = ds.reset_index()
                ds.rename({"EU_simpl": "zones"}, axis=1, inplace=True)

                ds = pd.concat((ds, uk))

            cells = ds["zones"].values
            zones = ds["zones"].values

            ds = ds.set_index("zones").stack().reset_index()

            ds = ds.rename(columns={"level_1": "h", 0: "Value"})
            ds.insert(0, "Tech", np.repeat(tech, ds.shape[0]))

            ds = ds[["h", "Tech", "zones", "Value"]]

            cf_out.append(ds)

            area_out.append(
                data2dd(
                    np.zeros_like(zones) + np.inf,
                    [np.array([tech]), zones, cells],
                    rounddp=0,
                )
            )

            cells_out.append(cells)
            zones_out.append(zones)

            continue

        cell_areas = pd.read_csv(cell_areas)

        try:
            cell_areas = cell_areas.rename({"Name_1": "zone"}, axis=1)
        except KeyError:
            pass
        try:
            cell_areas = cell_areas.rename({"zones": "zone"}, axis=1)
        except KeyError:
            pass

        if sel_zones is not None:
            cell_areas = cell_areas.loc[cell_areas["zone"].isin(sel_zones), :]

        # build multi year data set

        if dstart.year != dend.year:
            for i, f in enumerate(flist):
                if i == 0:
                    ds = xr.open_dataset(data_path / f)
                else:
                    temp = xr.open_dataset(data_path / f)
                    ds = xr.concat((ds, temp), dim="time")
        else:
            ds = xr.open_dataset(data_path / flist.iloc[0])

        # cut out desired timeslice

        ds = ds.sel(time=slice(dstart, dend))

        nlat = ds.lat.shape[0]
        nlon = ds.lon.shape[0]
        ntime = ds.time.shape[0]

        try:
            gen = ds["gen"].values
        except KeyError:
            gen = ds["capacity factor"].values
        gen[np.isnan(gen)] = 0.0

        #        gen=ds.gen
        #        gen.data[np.isnan(gen)]=0.
        #
        #        gen.data[0,:,:]=np.mean(gen.data,axis=0)
        #
        #        mask=np.mean(gen.data,axis=0)
        #
        #        gen.data[0,mask <cf_cutoff]=0.
        #
        #        gen[0,:,:].salem.quick_map()
        #
        #        plt.show()

        # reshape data array into ntime by (nlat*nlon)

        arr = gen.reshape(-1, nlat * nlon).T
        arr_id = np.arange(arr.shape[0])

        # remove low cf cells

        arr_id = arr_id[(np.mean(arr, axis=1) > cf_cutoff)]
        arr = arr[(np.mean(arr, axis=1) > cf_cutoff), :]

        gen = pd.DataFrame(arr, columns=np.arange(ntime))
        gen.index = arr_id
        gen.index.name = "cell_id"
        gen = gen.reset_index()

        #        print(gen.loc[gen["cell_id"]==33484,:])
        #
        #        print(gen.shape)

        nmerge = cell_areas.shape[0]

        gen = cell_areas.merge(gen, how="inner", on="cell_id")

        # check that the number of cells in gen is less than or equal to that
        # in cell_areas
        # -> every zone/cell pair in cell_areas must have a time series unless
        # than cell has a low CF

        if gen.shape[0] > nmerge:
            raise ValueError(
                "Too many rows in generation dataframe, check merge method"
            )

        gen = gen.loc[gen["cell_area"] > min_area, :]

        if zonal:
            # seems to be the fastest was to do a weighted mean -> .loc
            # reassignment
            # into a DF is slow

            cf = gen.loc[:, np.arange(ntime)].values

            areas = gen["cell_area"].values.reshape(-1, 1)
            zones = gen["zone"].values

            cf = cf * areas

            gen = pd.DataFrame(cf, columns=np.arange(ntime))

            gen.insert(0, "zone", zones)
            gen.insert(1, "cell_area", areas)

            groups = gen.groupby("zone")

            gen = (
                groups.sum()
                .loc[:, np.arange(ntime)]
                .div(groups["cell_area"].sum(), axis=0)
            )

            # old unweighted mean

            # gen=groups.mean()

            gen.insert(0, "cell_area", groups["cell_area"].sum())

            gen = gen.reset_index()

            cells = gen["zone"].values
            zones = gen["zone"].values

            if write_csv:
                csv = gen.copy()

                csv.insert(0, "Tech", np.repeat(tech, csv.shape[0]))

                csv.insert(2, "mean_cf", csv.loc[:, np.arange(ntime)].mean(axis=1))

                csv_out.append(csv)

            # it is faster to do the below to the df than to feed into data2dd

            cf = gen.copy()

            cf = cf.drop("cell_area", axis=1).set_index("zone").stack().reset_index()
            cf = cf.rename(columns={"level_1": "h", "Name_1": "r", 0: "Value"})
            cf.insert(0, "Tech", np.repeat(tech, cf.shape[0]))
            cf = cf[["h", "Tech", "zone", "Value"]]

            cf_out.append(cf)

        else:
            print(gen)

            stop

            cells = gen["cell_id"].values
            zones = gen["zones"].values

            cf = (
                gen.drop_duplicates("cell_id")
                .set_index("cell_id")
                .drop(["cell_area", "zone"], axis=1)
                .stack()
                .reset_index()
            )

            cf = cf.rename(columns={"level_1": "h", "cell_id": "r", 0: "Value"})

            cf.insert(0, "Tech", np.repeat(tech, cf.shape[0]))

            cf_out.append(cf)

        area_out.append(
            data2dd(gen["cell_area"], [np.array([tech]), zones, cells], rounddp=0)
        )

        cells_out.append(cells)

        zones_out.append(zones)

        del ds

    area_out = np.vstack(area_out)

    cells_out = np.hstack(cells_out)
    zones_out = np.hstack(zones_out)

    if dstart.year != dend.year:
        yr = str(dstart.year) + "-" + str(dend.year)
    else:
        yr = str(dstart.year)

    print(ntime)

    np.savetxt(out / (fmod + "_regions.dd"), np.unique(cells_out), fmt="%s")
    np.savetxt(out / "zones.dd", np.unique(zones_out), fmt="%s")

    cf_out = np.vstack(cf_out)
    # cf_out[:,-1]=np.round(cf_out[:,-1].astype(float),4)

    # cf_out=pd.DataFrame(cf_out,columns=["vre","r","h","Value"])
    cf_out = pd.DataFrame(cf_out, columns=["h", "vre", "r", "Value"])
    cf_out.loc[:, "Value"] = np.round(cf_out.loc[:, "Value"].astype(float), 3)

    # cf_out.loc[:,"Value"]=np.round(cf_out.loc[:,"Value"].astype(float),3)

    cf_file = out / ("vre_" + yr + "_" + fmod)
    area_file = out / ("vre_areas_" + yr + "_" + fmod)

    if zonal and write_csv:
        csv_out = pd.concat(csv_out)
        csv_out.to_csv(csv_file + "vre_gen_" + yr + "_" + fmod + ".csv", index=False)

        # cf_file=out/("vre_"+str(yr)+"_agg_"+fmod)
        # area_file=out/("vre_areas_"+str(yr)+"_agg_"+fmod)

    # else:

    # cf_file=out/("vre_"+str(yr)+"_"+fmod)
    # area_file=out/("vre_areas_"+str(yr)+"_"+fmod)

    # cf_out = (cf_out
    #           .loc[cf_out["r"].str.contains("UK") | cf_out["r"]
    #                .str.contains("BNL"),:])

    # area_out = pd.DataFrame(area_out,columns=["a","b"])

    # area_out = (area_out
    #             .loc[area_out["a"].str.contains("UK") | area_out["a"]
    #                  .str.contains("BNL"),:].values)

    cf_out.to_csv(cf_file.parent / (cf_file.name + ".csv"), index=False)

    wrapdd(area_out, "area", "parameter", outfile=area_file + ".dd")

    os.system(
        "csv2gdx "
        + str(cf_file.parent / (cf_file.name + ".csv"))
        + " output="
        + str(cf_file.parent / (cf_file.name + ".gdx"))
        + " ID=vre_gen Index=(1,2,3) Value=(4) UseHeader=True"
        "StoreZero=True"
    )
    (cf_file.parent / (cf_file.name + ".csv")).unlink()


def euro_add_roofpv(ddfile, roofareas):
    roofareas = pd.read_csv(roofareas)

    dd = pd.read_csv(ddfile, skiprows=1, delimiter=" ", skipfooter=2, engine="python")

    temp = dd.copy()

    temp[["tech", "zone1", "zone2"]] = temp["area"].str.split(".", expand=True)

    temp = temp.loc[temp["tech"] == "Solar", :].merge(
        roofareas, left_on="zone1", right_on="zone"
    )

    temp["/"] = temp["/"] + temp["useable roof surface (km2)"]

    for a in temp["area"]:
        dd.loc[dd["area"] == a, "/"] = temp.loc[temp["area"] == a, "/"]

    wrapdd(dd[["area", "/"]].values, "area", "parameter", outfile=ddfile)


def gb_add_roofpv(ddfile, roofareas):
    roofareas = pd.read_csv(roofareas)[["Name_1", "Sum_area"]]

    dd = pd.read_csv(ddfile, skiprows=1, delimiter=" ", skipfooter=2, engine="python")

    temp = dd.copy()

    temp[["tech", "zone1", "zone2"]] = temp["area"].str.split(".", expand=True)

    temp = temp.loc[temp["tech"] == "Solar", :].merge(
        roofareas, left_on="zone1", right_on="Name_1"
    )

    # assume 43% of roof area can achieve modelled azimuth (south) and tilt

    temp["/"] = temp["/"] + (temp["Sum_area"] * 0.43 / 1e6).round(0)

    for a in temp["area"]:
        dd.loc[dd["area"] == a, "/"] = temp.loc[temp["area"] == a, "/"]

    wrapdd(dd[["area", "/"]].values, "area", "parameter", outfile=ddfile)
