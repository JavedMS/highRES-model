#!/usr/bin/env bash

# TODO: export environment variable to then have gams read the path from there.
# TODO: also call gams from here and then set the log file path here too


cd ${snakemake_params[modelpath]}
pwd -P
${snakemake_params[gamspath]}gams highres2.gms logOption=2 gdxCompress=1 \
--hydro_res_min "${snakemake_params[hydroresmin]}"