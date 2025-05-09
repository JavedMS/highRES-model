# Engaging young people in modelling for a more inclusive energy transition at the country scale.

## What is highRES-Norway?

Welcome to the repository for the Norwegian version of the high temporal and spatial resolution electricity system model (highRES-Norway). The highRES-Norway is specifically modified from highRES European version (<https://github.com/highRES-model>) to assess the impact of incorporating the stakeholders (which is youth in our study) impact on achieving net-zero electiricty system for Norway. The preferences and perspectives of young people are incorporated by introducing wildcards in the sanekmake workflow (trans, import_xxx, varnewpcapQ, corine_solar, corine_onshore, fylke_tech_limit) and the new variables in the core electicity system model developed in GAMS. These new variables are introduced in GAMS codes via snakemake rules (build_gams, add_transmission_type, import_export_changes). 

The pupil preferences related to electricity trading, transmission expansion, landscape choices, renewable technology choices, and regional choices for new renewable energy installations are incorporated by calculating the preference coefficients, reflecting the percent of participants choose the particular option. The detailed quantitative assessment of workshops data along with methodological details is available at <https://github.com/JavedMS/Mapping_Workshops_EFF>.  

Moreover, taking the benefit of snakemake scenario generation ability via wildcards, a range of scenarios explored based on the prioritization and incremental integration of pupil choices alongside step-wise exclusions of disagreed landscapes. We find that, given pupil priorities regarding certain power system elements and their cumulative impact, substantial shifts occur in national capacity potentials (approximately ±50%), cost projections (-7% to +25%), capacity mixes (notably from 40% to 0% onshore wind), and regional equity assessments—where high costs and youth-driven pathways do not necessarily guarantee equitable systems. Although applied to young people in Norway, the proposed workshop-informed modelling framework serves as a tool to meaningfully engage and empower diverse groups while understanding and incorporating their localised preferences in energy system planning.

highRES model is written in GAMS and its objective is to minimise power system investment and operational costs to meet hourly demand, subject to a number of system constraints. The transmission grid is represented using a linear transport model. To realistically model variable renewable supply, the model uses spatially and temporally-detailed renewable generation time series that are based on weather data. The further documentation details about mathematical formulation, nomenclatures, abbreviations, setting-up the configuration file, and snakemake workflow rules can be found here: <https://highres-europe-wf.readthedocs.io/en/latest/introduction.html>.  

## How to run the model

GAMS must be installed and licensed. This version was tested/developed with GAMS version 27.2.0. To run the full workflow, two datapackages are needed they can be downloaded from:

1. (~80MB compressed, ~300MB uncompressed) <https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/Eftsg10mEK9Mpi4TSN8aS9kBWlooGJ_99YDDaYcGiQvrYQ?e=xeu9Lk&download=1>.
2. (~10GB) <https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/EdEmFkUQoL5Imy3-OumK_o0BcFqilpjB3CQOCbUwi_1T8g?e=O0kq50&download=1>

## Windows
1. Clone the repository
2. Install snakemake
    - Download miniforge windows exe <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe>
    - Install Minforge
    - Run the minimal install of the snakemake environment `mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal pandas zstd`
3. Activate the snakemake environment
4. Navigate to the repository in your snakemake conda environment shell
4. Get the required input files
    ```
   curl -o shared_input.tar -L -b cookies.txt "https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/EdEmFkUQoL5Imy3-OumK_o0BcFqilpjB3CQOCbUwi_1T8g?e=O0kq50&download=1" -o resources.tar.zst -L -b cookies.txt "https://uio-my.sharepoint.com/:u:/g/personal/tobiasvh_uio_no/Eftsg10mEK9Mpi4TSN8aS9kBWlooGJ_99YDDaYcGiQvrYQ?e=xeu9Lk&download=1"
   ```
5. Extract the required input files
    ```
    zstd -d resources.tar.zst
    mkdir resources
    tar xf resources.tar -C resources
    mkdir shared_input
    tar xf shared_input.tar -C shared_input
    ```
5. Run snakemake -c --use-conda
 

## References

This model uses the following data sources:

- ERA5 for solar, on and offshore wind capacity factors and runoff data (<https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>)

- Demand data is from ENTSO-E for the year 2010

- Cost and technical data is taken from UKTM (<https://www.ucl.ac.uk/energy-models/models/uk-times>) with data from the JRC report "Cost development of low carbon energy technologies: Scenario-based cost trajectories to 2050" (<https://publications.jrc.ec.europa.eu/repository/handle/JRC109894>) used to update some areas.

- Data on run-of-river, reservoir and pumped hydro power capacities is taken from <https://transparency.entsoe.eu/>, <https://www.entsoe.eu/data/power-stats/> and <https://github.com/energy-modelling-toolkit/hydro-power-database>

- Energy storage capacities for reservoir and pumped storage are taken from Schlachtberger et al. (2017) and Geth et al. (2015) respectively (see <https://doi.org/10.1016/j.energy.2017.06.004> and <https://doi.org/10.1016/j.rser.2015.07.145>).


The initial version of highRES electricity model was published in:

- Price, James, and Marianne Zeyringer. 2022. ‘HighRES-Europe: The High Spatial and Temporal Resolution Electricity System Model for Europe’. SoftwareX 17 (January): 101003. <https://doi.org/10.1016/j.softx.2022.101003>.

