# Engaging young people in modelling for a more inclusive energy transition at the country scale.

## What is highRES-Norway?

Welcome to the repository for the Norwegian version of the high temporal and spatial resolution electricity system model (highRES-Norway). The highRES-Norway is specifically modified from [highRES European version](https://github.com/highRES-model) to assess the impact of incorporating the stakeholders (which is youth in our study) on acheiving net-zero electiricty system for Norway. The preferences and perspectives of young people are incorporated by introducing wildcards in the sanekmake workflow (``trans``, ``import_xxx``, ``varnewpcapQ``, ``corine_solar``, ``corine_onshore``, ``fylke_tech_limit``) and the new variables in the core electicity system model developed in ``GAMS``. These new variables are introduced in ``GAMS`` codes via snakemake rules (``build_gams``, ``add_transmission_type``, ``import_export_changes``). 

The pupil preferences related to ``electricity trading``, ``transmission expansion``, ``landscapes``, ``renewable technology``, and ``regions`` for net-zero power system are incorporated by calculating the preference coefficients, reflecting the percent of participants choose the particular option. The detailed quantitative assessment of workshops data along with methodological details is openly available at [github](https://github.com/JavedMS/Mapping_Workshops_EFF).  

Moreover, taking the benefit of ``snakemake`` scenario generation ability via wildcards, a range of scenarios explored based on the prioritization and incremental integration of pupil choices alongside step-wise exclusions of disagreed landscapes. We find that, given pupil priorities regarding certain power system elements and their cumulative impact, substantial shifts occur in national capacity potentials (approximately ±50%), cost projections (-7% to +25%), capacity mixes (notably from 40% to 0% onshore wind), and regional equity assessments—where high costs and youth-driven pathways do not necessarily guarantee equitable systems. Although applied to young people in Norway, the proposed workshop-informed modelling framework serves as a tool to meaningfully engage and empower diverse groups while understanding and incorporating their localised preferences in energy system planning.

``highRES`` core model is written in ``GAMS`` and its objective is to minimise power system investment and operational costs to meet hourly demand, subject to a number of system constraints. The transmission grid is represented using a linear transport model. To realistically model variable renewable supply, the model uses spatially and temporally-detailed renewable generation time series that are based on weather data. The further documentation details about mathematical formulation, nomenclatures, abbreviations, setting-up the configuration file, and ``snakemake`` workflow rules can be found [here](https://highres-europe-wf.readthedocs.io/en/latest/introduction.html).  

## How to run the model

This repositry is structure around ``snakemake`` workflow; dependencies are managed using conda/mamba. ``GAMS`` must be installed and licensed. This version was tested/developed with GAMS version 27.2.0. To run the full workflow, following datapackages are needed to download:

1. Resources (488 MB compressed, 843 MB uncompressed) <https://zenodo.org/records/15401853/files/resources.zip?download=1>

2. Shared input (10.4 GB compressed, 10.8 GB uncompressed) <https://zenodo.org/records/15401853/files/shared_input.zip?download=1>

3. wind_bias_correction (23.7 GB compressed, 35.5 GB uncompressed)

## Windows

1. Clone the repository
2. Install snakemake
    - Download miniforge windows exe <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe>
    - Install Minforge
    - Run the minimal install of the snakemake environment `mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal pandas zstd`
3. Install and activate the ``snakemake_env.yml`` environment
4. Navigate to the repository in your snakemake conda environment shell
5. Get the required input files
    ```
   curl -L -b cookies.txt "https://zenodo.org/records/15401853/files/resources.zip?download=1" --output resources.zip
   curl -L -b cookies.txt "https://zenodo.org/records/15401853/files/shared_input.zip?download=1" --output shared_input.zip
   curl -L -b cookies.txt "https://zenodo.org/records/15401853/files/wind_bias_correction.zip?download=1" --output wind_bias_correction.zip   
 
   ```
6. Extract the required input files
    ```
    unzip resources.zip
    unzip shared_input.zip
    unzip wind_bias_correction.zip

    ```
7. set the paths for GAMS, input data, and results directory in config file ``cluster_system_os_config.yaml``  
8. Run snakemake -c --use-conda
   - Running snakemake will create and use the environments defined in ``highRES-Norway/workflow/envs``, which contains the specific version numbers of each dependency that are known to work.

The installation should not take more than a few minutes on a typical laptop or desktop computer. 
 
## Computational requirements

Running the complete workflow in this repository involves executing all possible combinations (i.e., scenarios) of wildcards specified in ``rule all``. A single execution of the entire workflow typically takes approximately 45–50 minutes, encompassing ``GAMS building``, ``technoeconomic_inputs``, ``land exclusions``, ``building_inputs``, and ``run_model``. With parallelization, a total of 960 scenarios were efficiently conducted on the University of Oslo cluster in under four days. The main computational demand lies with rule ``build_weather``, where exclusions occur at the grid cell level.

Depending on your available resources, you can modify the CPLEX solver options found in ``highRES-Norway/resources/cplex.opt``, which may affect the model solving time. To test the workflow with a reduced number of scenarios or specific scenarios, adjust by commenting or uncommenting wildcard values specified at the workflow's start. For instance, to run a single scenario, retain only one value for each wildcard in ``rule all``.

## Orgainsation

The ``highRES-Norway`` model introduces several important rules, enhancing and expanding upon the existing [highRES European version](https://github.com/highRES-model). Below is a summary of these rules, including newly added functionalities:

1. The rule ``build_gams`` integrates new variables and equations into the default ``GAMS`` model. The associated script ``build_gams.py`` is responsible for adding these equations and variables to specific model files. 

2. The rule ``add_transmission_type`` sets the upper limit of the overhead and subsurace transmission capacities based on the stakeholder chocies.

3. The rule ``import_export_changes`` manages the import and export limits according to stakeholder preferences using script ``import_input_change.py``, ensuring the model reflects desired trade scenarios accurately 

4. The rule ``build_weather`` perform technical, environmental, and stakeholders-led land exclusions using the CORINE and high resolution Norwegian landscape data.  Because exclusions occur at the grid cell level, this rule requires most computational time and resources.

5. The rule ``build_inputs`` constructs and finalizes the input data necessary for the core ``GAMS`` model, setting the foundation for the model’s operation.

6. The rule ``run_model`` run the model, simulating stakeholder-influenced scenarios to evaluate outcomes based on varied parameters and conditions.
   


## References

The initial version of highRES electricity model was published in:

- Price, James, and Marianne Zeyringer. 2022. ‘HighRES-Europe: The High Spatial and Temporal Resolution Electricity System Model for Europe’. SoftwareX 17 (January): 101003. <https://doi.org/10.1016/j.softx.2022.101003>.

