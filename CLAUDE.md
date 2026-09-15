# anatembea project context

## What this project is

`anatembea` is an R package for inferring temporal malaria transmission patterns from monthly prevalence data, especially data collected at first antenatal care visits (ANC1). The package combines a mechanistic malaria transmission model with particle Markov chain Monte Carlo (pMCMC) to reconstruct plausible seasonal transmission dynamics, EIR, clinical incidence, and related burden metrics from observed prevalence.

The package is intended for settings where routine surveillance is sparse or indirect, and the goal is to translate cross-sectional or monthly prevalence observations into a time-varying transmission signal. The package supports both general-population prevalence fitting and ANC-based comparisons using pregnant women stratified by gravidity (for example primigravidae and multigravidae).

## High-level purpose

This codebase implements a semi-stochastic version of the malaria transmission model used elsewhere in the `malariasimulation` ecosystem, but with a stochastic mosquito emergence/random-walk process to allow better fitting to monthly prevalence time series.

In practical terms, the package:

- loads and cleans prevalence data
- converts dates into the internal time scale used by the particle filter
- sets up mechanistic model parameters and seasonality assumptions
- estimates an initial equilibrium state from observed prevalence or user-specified target prevalence
- fits the model to monthly data using pMCMC
- produces posterior summaries of transmission-related outputs such as EIR, incidence, and prevalence

## Main package workflow

The central pipeline is centered on `run_pmcmc()`, with the following broad flow:

1. Data preparation
   - `data_process()` or `data_process_annual()` converts raw monthly or annual data into a time-indexed format for `mcstate::particle_filter_data()`.
   - `format_na()` ensures missing or zero-tested observations do not produce invalid likelihood contributions.
   - `admin_match()` can choose a seasonal profile based on country/admin unit if a seasonality model is used.

2. Model parameter setup
   - `model_param_list_create()` creates the mechanistic parameter list used by the odin model.
   - This includes malaria biology parameters, immunity-related parameters, entomological parameters, and seasonality configuration.
   - The function also computes derived quantities such as initial equilibrium values and age-specific comparison parameters.

3. Initial state creation
   - `initialise()` estimates the starting state of the dynamical system from either a user-specified EIR or a target prevalence.
   - `user_informed()` and `data_informed()` wrap transformation logic for the pMCMC parameter set.
   - `get_init_EIR()` and related helpers can convert prevalence targets into initial transmission estimates.

4. pMCMC fitting
   - `run_pmcmc()` sets up the stochastic model, particle filter, proposal matrix, and posterior simulation.
   - It merges ANC datasets when comparing primigravidae and multigravidae, or fits to a non-pregnant general population when `comparison = 'u5'` or a flex age window like `'2to10'`.
   - It chooses the comparison function (`compare_u5`, `compare_pg`, `compare_mg`, `compare_pgmg`, `compare_ancall`, etc.) and the state index for the pMCMC outputs.
   - It calls the underlying `odin`/`dust` stochastic model and runs `mcstate` pMCMC.

5. Output handling
   - The returned object contains posterior samples and/or trajectory histories, depending on `save_state` and `save_trajectories`.
   - `pmcmc_trajectories()` extracts the time and state histories when they are available.
   - Output is typically summarized over time for key quantities: prevalence (`prev_*`), EIR, `betaa`, and clinical incidence.

## Important data expectations

The expected user data is a data frame with at least:

- `month` as a `yearmon` value from the `zoo` package
- `tested`: number tested
- `positive`: number positive

For ANC-specific analyses, the package may also use columns such as:

- `positive.pg`, `tested.pg` for primigravidae
- `positive.mg`, `tested.mg` for multigravidae
- `positive`, `tested` for general-population comparisons

The package supports comparison modes including:

- `u5`: under-5 prevalence
- `pg`: primigravid prevalence
- `mg`: multigravid prevalence
- `pgmg`: combined primigravid and multigravid comparison
- `pgsg`: primigravid and secondary gravidae style comparison
- `ancall`: all ANC data comparison
- age-range strings such as `2to10` for general population age windows

## Core files and responsibilities

### R/

- `R/run_pmcmc.R` — main pMCMC function and output-index logic
- `R/run_pmcmc_annual.R` — annualized alternative fitting workflow
- `R/model_parameters.R` — mechanistic parameter list and model setup
- `R/utils.R` — data cleaning, time processing, model comparison helpers, transformations, and probability conversions
- `R/equilibrium-init-create-stripped.R` — equilibrium initialization logic for the stripped model

### Package metadata

- `DESCRIPTION` — package name, R dependencies, imported packages, and test/vignette setup
- `NAMESPACE` — exported package functions
- `README.md` — high-level overview and installation instructions
- `vignettes/` — usage guides and examples
- `tests/testthat/` — project tests for data processing, likelihood comparisons, and parameter setup

## Model and inference approach

The package uses a stochastic random-walk mosquito emergence process rather than a purely deterministic seasonal model for monthly data. The random-walk process is embedded in the mechanistic disease model and then inferred by comparing simulated prevalence trajectories to observed prevalence using pMCMC.

The key inferred quantities are usually:

- prevalence in different populations/age groups
- mosquito emergence rate or effective transmission intensity (`betaa`)
- EIR over time
- clinical incidence (`clininc_*`)
- state variables and model diagnostics for internal checks

## Typical usage pattern

A typical call looks like this:

```r
library(anatembea)

result <- anatembea::run_pmcmc(
  data_raw = tanga_data_slim,
  target_prev = 0.4,
  n_particles = 200,
  n_steps = 1000,
  comparison = "u5"
)
```

This call:

- prepares the dataset
- estimates the initial state
- sets up the pMCMC and particle filter
- returns posterior samples and optionally trajectory history

## Important helper functions

A few helper functions are central to package behavior:

- `load_file()` — reads package resources from `inst/extdata`
- `match_clean()` — fuzzy string matching for country/admin unit lookup
- `admin_match()` — matches admin units and country for seasonality profile selection
- `format_na()` — fills missing values with zero for binomial likelihood calculations
- `data_process()` / `data_process_annual()` — transforms observations into particle filter input
- `transform_init()` — transforms equilibrium states for fitting
- `get_prev_from_log_odds()` and `get_odds_from_prev()` — convert odds and prevalence
- `get_anc_from_u5()` — maps under-5 prevalence to ANC strata prevalence
- `check_seasonality()` — validates seasonality profile assumptions

## Output conventions

`run_pmcmc()` returns a list-like object containing posterior and simulation outputs. Typical entries include:

- `mcmc`: MCMC samples
- `history`: posterior trajectory summaries for model states over time
- `trajectories`: optional saved particle-filter trajectories
- `output`: plotting or summary artefacts created by downstream analysis

The `output_level` argument controls how much state information is retained:

- `minimal`: only the values needed for the likelihood calculation
- `standard`: key fitted prevalence/EIR/betaa/incidence outputs
- `diagnostic`: full state index retained for debugging and visual diagnostics

## Testing and validation

The package includes `testthat` tests under `tests/testthat/` covering:

- data formatting and missing-data handling
- comparison likelihood functions
- odss/prevalence transformations
- admin-match fallback behavior
- pmCMC index creation and trajectory extraction
- model-parameter construction and derived values

The standard verification approach is:

```r
devtools::test()
```

or a more focused run through the package test folder.

## Repository notes for agents

- This is an R package, not a script repository. Prefer working within `R/` and keeping functions exported through `NAMESPACE` when public API changes are intended.
- The project is strongly dependent on the `dust`, `odin`, and `mcstate` toolchain; any model changes should preserve compatibility with those packages.
- The package uses roxygen documentation style. If adding exported functions or changing signatures, update the roxygen comments and regenerate documentation if needed.
- The vignettes are the best source of concrete usage examples and expected workflows.
- The central routine is `run_pmcmc()`, but `run_pmcmc_annual()` is the annual variant that uses a yearly time scale.

## Summary

This repository implements a malaria transmission reconstruction framework that uses monthly prevalence data to fit a mechanistic transmission model with pMCMC. It is designed to recover seasonal transmission intensity, burden, and the hidden time-varying drivers of infections from noisy observational prevalence data, particularly in ANC-based malaria settings.
