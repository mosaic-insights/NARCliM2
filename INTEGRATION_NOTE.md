# Integration note

This repository contains the NARCliM2 climate-data workflow together with
climate uncertainty, spatial-summary, rainfall and temperature time-series,
storm-hazard, and flood-hazard analysis tools.

Before pushing an updated version of the repository to GitHub:

1. Ensure the package modules are located under:

   `src/narclim_workflow/`

2. Confirm that `src/narclim_workflow/__init__.py` exposes the intended
   public API, including:

   - `run_workflow`
   - `run_uncertainty_workflow`
   - `run_spatial_summary`
   - `plot_ensemble_maps`
   - `plot_ensemble_raster`
   - `run_rain_temp_timeseries`
   - `run_flood_hazard`

3. Keep the core NARCliM workflow modules, including:

   - `catalog.py`
   - `config.py`
   - `domain.py`
   - `spatial.py`
   - `workflow.py`
   - `uncertainty.py`
   - `spatial_summary.py`
   - `ensemble_maps.py`

4. Include the additional analysis modules:

   - `rain_temp_timeseries.py`
   - `flood.py`

5. Ensure `config.py` contains the current generic `StormDays`
   `VariableSpec` used by the configurable storm workflow.

6. Ensure `pyproject.toml` contains all dependencies required by the
   package, including dependencies used by the rainfall/temperature
   time-series and flood-hazard modules.

7. Update the repository-level documentation where required:

   - `README.md`
   - `example.py`
   - `NARCliM_Example.ipynb`
   - `CHANGELOG.md`
   - `CITATION.cff`

8. Keep downloaded climate data and generated outputs out of the Git
   repository using `.gitignore`.

9. Before pushing to GitHub, run:

   - a small NARCliM download/subsetting test;
   - the uncertainty workflow;
   - the rainfall/temperature time-series workflow;
   - the flood-hazard workflow;
   - the example notebook.

10. Check that examples and documentation do not contain local paths,
    credentials, API keys, or project-specific confidential information.

11. Choose and include an appropriate `LICENSE` before making the
    repository public.