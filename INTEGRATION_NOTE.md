# Integration note

This bundle contains the repository-level files revised from the supplied
README, example, notebook, and pyproject configuration, plus the latest
changed workflow modules produced in this conversation.

Before pushing your existing repository to GitHub:

1. Copy `README.md`, `example.py`, `NARCliM_Example.ipynb`, `pyproject.toml`,
   `.gitignore`, `CITATION.cff`, and `CHANGELOG.md` to the repository root.
2. Copy `__init__.py` to `src/narclim_workflow/__init__.py` after checking
   that it does not remove any additional public exports you already use.
3. Copy files under `latest_changed_modules/` into
   `src/narclim_workflow/`.
4. Keep your existing `catalog.py`, `config.py`, `domain.py`, `spatial.py`,
   and any other package modules not included in this bundle.
5. Ensure `config.py` contains the current generic `StormDays` VariableSpec
   used by the configurable storm workflow.
6. Run your example notebook and a small one-model test before a public push.
7. Choose an appropriate LICENSE before making the repository public.
