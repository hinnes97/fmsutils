# fmsutils
Diagnostic tools for ExoFMS output

# Installation guide:

- Clone this repository (for the tutorial to work you will need access to shared storage/a generic ExoFMS output file, so keep this in mind)
- Move into the top level directory (with setup.py and fmuenv.yml)
- Create a conda environment with the miniumum necessary prequisites using `conda env create --file fmuenv.yml` (these can also be installed by just manually installing xarray, cartopy, jupyter and windspharm (using the conda-forge channel) to a fresh environment)
- Activate this environment using `conda activate fmuenv`
- In the same directory, run `pip install -e .` which installs the package in 'editable' mode, i.e. you can change the source files and changes should take effect immediately
- Look at the notebook `tutorial.ipynb` for examples of how to use the package
