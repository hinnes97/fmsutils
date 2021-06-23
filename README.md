# fmsutils
Diagnostic tools for ExoFMS output

# Installation guide:

- Clone this repository (for the tutorial to work you will need access to shared storage/a generic ExoFMS output file, so keep this in mind)
- Move into the top level directory (with setup.py and fmuenv.yml)
- Create a conda environment with the miniumum necessary prequisites using `conda env create --file fmuenv.yml` (these can also be installed by just manually installing xarray, cartopy and jupyter to a fresh environment)
- Move into the top level directory (with setup.py and fmsutils.txt)
- Create a conda environment with the miniumum necessary prequisites using `conda env create --file fmuenv.txt -n fmuenv` (these can also be installed by just manually installing xarray, cartopy, windspharm (requires `-c conda-forge` channel) and jupyter to a fresh environment)
- In the same directory, run `pip install -e .` which installs the package in 'editable' mode, i.e. you can change the source files and changes should take effect immediately
- Look at the notebook `tutorial.ipynb` for examples of how to use the package
