# Command line scripts for performing RDKit-based cheminformatics operations

Questions? email danrusso@rutgers.edu

## Requirements

The scripts used in this package require Python and Python libraries. Fortunately, ([Anaconda](https://www.anaconda.com/download/)
provides a version of Python and an easy way to load all the required libraries for this project  using the `environment.yml` file.
Please follow the directions on the installer to set up Anaconda.

## Getting set up

1. After installing Anaconda a python environment with all the required packages can be created using the `environment.yml`
file by running the following command in the directory containing the file:

`$ conda env create -f environment.yml` in this directory.

2) After successful installation of the python environment it can be loaded on by running the following command:

`$ activate rdkit` for windows or `$ source activate rdkit` for linux/mac

3) With the rdkitcmd environment loaded, all scripts can be run as directed in the script from the command line.

