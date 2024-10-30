# README

## Description

This repository contains the analysis code for the reproducing the cacoa analysis
on the tabula muris dataset using snakemake

## Usage

We assume that the user has already cloned the repository, has some version of python3 installed with pip and navigated to this folder

1. Install Required software 

    ```sh
    pip install -r requirements.txt
    ```
2. Execute workflow (this will download the data and run the analysis). Note these steps can you excess of 32 GB of RAM so it may be required to run analysis on high performance compute hardware with sufficient RAM

```sh
    snakemake --cores 4 --configfile tabulamuris.yaml --use-conda 
```

## Results

results can be found in the data/tabulamurs/report folder structure
