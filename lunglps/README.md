# README

Lunglps analysis run -> 

## Description


## Usage

We assume that you have already cloned the repository and that some version of python3 is setup to work with pip. 

1. Install the required software (snakemake and conda).:
```sh
pip install -r requirements.txt
```

2. Download data with provided bash script make sure your working directory is in the lunglps folder
```sh
sh download_data.sh
```

3. Execute the workflow, from the lunglps folder. 
```sh
snakemake --cores 4 --configfile lung_lps.yaml --snakefile rules/process_10x.smk
```

This analysis can take a while to finish if on a low powered machine. A larger amount of RAM is recommend else the analysis may fail. 

## Results

The results folder all results. The main part is a html file in results/merged/report/report.html 
which contains a variety of figures derived from this analysis.  