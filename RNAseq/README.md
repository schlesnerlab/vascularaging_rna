# README

## Description

This repository contains the analysis code for the paper
titled "Vascular Aging causes Metastasis".
This repository includes analysis code for bulk
RNAseq analysis and single cell Analysis of Tabula Muris Senis
and LPS treated mice.

## Structure

Each folder holds the analysis to reproduce the analysis and figures.
Note some analysis require larger amounts of computational power and RAM to
be successfully executed.

## Usage

To run the analysis, follow these steps:

1. Clone the repository:

    ```sh
    git clone https://github.com/schlesnerlab/vascularaging_rna.git
    cd vascularaging_rna/RNAseq
    ```

2. install Required software 

    ```sh
    pip install -r requirements.txt
    ```

3. Run snakedeploy to deploy workflow:
```sh
snakedeploy deploy-workflow https://github.com/schlesnerlab/multicondition-deseq2-enrichment multicondition-deseq2-enrichment --branch paper_freeze 
```

4. Download Data:

    ```sh
    cd data
    bash download_data.sh
    ```

5. Execute Workflow
    
    ```sh
    cd multicondition-deseq2-enrichment
    snakemake --cores 4 --configfile config/VascAge.yaml
    ```
## Data

The data used in this analysis are the following:

- RNAseq: - GSE********
- Tabula Muris Senis: "<https://figshare.com/ndownloader/files/15467792>"
- LungLPS: GSE14849

Exact instrutions are given in the README

Divided into the following steps:

## Results

The results of the analysis, including tables and figures, are saved in the `results/` directory.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
