#!/bin/bash

# Set the GSE accession number
GSE_ACC="GSE148499"

# Create a directory to store the downloaded data
mkdir -p ./data

# Change to the directory
cd ./data

# Download the data using wget
wget -r -nH --cut-dirs=3 --no-parent --reject="index.html*" ftp://ftp.ncbi.nlm.nih.gov/geo/series/"${GSE_ACC:0:6}nnn/${GSE_ACC}/suppl/"

# Optional: Unzip the downloaded files if they are compressed
# Uncomment the following line if needed
# unzip "*.zip"
tar -xvf ./GSE148499/suppl/GSE148499_RAW.tar 

# Print a success message
echo "Data downloaded successfully!"