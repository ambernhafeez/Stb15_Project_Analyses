# Scripts for analyses used in the _Stb15_ project
R scripts used for the analysis of Septoria tritici blotch pathology data, plot generation, haplotype analysis adn Principle Coordinates Analysis for the [publication detailing the cloning of_Stb15_using the Watkins collection of wheat landraces](https://doi.org/10.1101/2023.09.11.557217). 

The steps each script/notebook was used for are detailed below. 

## Step 1: Analysis of Septoria pathology data in the core 300 Watkins landraces
R script: `Watkins_core300_Septoria_analysis.R`
This script contains several steps as follows:
1. Read and format the data
2. Calculate the % maximum dAUDPC
3. Fit linear mixed models to pycnidia and damage datasets
4. Generate residuals plots
5. Calculate estimated marginal means from the models
6. Final formatting of datasets

These steps were used to analyse the three initial screens carried out on the Watkins core 300 collection. 

## Step 2: GWAS
The estimated means calculated as described above were used as input for GWAS as part of the [WatSeq project](https://doi.org/10.1038/s41586-024-07682-9). 
[See GWAS pipeline here](https://github.com/ShifengCHENG-Laboratory/WWWG2B). 

## Step 3: Generation of SNP distance matrix 
Python script: `VCFparse_distance` 
Written by Amber Hafeez and Burkhard Steuernagel

The candidate regions from the GWAS were analysed with the above script to later define haplotype groups. 
The script parses VCF and GFF files and calculates a SNP distance matrix within a specified region.

Input files required:

## Step 4: Determination of haplotype groups and their association with resistance
R Markdown notebook:

## Step 5: Generation of informative plots
R Markdown notebook:

## Step 6: Evaluation of transgenic lines carrying _Stb15_
R Markdown notebook:

