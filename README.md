# Scripts for analyses used in the _Stb15_ project
These scripts were used for analysis and generation of figures in the [publication detailing the cloning of Stb15 using the Watkins collection of wheat landraces](https://doi.org/10.1101/2023.09.11.557217) (Hafeez et al., 2023). 

R scripts/notebooks include the statistical analysis of Septoria tritici blotch pathology data, plot generation, haplotype analysis and Principal Coordinate Analysis. R markdown notebooks are also provided as `.html` files. These can be downloaded and viewed in a browser to allow the code alongisde outputs to be visualised.

Calculation of a SNP distance matrix was done using a Python script.

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

These steps were used to analyse the three initial screens carried out on the Watkins core 300 collection. Similar steps were undertaken in **Step 5**, wherein a more detailed notebook detailing the analysis process is provided.

## Step 2: GWAS
The estimated means calculated as described above were used as input for GWAS as part of the [WatSeq project](https://doi.org/10.1038/s41586-024-07682-9). 
[See GWAS pipeline here](https://github.com/ShifengCHENG-Laboratory/WWWG2B). 

## Step 3: Generation of SNP distance matrix 

Python script: `vcf_parse_distance.py` 

Written by Amber Hafeez and Burkhard Steuernagel

The candidate regions from the GWAS were analysed with the above script to later define haplotype groups. 
The script parses VCF and GFF files and calculates a SNP distance matrix within a specified region. The SNP distance matrix is the number of SNPs which differ between each pair of samples in the dataset. E.g.

```
                ArinaLrFor  ChineseSpring  Watkins1  Watkins2
ArinaLrFor      0           12             4         12
ChineseSpring   12          0              12        0  
Watkins1        4           12             0         12
Watkins2        12          0              12        0            
```

In the above example, Watkins2 may have the same allele as Chinese Spring (SNP distance = 0 between these lines), whilst Watkins1 differs from the other lines (SNP distance >0) so may have a unique allele.

Usage: `vcf_parse_distance.py gff_file.gff vcf_file.vcf output_file.tsv chromosome start_position end_position`

Further details about input files required:
- **gff_file.gff**: this is a GFF file for the reference genome.
- **vcf_file.vcf**: this is the VCF file containing SNPs called against the reference genome.
- **output_file.tsv**: this is the name of the desired output file. It is in TSV format.
- **chromosome**: this is the chromosome of the locus of interest. It must exactly match the chromosome name as specified in the GFF file.
- **start_position** and **end_position**: the start and end position of the locus of interest.


## Step 4: Determination of haplotype groups and their association with resistance

R Markdown notebook: `Haplotype_heatmaps.Rmd`

This notebook details the process of using a SNP distance matrix of a gene/locus (**Step 3**) and phenotype data (**Step 1**) to generate a heatmap using `pheatmap` and subsequently identify alleles of interest.

### Step 4.5: Generation of informative plots

R Markdown notebook: `Pycnidia_plots_Stb_genes.Rmd`

This script shows how the plots in figures 1b, 3a and 3b were generated. These figures all have the alleles of _Stb6_ and _Stb15_ indicated by the colour and shape of the points. 
- 1b: the pycnidia scores from assays of _Z. tritici_ isolates IPO323 and IPO88004 plotted against eachother to allow interpretation of isolate-specific resistance conferred by _Stb_ gene alleles.
- 3a: landraces plotted on a world map using the coordinates of where they were originally collected.
- 3b: a Principal Coordinate Analysis of Axiom array SNP data for the landraces.

The PCA plot was generated following the steps of this script for data preparation: https://rpubs.com/liziqi961/980723  
And this script for data analysis: https://rpubs.com/liziqi961/981419

## Step 5: Evaluation of transgenic lines carrying _Stb15_

R Markdown notebook: `Stb15_transgenics_analysis.Rmd`

This notebook contains a step-by-step explanation of the steps undertaken to fit linear mixed models to pathology results from testing _Z. tritici_ isolates on transgenic lines carrying _Stb15_ and controls. The process for generation of plots from Supplementary Figure 7 and their interpretation is also detailed. 

## Step 6: Evolutionary analysis

The protein predicted from the ArinaLrFor allele of Stb15 was used to perform blastp searches against 32 plant species and the top hits were used to generate a phylogenetic tree. The clade containing Stb15 was then extracted from the tree and the gene annotations corresponding to these proteins were obtained from the GFF files of their respective genomes. 

The intron/exon structure of each gene in the clade most closely related to _Stb15_ was plotted in R. The code used to do this is presented in `draw_gene_annotations.R`. The script was written with help from [this tutorial](https://dbsloan.github.io/TS2019/exercises/r_figure_drawing.html).

The main input file required is a text file containing gene annotations extracted from GFF files - the gene name, start position, end position, strand and feature (exon, CDS, gene). Additionally, information about the clade position of each gene can be used to plot each gene in the correct order so it can later be aligned with a tree. 

This method was used to generate Figure 3c and Supplementary Figure 8 in Hafeez et al. (2023).
