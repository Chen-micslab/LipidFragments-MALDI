# R scripts for processing and plotting mass spectra

## Introduction

This is a collection of R scripts used for processing mass spectra generated along the research **"Uncover the Interference of Lipid Fragments on the Qualification and Quantification of Serum Metabolites in MALDI-TOF MS Analysis"**.

Four directories are included in this collection: 

- "lib": collections of utility functions that are used by other scripts (adapted from the [RealMass](https://github.com/zwpwjwtz/RealMass) project)

- "process": scripts for processing (e.g. smoothing and peak picking) raw mass spectra

- "summary": scripts that summarize (e.g. cluster and average) processed spectra and peak lists, and generate averaged peak list among replicates of spectra

- "plot": scripts for visualize spectra and peaks

## How to use

### Paths and working directory

Since the scripts use relative path to reference each other, it is advised to keep the structure of directory unchanged.

We also assumed that the current working directory (check it with ```getwd()``` in R) is the root directory of this collection. It is then advised to set the working directory (using ```setwd()```) to the root of this collection before run any script.

### Source data

Source data can be found on [iProx](http://www.iprox.org/page/project.html?id=IPX0003892000). Please extract the downloaded archives into a subdirectory named "spectra", under a directory named "data" alongside the project root, so that processing scripts will find them. Two additional subdirectories under the "data" directory, named "summary" and "plot", will be necessary if further statistics and plotting are required. The final directory structure is expected as follows:

```
-- (upper level directory)
	|-- data
	|	|-- plot
	|	|-- spectra
	|	|	|-- MALDI-TOF
	|	|	|	|-- 2021-09-10
	|	|	|	|-- (other subdirectories)
	|	|	|-- HPLC-MS
	|	|		|-- 20210917
	|	|		|-- (other subdirectories)
	|	|-- summary
	|		|-- MALDI-TOF
	|		|-- HPLC-MS
	|-- scripts (the project root)
		|-- lib
		|-- plot
		|-- process
		|-- summary
```

Alternatively, you may use the [init](./init.sh) script provided in this project to create the above directories for you.

###  Dependencies

Following R packages (and their versions) are required:

- alsace (1.14.0)
- ggplot2 (3.1.0)
- gridExtra (3.4.4)
- MALDIquant (1.19.3)
- mzR (3.5.2)

Packages of other versions might also work, but without guarantee.

### Have an overview

A demo script that demonstrate the basic workflow of metabolite quantification can be found at [main.R](./main.R). You may also choose to run scripts dedicating to a specific step, but be sure that the input data of that step has been generated before (e.g. peak lists must be generated before metabolite quantification can be made).

### Process mass spectra

To process and process a series of mass spectra, run 
```
Rscript ./process/processMZ_XXX.R
```
in a terminal, or
```
source('./process/processMZ_XXX.R')
```
in a R console (or the console in RStudio), where "XXX" is the suffix of a specific script file (usually indicating the date when the spectra were first processed).

### Preprocess (2-D) mass spectra

Data generated by a LC-MS system are two-dimensional spectra. To focus on the analyte of interest and to simplify the workflow, one may first need to reduce them to 1-D spectra. The idea here is relatively simple: for each target analyte, choose the time segment where it is eluted, then average all sub-spectra ("scans") acquired within this time segment. The peak intensity (or area) in the averaged spectrum will be an approximation to that in an integral spectrum (given that time spans between scans are not diverse too much).

To preprocess a series of 2-D mass spectra, run 
```
Rscript ./process/processLCMZ_XXX.R
```
in a terminal, or
```
source('./process/processLCMZ_XXX.R')
```
in a R console (or the console in RStudio), where "XXX" is the suffix of a specific script file.

The preprocessed 2-D mass spectra is then 1-D, and can serve as the input of processing scripts (the "processMZ_" series).

### Summary mass spectra

To summarize a series of mass spectra, run 
```
Rscript ./process/Summary-CommonMZ_XXX.R
```
in a terminal, or
```
source('./process/Summary-CommonMZ_XXX.R')
```
in a R console (or the console in RStudio), where "XXX" is the suffix of a specific script file.

### Plot mass spectra

To plot a series of processed mass spectra, run 
```
Rscript ./plot/Plot-Intensity_XXX.R
```
in a terminal, or
```
source('./plot/Plot-Intensity_XXX.R')
```
in a R console (or the console in RStudio), where "XXX" is the suffix of a specific script file.

## License

All files included in this collection are licensed under GNU General Public License Version 3 (GPLv3). The licence text can be found in file [LICENSE](./LICENSE).
