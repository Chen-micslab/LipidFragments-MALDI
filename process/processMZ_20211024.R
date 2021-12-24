# R script for pre-processing 1-D MS spectra

source('lib/PathUtils.R')
source('lib/MassSpectraProcessor.R')


mzDataDir <- '../data/spectra/MALDI-TOF'

mzOrigFileName <- c(
    '2021-10-24/25DHB+d3-Carnitine+SPCS/0_J1/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+SPCS/0_J2/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+SPCS/0_J3/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+SPCS/0_J4/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+SPCS/0_J5/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(5uM)/0_O1/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(5uM)/0_O2/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(5uM)/0_O3/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(5uM)/0_O4/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(5uM)/0_O5/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(10uM)/0_K1/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(10uM)/0_K2/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(10uM)/0_K3/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(10uM)/0_K4/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(10uM)/0_K5/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(20uM)/0_L1/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(20uM)/0_L2/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(20uM)/0_L3/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(20uM)/0_L4/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(20uM)/0_L5/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(30uM)/0_M1/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(30uM)/0_M2/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(30uM)/0_M3/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(30uM)/0_M4/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d3-Carnitine+Carnitine(30uM)/0_M5/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d9-Choline+HomoPhe+SerumExtract-100%EtOH-20210908/0_K7/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d9-Choline+HomoPhe+SerumExtract-100%EtOH-20210908/0_K8/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d9-Choline+HomoPhe+SerumExtract-100%EtOH-20210908/0_K9/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d9-Choline+HomoPhe+SerumExtract-100%EtOH-20210908/0_K10/1/Spectrum.mzXML',
    '2021-10-24/25DHB+d9-Choline+HomoPhe+SerumExtract-100%EtOH-20210908/0_K11/1/Spectrum.mzXML'
)
mzOrigFileName <- addParentDirectory(mzOrigFileName, mzDataDir)

mzDestFileName <- changeFilenameSuffix(mzOrigFileName, 
                                       '_processed.rds', '.mzXML')

mzPeakFileName <- changeFilenameSuffix(mzDestFileName, '-peaks.csv', '.rds')

mzRanges <- list(
    c(0, 1000)  # General
)

mzPeakQuantificationSNR <- 3  # General quantification
mzPeakMaxWidth <- 0.15

mzPeakTolerance <- 0.1  # General
mzPeakRelativeTolerance <- FALSE  # General


### MAIN ENTRY ###
# Parse the original data file, applying smoothing algorithm, etc.
processMassSpectra(mzOrigFileName, mzDestFileName,
                   mzRanges = mzRanges,
                   baseline = FALSE, smoothing = TRUE)

# Peak picking with specified criteria
pickPeaks(mzDestFileName, mzPeakFileName, 
          mzRanges, mzPeakQuantificationSNR, mzPeakMaxWidth)