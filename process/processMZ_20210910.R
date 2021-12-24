# R script for pre-processing 1-D MS spectra

source('lib/PathUtils.R')
source('lib/MassSpectraProcessor.R')


mzDataDir <- '../data/spectra/MALDI-TOF'

mzOrigFileName <- c(
    '2021-09-10/25DHB+d9-Choline+HomoPhe+SPCS/0_L13/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+SPCS/0_L14/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+SPCS/0_L15/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+SPCS/0_L16/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+SPCS/0_L17/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(10uM)/0_K1/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(10uM)/0_K2/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(10uM)/0_K3/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(10uM)/0_K4/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(10uM)/0_K5/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(20uM)/0_K7/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(20uM)/0_K8/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(20uM)/0_K9/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(20uM)/0_K10/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(20uM)/0_K11/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(30uM)/0_K13/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(30uM)/0_K14/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(30uM)/0_K15/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(30uM)/0_K16/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(30uM)/0_K17/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(40uM)/0_L1/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(40uM)/0_L2/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(40uM)/0_L3/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(40uM)/0_L4/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(40uM)/0_L5/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(50uM)/0_L7/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(50uM)/0_L8/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(50uM)/0_L9/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(50uM)/0_L10/1/Spectrum.mzXML',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(50uM)/0_L11/1/Spectrum.mzXML'
)
mzOrigFileName <- addParentDirectory(mzOrigFileName, mzDataDir)

mzDestFileName <- changeFilenameSuffix(mzOrigFileName, 
                                       '_processed.rds', '.mzXML')

mzPeakFileName <- changeFilenameSuffix(mzDestFileName, '-peaks.csv', '.rds')

mzRanges <- list(
    c(0, 1000)  # General
)

peakQuantificationSNR <- 6  # General quantification
peakMaxWidth <- 0.2


### MAIN ENTRY ###
# Parse the original data file, applying smoothing algorithm, etc.
processMassSpectra(mzOrigFileName, mzDestFileName,
                   mzRanges = mzRanges,
                   baseline = FALSE, smoothing = TRUE)

# Peak picking with specified criteria
pickPeaks(mzDestFileName, mzPeakFileName, 
          mzRanges, peakQuantificationSNR, peakMaxWidth)