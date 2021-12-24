# R script for pre-processing 1-D MS spectra

source('lib/PathUtils.R')
source('lib/MassSpectraProcessor.R')


mzDataDir <- '../data/spectra/MALDI-TOF'

mzOrigFileName <- c(
    '2021-10-25/25DHB+d9-Choline+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I13/1/Spectrum.mzXML',
    '2021-10-25/25DHB+d9-Choline+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I14/1/Spectrum.mzXML',
    '2021-10-25/25DHB+d9-Choline+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I15/1/Spectrum.mzXML',
    '2021-10-25/25DHB+d9-Choline+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I16/1/Spectrum.mzXML',
    '2021-10-25/25DHB+d9-Choline+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I17/1/Spectrum.mzXML',
    '2021-10-25/25DHB+d3-Carnitine+SerumExtract-100%EtOH-20210908/0_J13/1/Spectrum.mzXML',
    '2021-10-25/25DHB+d3-Carnitine+SerumExtract-100%EtOH-20210908/0_J14/1/Spectrum.mzXML',
    '2021-10-25/25DHB+d3-Carnitine+SerumExtract-100%EtOH-20210908/0_J15/1/Spectrum.mzXML',
    '2021-10-25/25DHB+d3-Carnitine+SerumExtract-100%EtOH-20210908/0_J16/1/Spectrum.mzXML',
    '2021-10-25/25DHB+d3-Carnitine+SerumExtract-100%EtOH-20210908/0_J17/1/Spectrum.mzXML',
    '2021-10-25/25DHB+d3-Carnitine-SerumExtract-50%CHCl3-50%MeOH-20210908/0_L13/1/Spectrum.mzXML',
    '2021-10-25/25DHB+d3-Carnitine-SerumExtract-50%CHCl3-50%MeOH-20210908/0_L14/1/Spectrum.mzXML',
    '2021-10-25/25DHB+d3-Carnitine-SerumExtract-50%CHCl3-50%MeOH-20210908/0_L15/1/Spectrum.mzXML',
    '2021-10-25/25DHB+d3-Carnitine-SerumExtract-50%CHCl3-50%MeOH-20210908/0_L16/1/Spectrum.mzXML',
    '2021-10-25/25DHB+d3-Carnitine-SerumExtract-50%CHCl3-50%MeOH-20210908/0_L17/1/Spectrum.mzXML'
)
mzOrigFileName <- addParentDirectory(mzOrigFileName, mzDataDir)

mzDestFileName <- changeFilenameSuffix(mzOrigFileName, 
                                       '_processed.rds', '.mzXML')

mzPeakFileName <- changeFilenameSuffix(mzDestFileName, '-peaks.csv', '.rds')

mzRanges <- list(
    c(0, 1000)  # General
)

mzPeakQuantificationSNR <- 3  # General
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