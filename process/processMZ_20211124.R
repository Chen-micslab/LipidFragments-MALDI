# R script for pre-processing 1-D MS spectra

source('lib/PathUtils.R')
source('lib/MassSpectraProcessor.R')


mzDataDir <- '../data/spectra/MALDI-TOF'

mzOrigFilenameList <- list(
    c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/135.1000.LIFT/Spectrum.mzXML',
      '2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/173.0000.LIFT/Spectrum.mzXML',
      '2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/177.0000.LIFT/Spectrum.mzXML',
      '2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/198.1000.LIFT/Spectrum.mzXML',
      '2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/217.0000.LIFT/Spectrum.mzXML',
      '2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/339.2000.LIFT/Spectrum.mzXML',
      '2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/463.2000.LIFT/Spectrum.mzXML'),
    c('2021-11-24/TiO2NP25+NaCl+PC(18(0)-18(2))-MS2/0_E1/1/198.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-24/TiO2NP25+NaCl+PE(18-18)-MS2/0_B4/1/173.0000.LIFT/Spectrum.mzXML',
      '2021-11-24/TiO2NP25+NaCl+PG(18(0)-18(1))-MS2/0_D4/1/173.0000.LIFT/Spectrum.mzXML'),
    c('2021-11-24/TiO2NP25+NaCl+PE(18-18)-MS2/0_B5/1/339.2000.LIFT/Spectrum.mzXML'),
    c('2021-11-24/TiO2NP25+NaCl+PE(18-18)-MS2/0_D1/1/463.2000.LIFT/Spectrum.mzXML'),
    c('2021-11-24/TiO2NP25+NaCl+PG(18(0)-18(1))-MS2/0_D4/1/135.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-24/TiO2NP25+NaCl+PG(18(0)-18(1))-MS2/0_D4/1/177.0000.LIFT/Spectrum.mzXML'),
    c('2021-11-24/TiO2NP25+NaCl+PG(18(0)-18(1))-MS2/0_D4/1/217.0000.LIFT/Spectrum.mzXML')
)
mzOrigFilenameList <- addParentDirectory(mzOrigFilenameList, mzDataDir)

mzDestFilenameList <- changeFilenameSuffix(mzOrigFilenameList, 
                                           '_processed.rds', '.mzXML')

mzPeakFilenameList <- changeFilenameSuffix(mzDestFilenameList, '-peaks.csv', '.rds')
mzPeakFilenameList2 <- changeFilenameSuffix(mzPeakFilenameList, '_noBG.csv', '.csv')

mzBackgroundPeakFilenameList <- list(
    c(),
    c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/198.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/173.0000.LIFT/Spectrum.mzXML'),
    c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/339.2000.LIFT/Spectrum.mzXML'),
    c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/463.2000.LIFT/Spectrum.mzXML'),
    c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/135.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/177.0000.LIFT/Spectrum.mzXML'),
    c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/217.0000.LIFT/Spectrum.mzXML')
)
mzBackgroundPeakFilenameList <- 
            changeFilenameSuffix(mzBackgroundPeakFilenameList, 
                                 '_processed-peaks.csv', '.mzXML')
mzBackgroundPeakFilenameList <- 
            addParentDirectory(mzBackgroundPeakFilenameList, mzDataDir)


mzRanges <- list(
    c(20, 1000)  # General
)

mzPeakSNR <- 3  # General
mzPeakMaxWidth <- 0.2  # General

mzPeakTolerance <- 0.1  # General
mzPeakRelativeTolerance <- FALSE  # General
mzPeakMinimumReplicate <- 1


### MAIN ENTRY ###
# Parse the original data file, applying smoothing algorithm, etc.
processMassSpectra(unlist(mzOrigFilenameList), unlist(mzDestFilenameList),
                   mzRanges = mzRanges,
                   baseline = FALSE, smoothing = TRUE)

# Peak picking with specified criteria
pickPeaks(unlist(mzDestFilenameList), unlist(mzPeakFilenameList), 
          mzRanges, mzPeakSNR, mzPeakMaxWidth)

# Background substraction
for (i in seq(1, length(mzPeakFilenameList)))
{
    if (i > length(mzBackgroundPeakFilenameList))
        break()
    if (length(mzBackgroundPeakFilenameList[[i]]) < 1)
        next()
    
    # Calculate an averaged background peak list
    bgPeakList <- list()
    for (bgFilename in mzBackgroundPeakFilenameList[[i]])
    {
        bgPeaks <- read.csv(bgFilename, as.is = TRUE)
        bgPeakList <- c(bgPeakList, list(bgPeaks[,1]))
    }
    bgPeakList <- averageMZ(bgPeakList, 
                            tolerance = mzPeakTolerance, 
                            relativeTolerance = mzPeakRelativeTolerance,
                            minFrequency = mzPeakMinimumReplicate,
                            algorithm = median)
    
    # Substract background peaks from each peak file
    peakList <- list()
    for (j in seq(1, length(mzPeakFilenameList[[i]])))
    {
        peaks <- read.csv(mzPeakFilenameList[[i]][j], as.is = TRUE)
        mzMasks <- !selectMZ(peaks[,1], bgPeakList, 
                             tolerance = mzPeakTolerance,
                             relativeTolerance = mzPeakRelativeTolerance)
        write.csv(peaks[mzMasks,], mzPeakFilenameList2[[i]][j], 
                  row.names = FALSE)
    }
}