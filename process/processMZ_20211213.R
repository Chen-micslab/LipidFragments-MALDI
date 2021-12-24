# R script for pre-processing 1-D MS spectra

source('lib/PathUtils.R')
source('lib/MassSpectraProcessor.R')


mzDataDir <- '../data/spectra/MALDI-TOF'

mzOrigFilenameList <- list(
    c('2021-12-13/TiO2NP25+SPCS-100%EtOH-20210908/0_E20/1/Spectrum.mzXML',
      '2021-12-13/TiO2NP25+SPCS-100%EtOH-20210908/0_E21/1/Spectrum.mzXML',
      '2021-12-13/TiO2NP25+SPCS-100%EtOH-20210908/0_E22/1/Spectrum.mzXML',
      '2021-12-13/TiO2NP25+SPCS-100%EtOH-20210908/0_E23/1/Spectrum.mzXML',
      '2021-12-13/TiO2NP25+SPCS-100%EtOH-20210908/0_E24/1/Spectrum.mzXML'),
    c('2021-12-13/TiO2NP25+SerumExtract-100%EtOH-20210908/0_F20/1/Spectrum.mzXML',
      '2021-12-13/TiO2NP25+SerumExtract-100%EtOH-20210908/0_F21/1/Spectrum.mzXML',
      '2021-12-13/TiO2NP25+SerumExtract-100%EtOH-20210908/0_F22/1/Spectrum.mzXML',
      '2021-12-13/TiO2NP25+SerumExtract-100%EtOH-20210908/0_F23/1/Spectrum.mzXML',
      '2021-12-13/TiO2NP25+SerumExtract-100%EtOH-20210908/0_F24/1/Spectrum.mzXML')
)
mzOrigFilenameList <- addParentDirectory(mzOrigFilenameList, mzDataDir)

mzDestFilenameList <- changeFilenameSuffix(mzOrigFilenameList, 
                                           '_processed.rds', '.mzXML')

mzPeakFilenameList <- changeFilenameSuffix(mzDestFilenameList, '-peaks.csv', '.rds')
mzPeakFilenameList2 <- changeFilenameSuffix(mzPeakFilenameList, '_noBG.csv', '.csv')

mzBackgroundPeakFilenameList <- list(
    c(), 
    c('2021-12-13/TiO2NP25+SPCS-100%EtOH-20210908/0_E20/1/Spectrum.mzXML',
      '2021-12-13/TiO2NP25+SPCS-100%EtOH-20210908/0_E21/1/Spectrum.mzXML',
      '2021-12-13/TiO2NP25+SPCS-100%EtOH-20210908/0_E22/1/Spectrum.mzXML',
      '2021-12-13/TiO2NP25+SPCS-100%EtOH-20210908/0_E23/1/Spectrum.mzXML',
      '2021-12-13/TiO2NP25+SPCS-100%EtOH-20210908/0_E24/1/Spectrum.mzXML')
)
mzBackgroundPeakFilenameList <- 
            changeFilenameSuffix(mzBackgroundPeakFilenameList, 
                                 '_processed-peaks.csv', '.mzXML')
mzBackgroundPeakFilenameList <- 
            addParentDirectory(mzBackgroundPeakFilenameList, mzDataDir)


mzRanges <- list(
    c(0, 1000)  # General
)

mzPeakSNR <- 3  # General
mzPeakMaxWidth <- 0.15  # General

mzPeakTolerance <- 0.05  # General
mzPeakRelativeTolerance <- FALSE  # General
mzPeakMinimumReplicate <- 4/5


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
