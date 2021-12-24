# R script for pre-processing 1-D MS spectra

source('lib/PathUtils.R')
source('lib/MassSpectraProcessor.R')


mzDataDir <- '../data/spectra/MALDI-TOF'

mzOrigFilenameList <- list(
    c('2021-11-07/TiO2NP25+EtOH+NaCl/0_C13/1/Spectrum.mzXML',
      '2021-11-07/TiO2NP25+EtOH+NaCl/0_C14/1/Spectrum.mzXML',
      '2021-11-07/TiO2NP25+EtOH+NaCl/0_C15/1/Spectrum.mzXML',
      '2021-11-07/TiO2NP25+EtOH+NaCl/0_C16/1/Spectrum.mzXML',
      '2021-11-07/TiO2NP25+EtOH+NaCl/0_C17/1/Spectrum.mzXML'),
    c('2021-11-07/TiO2NP25+NaCl+PE(18-18)/0_E13/1/Spectrum.mzXML',
      '2021-11-07/TiO2NP25+NaCl+PE(18-18)/0_E14/1/Spectrum.mzXML',
      '2021-11-07/TiO2NP25+NaCl+PE(18-18)/0_E15/1/Spectrum.mzXML',
      '2021-11-07/TiO2NP25+NaCl+PE(18-18)/0_E16/1/Spectrum.mzXML',
      '2021-11-07/TiO2NP25+NaCl+PE(18-18)/0_E17/1/Spectrum.mzXML')
)
mzOrigFilenameList <- addParentDirectory(mzOrigFilenameList, mzDataDir)

mzDestFilenameList <- changeFilenameSuffix(mzOrigFilenameList, 
                                           '_processed.rds', '.mzXML')

mzPeakFilenameList <- changeFilenameSuffix(mzDestFilenameList, '-peaks.csv', '.rds')
mzPeakFilenameList2 <- changeFilenameSuffix(mzPeakFilenameList, '_noBG.csv', '.csv')

mzBackgroundPeakFilenameList <- c(
    list(c()),
    list(c()),
    list(c()),
    rep(list(c('2021-11-07/TiO2NP25+EtOH+NaCl/0_C13/1/Spectrum_processed-peaks.csv',
               '2021-11-07/TiO2NP25+EtOH+NaCl/0_C14/1/Spectrum_processed-peaks.csv',
               '2021-11-07/TiO2NP25+EtOH+NaCl/0_C15/1/Spectrum_processed-peaks.csv',
               '2021-11-07/TiO2NP25+EtOH+NaCl/0_C16/1/Spectrum_processed-peaks.csv',
               '2021-11-07/TiO2NP25+EtOH+NaCl/0_C17/1/Spectrum_processed-peaks.csv')),
        3)
)
mzBackgroundPeakFilenameList <- addParentDirectory(mzBackgroundPeakFilenameList, 
                                                   mzDataDir)


mzRanges <- list(
    c(0, 1000)  # General
)

mzPeakQuantificationSNR <- 5  # General
mzPeakMaxWidth <- 0.15  # General

mzPeakTolerance <- 0.1  # General
mzPeakRelativeTolerance <- FALSE  # General
mzPeakMinimumReplicate <- 0.8


### MAIN ENTRY ###
# Parse the original data file, applying smoothing algorithm, etc.
processMassSpectra(unlist(mzOrigFilenameList), unlist(mzDestFilenameList),
                   mzRanges = mzRanges,
                   baseline = FALSE, smoothing = TRUE)

# Peak picking with specified criteria
pickPeaks(unlist(mzDestFilenameList), unlist(mzPeakFilenameList), 
          mzRanges, mzPeakQuantificationSNR, mzPeakMaxWidth)

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