# R script for pre-processing 1-D MS spectra

source('lib/PathUtils.R')
source('lib/MassSpectraProcessor.R')


mzDataDir <- '../data/spectra/MALDI-TOF'

mzOrigFilenameList <- list(
    c('2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J16/1/104.1000.LIFT/Spectrum.mzXML',
      '2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J16/1/125.0000.LIFT/Spectrum.mzXML',
      '2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J16/1/166.1000.LIFT/Spectrum.mzXML',
      '2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J16/1/184.1000.LIFT/Spectrum.mzXML',
      '2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J17/1/104.1000.LIFT/Spectrum.mzXML',
      '2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J17/1/147.0000.LIFT/Spectrum.mzXML',
      '2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J17/1/184.1000.LIFT/Spectrum.mzXML',
      '2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J17/1/198.1000.LIFT/Spectrum.mzXML'),
    c('2021-12-25/TiO2NP25+NaCl+PC(18(0)-18(2))-MS2/0_K16/1/104.1000.LIFT/Spectrum.mzXML',
      '2021-12-25/TiO2NP25+NaCl+PC(18(0)-18(2))-MS2/0_K17/1/104.1000.LIFT/Spectrum.mzXML'),
    c('2021-12-25/TiO2NP25+NaCl+PC(18(0)-18(2))-MS2/0_K16/1/125.0000.LIFT/Spectrum.mzXML'),
    c('2021-12-25/TiO2NP25+NaCl+PC(18(0)-18(2))-MS2/0_K16/1/166.1000.LIFT/Spectrum.mzXML'),
    c('2021-12-25/TiO2NP25+NaCl+PC(18(0)-18(2))-MS2/0_K16/1/184.1000.LIFT/Spectrum.mzXML',
      '2021-12-25/TiO2NP25+NaCl+PC(18(0)-18(2))-MS2/0_K17/1/184.1000.LIFT/Spectrum.mzXML'),
    c('2021-12-25/TiO2NP25+NaCl+PC(18(0)-18(2))-MS2/0_K17/1/147.0000.LIFT/Spectrum.mzXML'),
    c('2021-12-25/TiO2NP25+NaCl+PC(18(0)-18(2))-MS2/0_K17/1/198.1000.LIFT/Spectrum.mzXML')
)
mzOrigFilenameList <- addParentDirectory(mzOrigFilenameList, mzDataDir)

mzDestFilenameList <- changeFilenameSuffix(mzOrigFilenameList, 
                                           '_processed.rds', '.mzXML')

mzPeakFilenameList <- changeFilenameSuffix(mzDestFilenameList, '-peaks.csv', '.rds')
mzPeakFilenameList2 <- changeFilenameSuffix(mzPeakFilenameList, '_noBG.csv', '.csv')

mzBackgroundPeakFilenameList <- list(
    c(),
    c('2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J16/1/104.1000.LIFT/Spectrum.mzXML',
      '2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J17/1/104.1000.LIFT/Spectrum.mzXML'),
    c('2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J16/1/125.0000.LIFT/Spectrum.mzXML'),
    c('2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J16/1/166.1000.LIFT/Spectrum.mzXML'),
    c('2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J16/1/184.1000.LIFT/Spectrum.mzXML',
      '2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J17/1/184.1000.LIFT/Spectrum.mzXML'),
    c('2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J17/1/147.0000.LIFT/Spectrum.mzXML'),
    c('2021-12-25/TiO2NP25+EtOH+NaCl-MS2/0_J17/1/198.1000.LIFT/Spectrum.mzXML')
)
mzBackgroundPeakFilenameList <- 
            changeFilenameSuffix(mzBackgroundPeakFilenameList, 
                                 '_processed-peaks.csv', '.mzXML')
mzBackgroundPeakFilenameList <- 
            addParentDirectory(mzBackgroundPeakFilenameList, mzDataDir)


mzRanges <- list(
    c(20, 1000)  # MS/MS
)

mzPeakSNR <- 3  # General
mzPeakMaxWidth <- 0.2  # General

mzPeakTolerance <- 0.1  # General
mzPeakRelativeTolerance <- FALSE  # General
mzPeakMinimumReplicate <- 1  # MS/MS


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
