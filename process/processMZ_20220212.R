# R script for pre-processing 1-D MS spectra

source('lib/PathUtils.R')
source('lib/MassSpectraProcessor.R')


mzDataDir <- '../data/spectra/MALDI-TOF'

mzOrigFilenameList <- list(
    c('2022-02-11/TiO2NP25+EtOH+NaCl-MS2/0_O1/1/104.1000.LIFT/Spectrum.mzXML',
      '2022-02-11/TiO2NP25+EtOH+NaCl-MS2/0_O2/1/125.0000.LIFT/Spectrum.mzXML',
      '2022-02-12/TiO2NP25+EtOH+NaCl-MS2/0_A7/1/147.0000.LIFT/Spectrum.mzXML',
      '2022-02-12/TiO2NP25+EtOH+NaCl-MS2/0_A8/1/166.1000.LIFT/Spectrum.mzXML',
      '2022-02-10/TiO2NP25+EtOH+NaCl-MS2/0_J4/1/184.1000.LIFT/Spectrum.mzXML',
      '2022-02-10/TiO2NP25+EtOH+NaCl-MS2/0_J5/1/198.1000.LIFT/Spectrum.mzXML'),
    c('2022-02-11/TiO2NP25+NaCl+PC(18(1)-18(1))-MS2/0_M2/1/104.1000.LIFT/Spectrum.mzXML'),
    c('2022-02-11/TiO2NP25+NaCl+PC(18(1)-18(1))-MS2/0_M1/1/125.0000.LIFT/Spectrum.mzXML'),
    c('2022-02-12/TiO2NP25+NaCl+PC(18(1)-18(1))-MS2/0_B7/1/147.0000.LIFT/Spectrum.mzXML'),
    c('2022-02-12/TiO2NP25+NaCl+PC(18(1)-18(1))-MS2/0_B11/1/166.1000.LIFT/Spectrum.mzXML'),
    c('2022-02-10/TiO2NP25+NaCl+PC(18(1)-18(1))-MS2/0_K4/1/184.1000.LIFT/Spectrum.mzXML'),
    c('2022-02-10/TiO2NP25+NaCl+PC(18(1)-18(1))-MS2/0_K5/1/198.1000.LIFT/Spectrum.mzXML')
)
mzOrigFilenameList <- addParentDirectory(mzOrigFilenameList, mzDataDir)

mzDestFilenameList <- changeFilenameSuffix(mzOrigFilenameList, 
                                           '_processed.rds', '.mzXML')

mzPeakFilenameList <- changeFilenameSuffix(mzDestFilenameList, '-peaks.csv', '.rds')
mzPeakFilenameList2 <- changeFilenameSuffix(mzPeakFilenameList, '_noBG.csv', '.csv')

mzBackgroundFilenameList <- list(
    c(),
    c('2022-02-11/TiO2NP25+EtOH+NaCl-MS2/0_O1/1/104.1000.LIFT/Spectrum.mzXML'),
    c('2022-02-11/TiO2NP25+EtOH+NaCl-MS2/0_O2/1/125.0000.LIFT/Spectrum.mzXML'),
    c('2022-02-12/TiO2NP25+EtOH+NaCl-MS2/0_A7/1/147.0000.LIFT/Spectrum.mzXML'),
    c('2022-02-12/TiO2NP25+EtOH+NaCl-MS2/0_A8/1/166.1000.LIFT/Spectrum.mzXML'),
    c('2022-02-10/TiO2NP25+EtOH+NaCl-MS2/0_J4/1/184.1000.LIFT/Spectrum.mzXML'),
    c('2022-02-10/TiO2NP25+EtOH+NaCl-MS2/0_J5/1/198.1000.LIFT/Spectrum.mzXML')
)
mzBackgroundPeakFilenameList <- 
            changeFilenameSuffix(mzBackgroundFilenameList, 
                                 '_processed-peaks.csv', '.mzXML')
mzBackgroundPeakFilenameList <- 
            addParentDirectory(mzBackgroundPeakFilenameList, mzDataDir)


mzRanges <- list(
    c(20, 1000)  # General
)

mzPeakSNR <- 3  # General
mzPeakMaxWidth <- 0.2  # MS/MS

mzPeakTolerance <- 0.1  # MS/MS
mzPeakRelativeTolerance <- FALSE  # General
mzPeakMinimumReplicate <- 1/1  # MS/MS


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
