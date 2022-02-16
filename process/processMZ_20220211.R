# R script for pre-processing 1-D MS spectra

source('lib/PathUtils.R')
source('lib/MassSpectraProcessor.R')


mzDataDir <- '../data/spectra/MALDI-TOF'

mzOrigFilenameList <- list(
    c('2022-02-11/TiO2NP25+EtOH+NaCl/0_O4/1/Spectrum.mzXML',
      '2022-02-11/TiO2NP25+EtOH+NaCl/0_O5/1/Spectrum.mzXML',
      '2022-02-11/TiO2NP25+EtOH+NaCl/0_O6/1/Spectrum.mzXML',
      '2022-02-11/TiO2NP25+EtOH+NaCl/0_P4/1/Spectrum.mzXML',
      '2022-02-11/TiO2NP25+EtOH+NaCl/0_P5/1/Spectrum.mzXML'),
    c('2022-02-11/TiO2NP25+NaCl+PC(18(1)-18(1))/0_N1/1/Spectrum.mzXML',
      '2022-02-11/TiO2NP25+NaCl+PC(18(1)-18(1))/0_N2/1/Spectrum.mzXML',
      '2022-02-11/TiO2NP25+NaCl+PC(18(1)-18(1))/0_N3/1/Spectrum.mzXML',
      '2022-02-11/TiO2NP25+NaCl+PC(18(1)-18(1))/0_N4/1/Spectrum.mzXML',
      '2022-02-11/TiO2NP25+NaCl+PC(18(1)-18(1))/0_N5/1/Spectrum.mzXML')
)
mzOrigFilenameList <- addParentDirectory(mzOrigFilenameList, mzDataDir)

mzDestFilenameList <- changeFilenameSuffix(mzOrigFilenameList, 
                                           '_processed.rds', '.mzXML')

mzPeakFilenameList <- changeFilenameSuffix(mzDestFilenameList, '-peaks.csv', '.rds')
mzPeakFilenameList2 <- changeFilenameSuffix(mzPeakFilenameList, '_noBG.csv', '.csv')

mzBackgroundFilenameList <- list(
    c(),
    c('2022-02-11/TiO2NP25+EtOH+NaCl/0_O4/1/Spectrum.mzXML',
      '2022-02-11/TiO2NP25+EtOH+NaCl/0_O5/1/Spectrum.mzXML',
      '2022-02-11/TiO2NP25+EtOH+NaCl/0_O6/1/Spectrum.mzXML',
      '2022-02-11/TiO2NP25+EtOH+NaCl/0_P4/1/Spectrum.mzXML',
      '2022-02-11/TiO2NP25+EtOH+NaCl/0_P5/1/Spectrum.mzXML')
)
mzBackgroundPeakFilenameList <- 
            changeFilenameSuffix(mzBackgroundFilenameList, 
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
mzPeakMinimumReplicate <- 4/5  # General


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
