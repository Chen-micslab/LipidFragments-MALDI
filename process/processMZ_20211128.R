# R script for pre-processing 1-D MS spectra

source('lib/PathUtils.R')
source('lib/MassSpectraProcessor.R')


mzDataDir <- '../data/spectra/MALDI-TOF'

mzOrigFilenameList <- list(
    c('2021-11-28/25DHB@MeOH-MS2/0_B8/1/104.1000.LIFT/Spectrum.mzXML',
      '2021-11-28/25DHB@MeOH-MS2/0_B10/1/184.1000.LIFT/Spectrum.mzXML',
      '2021-11-28/25DHB+SPCS-100%EtOH-20210908-MS2/0_C7/1/104.1000.LIFT/Spectrum.mzXML',
      '2021-11-28/25DHB+SPCS-100%EtOH-20210908-MS2/0_C15/1/166.1000.LIFT/Spectrum.mzXML',
      '2021-11-28/25DHB+SPCS-100%EtOH-20210908-MS2/0_C19/1/184.1000.LIFT/Spectrum.mzXML',
      '2021-11-28/25DHB+SPCS-50%CHCl3-50%MeOH-20210908-MS2/0_E7/1/104.1000.LIFT/Spectrum.mzXML',
      '2021-11-28/25DHB+SPCS-50%CHCl3-50%MeOH-20210908-MS2/0_E18/1/166.1000.LIFT/Spectrum.mzXML',
      '2021-11-28/25DHB+SPCS-50%CHCl3-50%MeOH-20210908-MS2/0_E20/1/184.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB+Carnitine-MS2/0_B14/1/184.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB+CholineChloride-MS2/0_O6/1/104.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB+SerumExtract-100%EtOH-20210908-MS2/0_D7/1/104.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB+SerumExtract-100%EtOH-20210908-MS2/0_D15/1/166.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB+SerumExtract-100%EtOH-20210908-MS2/0_D19/1/184.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB+SerumExtract-50%CHCl3-50%MeOH-20210908-MS2/0_F7/1/104.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB+SerumExtract-50%CHCl3-50%MeOH-20210908-MS2/0_F18/1/166.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB+SerumExtract-50%CHCl3-50%MeOH-20210908-MS2/0_F21/1/184.1000.LIFT/Spectrum.mzXML')
)
mzOrigFilenameList <- addParentDirectory(mzOrigFilenameList, mzDataDir)

mzDestFilenameList <- changeFilenameSuffix(mzOrigFilenameList, 
                                           '_processed.rds', '.mzXML')

mzPeakFilenameList <- changeFilenameSuffix(mzDestFilenameList, '-peaks.csv', '.rds')
mzPeakFilenameList2 <- changeFilenameSuffix(mzPeakFilenameList, '_noBG.csv', '.csv')

mzBackgroundPeakFilenameList <- list(
    c(),
    c('2021-11-28/25DHB@MeOH-MS2/0_B10/1/184.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB@MeOH-MS2/0_B8/1/104.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB+SPCS-100%EtOH-20210908-MS2/0_C7/1/104.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB+SPCS-100%EtOH-20210908-MS2/0_C15/1/166.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB+SPCS-100%EtOH-20210908-MS2/0_C19/1/184.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB+SPCS-50%CHCl3-50%MeOH-20210908-MS2/0_E7/1/104.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB+SPCS-50%CHCl3-50%MeOH-20210908-MS2/0_E18/1/166.1000.LIFT/Spectrum.mzXML'),
    c('2021-11-28/25DHB+SPCS-50%CHCl3-50%MeOH-20210908-MS2/0_E20/1/184.1000.LIFT/Spectrum.mzXML')
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
