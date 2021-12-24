# R script for comparing untargeted peaks in MALDI spectra 
# obtained by different final concentration of analytes 

source('lib/DataFilter.R')
source('lib/MSUtils.R')
source('lib/PathUtils.R')

mzDataDir <- '../data/spectra/MALDI-TOF'
mzStatDir <- '../data/summary/MALDI-TOF'

mzPeakFileNameList <- list(
    c('2021-11-26/25DHB+SPCS-50%CHCl3-50%MeOH-20210908/0_H1/1/Spectrum_processed-peaks.csv',
      '2021-11-26/25DHB+SPCS-50%CHCl3-50%MeOH-20210908/0_H2/1/Spectrum_processed-peaks.csv',
      '2021-11-26/25DHB+SPCS-50%CHCl3-50%MeOH-20210908/0_H3/1/Spectrum_processed-peaks.csv',
      '2021-11-26/25DHB+SPCS-50%CHCl3-50%MeOH-20210908/0_H4/1/Spectrum_processed-peaks.csv',
      '2021-11-26/25DHB+SPCS-50%CHCl3-50%MeOH-20210908/0_H5/1/Spectrum_processed-peaks.csv'),
    c('2021-11-26/25DHB+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I1/2/Spectrum_processed-peaks_noBG.csv',
      '2021-11-26/25DHB+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I2/1/Spectrum_processed-peaks_noBG.csv',
      '2021-11-26/25DHB+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I3/1/Spectrum_processed-peaks_noBG.csv',
      '2021-11-26/25DHB+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I4/1/Spectrum_processed-peaks_noBG.csv',
      '2021-11-26/25DHB+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I5/1/Spectrum_processed-peaks_noBG.csv')
)
mzPeakFileNameList <- addParentDirectory(mzPeakFileNameList, mzDataDir)

mzDestPeakFileName <- c(
    '25DHB+SPCS-50%CHCl3-50%MeOH-20210908.csv',
    '25DHB+SerumExtract-50%CHCl3-50%MeOH-20210908_noBG.csv'
)
mzDestPeakFileName <- addParentDirectory(mzDestPeakFileName, mzStatDir)

mzRanges <- list(
    c(0, 1000)  # General
)
mzPeakTolerance <- 0.05
mzPeakRelativeTolerance <- FALSE
mzPeakIntraGroupMinimumReplicate <- 4/5
mzNormalizationReference <- c(
)


selectDataByIndexes <- function(index, data, indexList, threshold, 
                                algorithm, ...)
{
    cluster <- unlist(sapply(seq(1, length(indexList)),
                             probe = index,
                             probeList = indexList,
                             dataList = data,
                             function(i, probe, probeList, dataList)
                             {
                                 return(dataList[[i]][probeList[[i]] == probe])
                             }))
    if (length(cluster) / length(data) >= threshold)
        return(algorithm(cluster, ...))
    else
        return(NULL)
}
normalizeSpectra <- function(mzList, intensityList,
                             referenceMZ, referenceTolerance, relativeTolerance)
{
    spectraCount <- min(length(mzList), length(intensityList))
    if (spectraCount <= 0)
        return(NULL)
    return(
        lapply(seq(1, spectraCount),
               function(i)
               {
                   mzMasks <- selectMZ(mzList[[i]], referenceMZ,
                                       tolerance = referenceTolerance,
                                       relativeTolerance = relativeTolerance)
                   if (all(!mzMasks))
                       return(rep(NA, length(intensityList[[i]])))
                   refIntensity <- intensityList[[i]][mzMasks]
                   return(intensityList[[i]] / mean(refIntensity))
               }))
}


### MAIN ENTRY ###
for (i in seq(1, length(mzPeakFileNameList)))
{
    # Read m/z and intensities from data files
    mzList <- list()
    intensityList <- list()
    
    for (mzFile in mzPeakFileNameList[[i]])
    {
        mzPeakFile <- read.csv(mzFile, as.is = TRUE)
        mzIndexes <- rangeFilter.index(mzPeakFile[,1], mzRanges)
        mzList <- c(mzList, list(mzPeakFile[mzIndexes, 1]))
        intensityList <- c(intensityList, list(mzPeakFile[mzIndexes, 4]))
        # intensityList <- c(intensityList, list(mzPeakFile[mzIndexes, 5]))
    }
    
    # Normalization of intensities with respect to various references
    normalizedIntensityList <- lapply(mzNormalizationReference,
                                      referenceTolerance = mzPeakTolerance, 
                                      relativeTolerance = mzPeakRelativeTolerance,
                                      mzList = mzList, 
                                      intensityList = intensityList,
                                      normalizeSpectra)
    normalizedIntensityList <- 
        c(list(none = intensityList,
               TIC = lapply(intensityList, function(X) {return(X / sum(X))})), 
          normalizedIntensityList)
    
    # Clustering of m/z
    mzClusterIndexes <- clusterMZ(mzList, 
                                  mzPeakTolerance, 
                                  mzPeakRelativeTolerance)
    
    # Averaging the spectra
    allIndexes <- unique(unlist(mzClusterIndexes))
    averagedMZ <- unlist(lapply(allIndexes,
                         data = mzList,
                         indexList = mzClusterIndexes,
                         threshold = mzPeakIntraGroupMinimumReplicate,
                         algorithm = mean,
                         selectDataByIndexes))
    averagedIntensityList <- 
        lapply(normalizedIntensityList,
               function(X)
               {
                   unlist(sapply(allIndexes,
                                 data = X,
                                 indexList = mzClusterIndexes,
                                 threshold = mzPeakIntraGroupMinimumReplicate,
                                 algorithm = mean,
                                 na.rm = TRUE,
                                 selectDataByIndexes))
               })
    stdevIntensities <- 
        lapply(normalizedIntensityList,
               function(X)
               {
                   unlist(sapply(allIndexes,
                                 data = X,
                                 indexList = mzClusterIndexes,
                                 threshold = mzPeakIntraGroupMinimumReplicate,
                                 algorithm = sd,
                                 na.rm = TRUE,
                                 selectDataByIndexes))
               })
    
    # Write found m/z
    if (length(averagedMZ) > 0)
    {
        mzOrder <- order(averagedMZ)
        outputData <- data.frame(averagedMZ[mzOrder])
        for (j in seq(1, length(averagedIntensityList)))
        {
            outputData <- cbind(outputData,
                                averagedIntensityList[[j]][mzOrder],
                                stdevIntensities[[j]][mzOrder] / 
                                    averagedIntensityList[[j]][mzOrder])
        }
        colnames(outputData) <- c('m/z', 
                                  rep(c('Intensity', 'CV'), 
                                      length(averagedIntensityList)))
        write.csv(outputData, mzDestPeakFileName[i], row.names = FALSE)
    } else
        warning('Peak file ', mzDestPeakFileName[i], ' is empty; skip it.')
}