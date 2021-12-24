# R script for comparing untargeted peaks in MALDI spectra 
# obtained by different final concentration of analytes 

source('lib/DataFilter.R')
source('lib/MSUtils.R')
source('lib/PathUtils.R')


mzDataDir <- '../data/spectra/MALDI-TOF'
mzStatDir <- '../data/summary/MALDI-TOF'

mzPeakFileNameList <- list(
    c('2021-09-10/25DHB+d9-Choline+HomoPhe+SPCS/0_L13/1/Spectrum_processed-peaks.csv',
      '2021-09-10/25DHB+d9-Choline+HomoPhe+SPCS/0_L14/1/Spectrum_processed-peaks.csv',
      '2021-09-10/25DHB+d9-Choline+HomoPhe+SPCS/0_L15/1/Spectrum_processed-peaks.csv',
      '2021-09-10/25DHB+d9-Choline+HomoPhe+SPCS/0_L16/1/Spectrum_processed-peaks.csv',
      '2021-09-10/25DHB+d9-Choline+HomoPhe+SPCS/0_L17/1/Spectrum_processed-peaks.csv'),
    c('2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(10uM)/0_K1/1/Spectrum_processed-peaks.csv',
      '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(10uM)/0_K2/1/Spectrum_processed-peaks.csv',
      '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(10uM)/0_K3/1/Spectrum_processed-peaks.csv',
      '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(10uM)/0_K4/1/Spectrum_processed-peaks.csv',
      '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(10uM)/0_K5/1/Spectrum_processed-peaks.csv'),
    c('2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(20uM)/0_K7/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(20uM)/0_K8/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(20uM)/0_K9/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(20uM)/0_K10/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(20uM)/0_K11/1/Spectrum_processed-peaks.csv'),
    c('2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(30uM)/0_K13/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(30uM)/0_K14/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(30uM)/0_K15/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(30uM)/0_K16/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(30uM)/0_K17/1/Spectrum_processed-peaks.csv'),
    c('2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(40uM)/0_L1/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(40uM)/0_L2/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(40uM)/0_L3/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(40uM)/0_L4/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(40uM)/0_L5/1/Spectrum_processed-peaks.csv'),
    c('2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(50uM)/0_L7/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(50uM)/0_L8/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(50uM)/0_L9/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(50uM)/0_L10/1/Spectrum_processed-peaks.csv',
    '2021-09-10/25DHB+d9-Choline+HomoPhe+Choline(50uM)/0_L11/1/Spectrum_processed-peaks.csv'),
    c('2021-10-24/25DHB+d9-Choline+HomoPhe+SerumExtract-100%EtOH-20210908/0_K7/1/Spectrum_processed-peaks.csv',
      '2021-10-24/25DHB+d9-Choline+HomoPhe+SerumExtract-100%EtOH-20210908/0_K8/1/Spectrum_processed-peaks.csv',
      '2021-10-24/25DHB+d9-Choline+HomoPhe+SerumExtract-100%EtOH-20210908/0_K9/1/Spectrum_processed-peaks.csv',
      '2021-10-24/25DHB+d9-Choline+HomoPhe+SerumExtract-100%EtOH-20210908/0_K10/1/Spectrum_processed-peaks.csv',
      '2021-10-24/25DHB+d9-Choline+HomoPhe+SerumExtract-100%EtOH-20210908/0_K11/1/Spectrum_processed-peaks.csv'),
    c('2021-10-25/25DHB+d9-Choline+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I13/1/Spectrum_processed-peaks.csv',
      '2021-10-25/25DHB+d9-Choline+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I14/1/Spectrum_processed-peaks.csv',
      '2021-10-25/25DHB+d9-Choline+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I15/1/Spectrum_processed-peaks.csv',
      '2021-10-25/25DHB+d9-Choline+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I16/1/Spectrum_processed-peaks.csv',
      '2021-10-25/25DHB+d9-Choline+SerumExtract-50%CHCl3-50%MeOH-20210908/0_I17/1/Spectrum_processed-peaks.csv')
)
mzPeakFileNameList <- addParentDirectory(mzPeakFileNameList, mzDataDir)

mzDestPeakFileName <- c(
    '25DHB+SPCS+d9-Choline+HomoPhe.csv',
    '25DHB+SPCS+d9-Choline+HomoPhe+Choline(10uM).csv',
    '25DHB+SPCS+d9-Choline+HomoPhe+Choline(20uM).csv',
    '25DHB+SPCS+d9-Choline+HomoPhe+Choline(30uM).csv',
    '25DHB+SPCS+d9-Choline+HomoPhe+Choline(40uM).csv',
    '25DHB+SPCS+d9-Choline+HomoPhe+Choline(50uM).csv',
    '25DHB+d9-Choline+SerumExtract-100%EtOH-20210908.csv',
    '25DHB+d9-Choline+SerumExtract-50%CHCl3-50%MeOH-20210908.csv'
)
mzDestPeakFileName <- addParentDirectory(mzDestPeakFileName, mzStatDir)

mzRanges <- list(
    c(0, 1000)
)
mzPeakTolerance <- 0.05
mzPeakRelativeTolerance <- FALSE
mzPeakIntraGroupMinimumReplicate <- 4/5
mzNormalizationReference <- c(
    113.1640  # d9-Choline, [C5H5{2H}9NO]+
)


selectDataByIndexes <- function(index, data, indexList, threshold, algorithm)
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
        return(algorithm(cluster))
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