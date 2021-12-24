# R script for comparing untargeted peaks in MALDI spectra 
# obtained on replicates of different samples

source('lib/PathUtils.R')
source('lib/DataFilter.R')
source('lib/MSUtils.R')


mzDataDir <- '../data/spectra/HPLC-MS'
mzStatDir <- '../data/summary/HPLC-MS'

mzPeakFileNameList <- list(
    '20210917/SPCS@38%MeOH-positive-3_processed-peaks.csv',
    c('20210917/SPCS+d9-Choline+HomoPhe-positive_processed-peaks.csv',
      '20210917/SPCS+d9-Choline+HomoPhe-positive-2_processed-peaks.csv',
      '20210917/SPCS+d9-Choline+HomoPhe-positive-3_processed-peaks.csv',
      '20210917/SPCS+d9-Choline+HomoPhe-positive-4_processed-peaks.csv'),
    '20210918/SPCS+d9-Choline+HomoPhe+Choline(5uM)-positive2_processed-peaks.csv',
    '20210917/SPCS+d9-Choline+HomoPhe+Choline(10uM)-positive_processed-peaks.csv',
    '20210917/SPCS+d9-Choline+HomoPhe+Choline(20uM)-positive_processed-peaks.csv',
    '20210917/SPCS+d9-Choline+HomoPhe+Choline(30uM)-positive_processed-peaks.csv',
    '20210917/SPCS+d9-Choline+HomoPhe+Choline(40uM)-positive_processed-peaks.csv',
    c('20210917/SerumExtract-100%EtOH-20210908-positive_processed-peaks.csv',
      '20210917/SerumExtract-100%EtOH-20210908-positive-2_processed-peaks.csv',
      '20210917/SerumExtract-100%EtOH-20210908-positive-3_processed-peaks.csv'),
    c('20211025/SerumExtract-50%CHCl3-50%MeOH-20210908-positive-1_processed-peaks.csv',
      '20211025/SerumExtract-50%CHCl3-50%MeOH-20210908-positive-2_processed-peaks.csv',
      '20211025/SerumExtract-50%CHCl3-50%MeOH-20210908-positive-3_processed-peaks.csv')
)
mzPeakFileNameList <- addParentDirectory(mzPeakFileNameList, mzDataDir)
    
mzDestPeakFileName <- c(
    'SPCS@38%MeOH.csv',
    'SPCS+d9-Choline+HomoPhe.csv',
    'SPCS+d9-Choline+HomoPhe+Choline(5uM).csv',
    'SPCS+d9-Choline+HomoPhe+Choline(10uM).csv',
    'SPCS+d9-Choline+HomoPhe+Choline(20uM).csv',
    'SPCS+d9-Choline+HomoPhe+Choline(30uM).csv',
    'SPCS+d9-Choline+HomoPhe+Choline(40uM).csv',
    'SerumExtract-100%EtOH-20210908+d9-Choline+HomoPhe.csv',
    'SerumExtract-50%CHCl3-50%MeOH-20210908+d9-Choline+HomoPhe.csv'
)
mzDestPeakFileName <- addParentDirectory(mzDestPeakFileName, mzStatDir)

mzRanges <- list(
    c(0, 1000)
)

mzPeakIntraGroupTolerance <- 0.2
mzPeakIntraGroupRelativeTolerance <- FALSE
mzPeakIntraGroupMinimumReplicate <- 3/3
mzPeakInterGroupTolerance <- 0.1
mzPeakInterGroupRelativeTolerance <- FALSE

mzNormalizationReference <- c(
    113.1640 # d9-Choline, [C5H5{2H}9NO]+
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
    }
    
    # Normalization of intensities with respect to various references
    normalizedIntensityList <- 
        lapply(mzNormalizationReference,
               referenceTolerance = mzPeakIntraGroupTolerance,
               relativeTolerance = mzPeakIntraGroupRelativeTolerance,
               mzList = mzList, 
               intensityList = intensityList,
               normalizeSpectra)
    normalizedIntensityList <- c(
        list(none = intensityList,
             TIC = lapply(intensityList, function(X) {return(X / sum(X))})),
        normalizedIntensityList)
    
    # Clustering of m/z
    mzClusterIndexes <- clusterMZ(mzList, 
                                  mzPeakIntraGroupTolerance,
                                  mzPeakIntraGroupRelativeTolerance)
    
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
