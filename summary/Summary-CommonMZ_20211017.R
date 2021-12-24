# R script for comparing untargeted peaks in MALDI spectra 
# obtained by different final concentration of analytes 

source('lib/DataFilter.R')
source('lib/MSUtils.R')
source('lib/PathUtils.R')

mzDataDir <- '../data/spectra/HPLC-MS'
mzStatDir <- '../data/summary/HPLC-MS'


mzPeakFileNameList <- list(
    '20211019/SPCS+d3-Carnitine-positive_processed-peaks.csv',
    '20211019/SPCS+d3-Carnitine+Carnitine(10uM)-positive_processed-peaks.csv',
    '20211019/SPCS+d3-Carnitine+Carnitine(20uM)-positive_processed-peaks.csv',
    '20211019/SPCS+d3-Carnitine+Carnitine(30uM)-positive_processed-peaks.csv',
    '20211019/SPCS+d3-Carnitine+Carnitine(40uM)-positive_processed-peaks.csv',
    c('20211017/SerumExtract-100%EtOH-20210908-positive_processed-peaks.csv',
      '20211017/SerumExtract-100%EtOH-20210908-positive-2_processed-peaks.csv',
      '20211017/SerumExtract-100%EtOH-20210908-positive-3_processed-peaks.csv'),
    c('20211017/SerumExtract-50%CHCl3-50%MeOH-20210908-positive_processed-peaks.csv',
      '20211017/SerumExtract-50%CHCl3-50%MeOH-20210908-positive-2_processed-peaks.csv',
      '20211017/SerumExtract-50%CHCl3-50%MeOH-20210908-positive-3_processed-peaks.csv')
)
mzPeakFileNameList <- addParentDirectory(mzPeakFileNameList, mzDataDir)

mzDestPeakFileName <- c(
    'SPCS+d3-Carnitine.csv',
    'SPCS+d3-Carnitine+Carnitine(10uM).csv',
    'SPCS+d3-Carnitine+Carnitine(20uM).csv',
    'SPCS+d3-Carnitine+Carnitine(30uM).csv',
    'SPCS+d3-Carnitine+Carnitine(40uM).csv',
    'SerumExtract-100%EtOH-20210908+d3-Carnitine.csv',
    'SerumExtract-50%CHCl3-50%MeOH-20210908+d3-Carnitine.csv'
)
mzDestPeakFileName <- addParentDirectory(mzDestPeakFileName, mzStatDir)

mzRanges <- list(
    c(161, 167),
    c(183, 189)
)
mzPeakIntraGroupTolerance <- 0.2
mzPeakIntraGroupRelativeTolerance <- FALSE
# mzPeakIntraGroupMinimumReplicate <- 4/5
mzPeakIntraGroupMinimumReplicate <- 3/3
mzPeakInterGroupTolerance <- 0.1
mzPeakInterGroupRelativeTolerance <- FALSE

mzNormalizationReference <- c(
    165.13185, # d3-Carnitine, [M + H]+
    187.1132 # d3-Carnitine, [M + Na]+
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