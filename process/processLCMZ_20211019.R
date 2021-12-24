# R script for pre-processing 2-D LC-MS spectra
# Using executables from ThermoRawFileParser project
# Using functions provided by mzR package and MALDIquant package

library('mzR')
source('lib/PathUtils.R')
source('lib/MassSpectraProcessor.R')


mzDataDir <- '../data/spectra/HPLC-MS'

mzOrigFileName <- c(
    '20211019/SPCS+d3-Carnitine-positive.mzML',
    '20211019/SPCS+d3-Carnitine+Carnitine(10uM)-positive.mzML',
    '20211019/SPCS+d3-Carnitine+Carnitine(20uM)-positive.mzML',
    '20211019/SPCS+d3-Carnitine+Carnitine(30uM)-positive.mzML',
    '20211019/SPCS+d3-Carnitine+Carnitine(40uM)-positive.mzML'
)
mzOrigFileName <- addParentDirectory(mzOrigFileName, mzDataDir)

mzDestFileName <- changeFilenameSuffix(mzOrigFileName, 
                                       '_processed.rds', '.mzML')

mzPeakFileName <- changeFilenameSuffix(mzDestFileName, '-peaks.csv', '.rds')

mzExtractSpectrumIDList <- c(
    # TIC(max, positive) for 5 uL + 200 uL/min
    # With column
    # RT = 74.58 ~ 94.48 sec; Scan = 321 ~ 361 (Carnitine/d3-Carnitine@STD)
    rep(list(c(335:424)), 5)
)

mzRanges <- list(
)

mzPeakResolution <- 1e-3

mzPeakQuantificationSNR <- 2  # General quantification
mzPeakMaxWidth <- 0.3


# Utility functions
selectDataByIndexes <- function(index, data, indexList, threshold, algorithm, ...)
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


### MAIN ENTRY ###
mzFileNames <- c()
for (i in seq(1, length(mzOrigFileName)))
{
    # Extract spectra (with given IDs) for each spectra file    
    spectrumFile <- openMSfile(mzOrigFileName[i])
    mzPeaks <- peaks(spectrumFile, mzExtractSpectrumIDList[[i]])
    
    # Calculate an averaged spectrum
    mzList <- lapply(mzPeaks, function(X){ return(X[,1]) })
    intensityList <- lapply(mzPeaks, function(X){ return(X[,2]) })
    mzClusterIndexes <- clusterMZ(mzList, mzPeakResolution, TRUE)
    allIndexes <- unique(unlist(mzClusterIndexes))
    averagedMZ <- unlist(lapply(allIndexes,
                                data = mzList,
                                indexList = mzClusterIndexes,
                                threshold = 0,
                                algorithm = mean,
                                selectDataByIndexes))
    averagedIntensities <- unlist(sapply(allIndexes,
                                         data = intensityList,
                                         indexList = mzClusterIndexes,
                                         threshold = 0,
                                         count = length(intensityList),
                                         algorithm = function(X, count)
                                         {
                                             return(sum(X) / count)
                                         },
                                         selectDataByIndexes))
    averagedSpectrum <- cbind(averagedMZ, averagedIntensities)
    colnames(averagedSpectrum) <- c('m/z',' intensity')
    
    # Save the result as a temporary file
    mzFileName <- tempfile(fileext = '.rds')
    saveRDS(averagedSpectrum, mzFileName)
    mzFileNames <- c(mzFileNames, mzFileName)
}

processMassSpectra(mzFileNames, mzDestFileName,
                   mzRanges = mzRanges,
                   baseline = FALSE, smoothing = FALSE)

# Peak picking with specified criteria
pickPeaks(mzDestFileName, mzPeakFileName, 
          mzRanges, mzPeakQuantificationSNR, mzPeakMaxWidth)