# R script for pre-processing 1-D mass spectrum, 
# picking peaks and plotting spectrum & peaks
# Using functions provided by MALDIquant and alsace package

library('MALDIquant')
library('alsace')
source('lib/MSFileIO.R')
source('lib/MSUtils.R')
source('lib/DataFilter.R')


# Function for pre-processing 1-D mass spectrum
# Following steps are carried out in order:
# Smoothing, baseline substraction, background substraction
# sourceFiles: paths of original spectra to be processed
# targetFiles: paths of processed spectra to save to
# backgroundFiles: paths of background spectra to substract
# mzRanges: a list of m/z ranges to consider with; a empty list 
#           will lead to processing of full spectra
# refPeakList: a list of m/z that serve as 'reference', i.e. anchor points
# smoothing: whether or not to do peak smoothing
# baseline: whether or not to do baseline substraction
processMassSpectra <- function(sourceFiles,
                               targetFiles,  
                               backgroundFiles = c(),
                               mzRanges = list(),
                               refPeakList = c(),
                               refPeakSNR = 3,
                               refTolerance = 1E-3,
                               smoothing = TRUE,
                               baseline = TRUE)
{
    # Read m/z intensities from the source spectra
    mzSpectra <- list()
    mzSpectraCount <- c()
    for (source in sourceFiles)
    {
        mzData <- readSpectra(source)
        
        # Deal with files containing only one spectrum
        if (typeof(mzData) != 'list')
            mzData <- list(mzData)
        
        for (mzSpectrum in mzData)
        {
            if (length(mzRanges) > 0)
                mzSpectrum <- mzSpectrum[rangeFilter.index(mzSpectrum[,1], 
                                                           mzRanges),]
            
            mzSpectra <- c(mzSpectra,
                           list(createMassSpectrum(mzSpectrum[,1], 
                                                   mzSpectrum[,2])))
        }
        mzSpectraCount <- c(mzSpectraCount, length(mzData))
    }
    
    # Read background intensities
    bgSpectra <- list()
    if (length(backgroundFiles) > 0)
    {
        for (background in backgroundFiles)
        {
            mzData <- readSpectra(background)
            
            # Deal with files containing only one spectrum
            if (typeof(mzData) != 'list')
                mzData <- list(mzData)
            
            for (mzSpectrum in mzData)
            {
                if (length(mzRanges) > 0)
                    mzSpectrum <- mzSpectrum[rangeFilter.index(mzSpectrum[,1], 
                                                               mzRanges),]
                if (nrow(mzSpectrum) != length(mzSpectra[[1]]@mass))
                {
                    message('The number of features in background spectra 
                             is not equal to that in source spectra!')
                    return()
                }
                bgSpectra <- c(bgSpectra,
                               list(createMassSpectrum(mzSpectrum[,1], 
                                                       mzSpectrum[,2])))
            }
        }
    }
    
    ### BEGIN PROCESSING ###
    # Step 1: Smooth spectra
    if (smoothing)
    {
        mzSpectra <- smoothIntensity(mzSpectra)
        if (length(bgSpectra) > 0)
            bgSpectra <- smoothIntensity(bgSpectra)
    }
    
    # Step 2: Baseline correction
    if (baseline)
    {
        mzSpectra <- removeBaseline(mzSpectra)
        if (length(bgSpectra) > 0)
            bgSpectra <- removeBaseline(bgSpectra)
    }
    
    # Step 3: Align source spectra with reference peaks
    # If no reference peak list is provided, use peaks found in background
    mzPeaks <- detectPeaks(mzSpectra, SNR = refPeakSNR)
    if (length(refPeakList) > 1)
        refPeakList <- lapply(seq(1, length(mzPeaks)), function(i)
        {
            createMassPeaks(refPeakList, 
                            rep(max(mzPeaks[[i]]@intensity), 
                                length(refPeakList)))
        })
    else if (length(bgSpectra) > 0)
        refPeakList <- detectPeaks(bgSpectra, SNR = refPeakSNR)
    else
        refPeakList <- list()
    if (length(refPeakList) > 0)
    {
        # Calculate a reference spectrum
        refPeaks <- referencePeaks(refPeakList)
        
        # Align source spectra
        mzWarpFunction <- determineWarpingFunctions(mzPeaks,
                                                    reference = refPeaks,
        											tolerance = refTolerance,
                                                    allowNoMatches = TRUE)
        suppressWarnings(
            mzSpectra <- warpMassSpectra(mzSpectra, mzWarpFunction)
        )
        
        # Align background spectra
        if (length(refPeakList) > 1 && length(bgSpectra) > 1)
        {
            bgPeaks <- detectPeaks(bgSpectra, SNR = refPeakSNR)
            bgWarpFunction <- determineWarpingFunctions(bgPeaks,
                                                        reference = refPeaks,
            											tolerance = refTolerance,
                                                        allowNoMatches = TRUE)
            suppressWarnings(
                bgSpectra <- warpMassSpectra(bgSpectra, bgWarpFunction)
            )
        }
    }
    
    # Step 4: Substract background intensity from source spectra
    if (length(bgSpectra) > 0)
    {
        # Calculate an average background spectrum
        bgSpectrum <- averageMassSpectra(bgSpectra)
        
        sapply(seq(1, length(mzSpectra)), function(i)
        {
            mzIntensity <- mzSpectra[[i]]@intensity - bgSpectrum@intensity
            mzIntensity[mzIntensity < 0] <- 0
            intensity(mzSpectra[[i]]) <- mzIntensity
        })
    }
    ### END PROCESSING ###
    
    # Write result to target file
    # Convert MALDIquant:MassSpectrum object back to data.frame 
    # for compatibility...
    totalSpectraCount <- 0
    for (i in seq(1, length(mzSpectraCount)))
    {
        for (j in seq(1, mzSpectraCount[i]))
        {
            totalSpectraCount <- totalSpectraCount + 1
            mzFile <- data.frame(mzSpectra[[totalSpectraCount]]@mass, 
                                 mzSpectra[[totalSpectraCount]]@intensity)
            colnames(mzFile) <- c('m/z', 'intensity')
            targetFileName <- targetFiles[i]
            if (j > 1)
                targetFileName <- addFilenameSuffix(targetFileName,
                                                    paste0('_', j))
            saveRDS(mzFile, targetFileName)
        }
    }
}


# Function for picking peaks from spectra, withi given criteria
# sourceFiles: paths of spectra to deal with
# targetFile: paths of CSV files to save peak information to
# mzRanges: a list of m/z ranges to consider with; a empty list 
#           will lead to processing of full spectra
# peakSNR: a signal-to-noise ratio used to filter peaks
# maxPeakWidth: the moving window size (absolute m/z) when fitting a peak
pickPeaks <- function(sourceFiles,
                      targetFiles, 
                      mzRanges = list(),
                      peakSNR = 6,
                      maxPeakWidth = 1)
{
    sapply(seq(1, length(sourceFiles)), function(i)
    {
        # Read source spectrum
        mzSpectrum <- readRDS(sourceFiles[i])
        
        # Deal with format of source file
        # If the source object is a data frame (with at least two columns), 
        # convert it to a MALDIquant::MassSpectrum object
        if (!is.null(nrow(mzSpectrum)) && nrow(mzSpectrum) >= 2)
            mzSpectrum <- createMassSpectrum(mzSpectrum[,1], mzSpectrum[,2])
        
        # Pick peaks with given signal-to-noise ratio
        mzWindowSize <- round(length(mzSpectrum@mass) * maxPeakWidth / 
                              (max(mzSpectrum@mass) - min(mzSpectrum@mass)))
        mzWindowSize <- max(mzWindowSize, 10)
                        
        mzPeaks <- detectPeaks(mzSpectrum, 
                               halfWindowSize = mzWindowSize,
                               SNR = peakSNR)
        
        # Fit found peaks and get relative peak information
        # Only peaks within the given m/z range are choosen
        if (length(mzRanges) > 0)
            mzIndexes <- rangeFilter.index(mzPeaks@mass, mzRanges)
        else
            mzIndexes <- seq(1, length(mzPeaks@mass))
        mzIndexes <- mapUniqueMZ(source = mzPeaks@mass[mzIndexes],
                                 target = mzSpectrum@mass)
        mzIndexes <- mzIndexes[mzIndexes > 1 & 
                               mzIndexes < length(mzSpectrum@mass)]
        mzPeakInfo <- fitpeaks(mzSpectrum@intensity, mzIndexes)
        
        # Eliminating peaks whose SNR is smaller than given threshold
        # Remove also lines with NA
        # Assuming the first column is a list of m/z, 
        # and the fourth column is a list of intensities
        mzPeakInfoColNames <- colnames(mzPeakInfo)
        mzPeakMasks <- !is.na(mzPeakInfo[,1])
        mzIndexes <- mzPeakInfo[mzPeakMasks, 1]
        mzNoiseWindowSize <- min(mzWindowSize * 10 / length(mzSpectrum@mass), 1)
        mzNoises <- 
            estimateNoise(mzSpectrum, method = 'SuperSmoother',
                          span = mzNoiseWindowSize)
        mzPeakInfo <- 
            mzPeakInfo[mzPeakMasks &
                       mzPeakInfo[,4] >= mzNoises[mzIndexes, 2] * peakSNR,]
        
        # Ugly hack: correct dimension of mzPeakInfo and restore column names
        # when mzPeakInfo contains only one entry (i.e. one peak found)
        if (length(mzPeakInfo) == length(mzPeakInfoColNames))
        {
            dim(mzPeakInfo) <- c(1, length(mzPeakInfo))
            colnames(mzPeakInfo) <- mzPeakInfoColNames
        }
        
        # Convert m/z index to m/z
        mzList <- mzSpectrum@mass
        mzPeakInfo[,1] <- mzList[mzPeakInfo[,1]] # centroid
        mzInterval <- (max(mzList) - min(mzList)) / (length(mzList) - 1)
        mzPeakInfo[,2] <- mzPeakInfo[,2] * mzInterval # standard deviation
        mzPeakInfo[,3] <- mzPeakInfo[,3] * mzInterval # FWHM
        
        # Rename columns
        colnames(mzPeakInfo)[1] <- 'm/z'
        colnames(mzPeakInfo)[2] <- 'standard deviation'
        colnames(mzPeakInfo)[3] <- 'FWHM'
        
        # Save peaks information
        write.csv(mzPeakInfo, targetFiles[i], row.names = FALSE)
    })
}


# Function for plotting spectra files and save as images
# spectraFileNames: paths of spectra to plot
# peakFileNames: paths of CSV files containing peak information 
#                to append to each spectrum
# removeMonoisotopic: Whether or not to find monoisotopic patterns
#                     and to plot only the most abundant peaks
# plotFileNames: paths of image files to save to
# titles: Titles for each spectrum
# mzRange: a range of m/z used as the boundary of x axis
# mzTick: intervals of m/z labels; use 0 for auto adjusting
# intensityRange: a range of intensity used as the boundary of y axis
# peakSNR: a signal-to-noise ratio used for filtering when labelling peaks
#          a value less or equal to 0 will lead all peaks to be labelled
# maxPeakWidth: the moving window size (absolute m/z) when fitting a peak
# labelSize: size of the annotation text over peaks
# resolution: image resolution (width and height)
plotSpectra <- function(spectraFileNames, 
                        peakFileNames = c(),
                        plotFileNames = c(),
                        removeMonoisotopic = FALSE,
                        titles = c(),
                        mzRange = c(NA, NA),
                        mzTick = 0,
                        intensityRange = c(),
                        peakSNR = 0,
                        maxPeakWidth = 1,
                        labelSize = 1,
                        resolution = c(2000, 1000))
{
    
    # Use blank text as default title
    if (length(titles) < length(plotFileNames))
        titles <- c(titles, rep('', length(plotFileNames) - length(titles)))
    
    sapply(seq(1, length(spectraFileNames)), function(i)
    {
        mzSpectrum <- readRDS(spectraFileNames[i])
        
        # Deal with MALDIquant::MassSpectrum object
        massList <- c()
        intensityList <- c()
        if (class(mzSpectrum)[] == 'MassSpectrum')
        {
            massList <- mzSpectrum@mass
            intensityList <- mzSpectrum@intensity
        }
        else
        {
            massList <- mzSpectrum[,1]
            intensityList <- mzSpectrum[,2]
        }
        
        # Deal with NA values in plot range of m/z (x axis)
        if (is.na(mzRange[1]))
            plotMZMin <- round(min(massList), max(-log10(mzTick), 1E-6))
        else
            plotMZMin <- mzRange[1]
        if (is.na(mzRange[2]))
            plotMZMax <- round(max(massList), max(-log10(mzTick), 1E-6))
        else
            plotMZMax <- mzRange[2]
        
        # Calculate label interval for m/z (x axis) if necessary
        if (mzTick <= 0)
        {
            plotMZScale <- pretty(massList[massList >= plotMZMin & 
                                           massList <= plotMZMax])
            mzTick <- plotMZScale[2] - plotMZScale[1]
        }
        
        # Calculate min/max and label scale for intensity (y axis)
        if (length(intensityRange) < 2)
        {
            plotIntMin <- 0
            if (length(intensityRange) < 1)
            {
                mzIndexes <- massList > plotMZMin & massList < plotMZMax
                plotIntMax <- max(pretty(intensityList[mzIndexes]))
            }
            else
                plotIntMax <- intensityRange
        }
        else
        {
            plotIntMin <- intensityRange[1]
            plotIntMax <- intensityRange[2]
        }
        
        # Read peak list from file
        peakList <- NULL
        if (i <= length(peakFileNames))
        {
            # Assuming the first column is a list of m/z,
            # and the fourth column is a list of peak intensities
            peakFile <- read.csv(peakFileNames[i], as.is = TRUE)
            peakList <- createMassPeaks(peakFile[,1], peakFile[,4])
            
            # Remove monoisotopic peaks if necessary
            if (removeMonoisotopic)
            {
                peakList <- monoisotopicPeaks(peakList)
            }
            
            # Filter peaks with given SNR
            if (peakSNR > 0)
            {
                mzWindowSize <- max(round(length(massList) * maxPeakWidth / 
                                          (max(massList) - min(massList))), 10)
                mzNoiseWindowSize <- min(mzWindowSize * 10 / length(massList), 
                                         1)
                mzNoiseLevel <- estimateNoise(
                                    createMassSpectrum(massList, intensityList), 
                                    method = 'SuperSmoother',
                                    span = mzNoiseWindowSize)
                mzIndexes <- mapUniqueMZ(source = peakList@mass, 
                                         target = massList)
                peakFilter <- peakList@intensity > 
                                    mzNoiseLevel[mzIndexes, 2] * peakSNR
                peakFilter[is.na(peakFilter)] <- FALSE
                peakList <- createMassPeaks(peakFile[peakFilter, 1], 
                                            peakFile[peakFilter, 4])
            }
        }
        
        # Plot the spectrum, then label the peaks
        if (!is.null(plotFileNames[i]) && !is.na(plotFileNames[i]))
            saveToFile <- TRUE
        else
            saveToFile <- FALSE
        if (saveToFile)
            png(plotFileNames[i], 
                width = resolution[1], 
                height = resolution[2],
                res = 100)
        plot(massList, intensityList, 
             type = 'h', xaxt = 'n', xaxs = 'i',
             xlim = c(plotMZMin, plotMZMax),
             ylim = c(plotIntMin, plotIntMax),
             xlab = expression(italic('m/z')),
             ylab = 'Intensity')
        axis(1, at = seq(plotMZMin, plotMZMax, mzTick))
        if (!is.null(peakList) && length(peakList@mass) > 0)
            labelPeaks(peakList, cex = labelSize)
        title(titles[i])
        if (saveToFile)
            dev.off()
    })
}