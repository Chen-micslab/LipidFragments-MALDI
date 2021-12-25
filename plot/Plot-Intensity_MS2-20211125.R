# R script for plotting 1-D mass spectra for publication

library('ggplot2')
library('gridExtra')
source('lib/PathUtils.R')
source('lib/MSFileIO.R')
source('lib/MSUtils.R')


mzDataDir <- '../data/spectra/MALDI-TOF'
mzPlotDir <- '../data/plot'

mzSpectrumFilenameList <- list(
    c('2021-11-14/TiO2NP25+NaCl+PE(18-18)-MS2/0_O20/1/111.1000.LIFT/Spectrum_processed.rds'),
    c('2021-11-08/TiO2NP25+NaCl+PE(18-18)-MS2/0_E14/1/164.0000.LIFT/Spectrum_processed.rds'),
    c('2021-11-24/TiO2NP25+NaCl+PE(18-18)-MS2/0_B4/1/173.0000.LIFT/Spectrum_processed.rds'),
    c('2021-11-24/TiO2NP25+NaCl+PE(18-18)-MS2/0_B5/1/339.2000.LIFT/Spectrum_processed.rds'),
    c('2021-11-24/TiO2NP25+NaCl+PE(18-18)-MS2/0_D1/1/463.2000.LIFT/Spectrum_processed.rds'),
    c('2021-11-08/TiO2NP25+NaCl+PE(18-18)-MS2/0_E15/1/603.5000.LIFT/Spectrum_processed.rds')
    # c('2021-11-14/TiO2NP25+NaCl+PG(18(0)-18(1))-MS2/0_P20/1/111.1000.LIFT/Spectrum_processed.rds'),
    # c('2021-11-14/TiO2NP25+NaCl+PG(18(0)-18(1))-MS2/0_P22/1/123.1000.LIFT/Spectrum_processed.rds'),
    # c('2021-11-24/TiO2NP25+NaCl+PG(18(0)-18(1))-MS2/0_D4/1/135.1000.LIFT/Spectrum_processed.rds'),
    # c('2021-11-24/TiO2NP25+NaCl+PG(18(0)-18(1))-MS2/0_D4/1/173.0000.LIFT/Spectrum_processed.rds'),
    # c('2021-11-24/TiO2NP25+NaCl+PG(18(0)-18(1))-MS2/0_D4/1/177.0000.LIFT/Spectrum_processed.rds'),
    # c('2021-11-08/TiO2NP25+NaCl+PG(18(0)-18(1))-MS2/0_F14/1/195.0000.LIFT/Spectrum_processed.rds'),
    # c('2021-11-24/TiO2NP25+NaCl+PG(18(0)-18(1))-MS2/0_D4/1/217.0000.LIFT/Spectrum_processed.rds'),
    # c('2021-11-08/TiO2NP25+NaCl+PG(18(0)-18(1))-MS2/0_F17/1/605.5000.LIFT/Spectrum_processed.rds')
)        
mzSpectrumFilenameList <- addParentDirectory(mzSpectrumFilenameList, mzDataDir)

mzPeakFileNameList <- changeFilenameSuffix(mzSpectrumFilenameList,
                                           '-peaks.csv', '.rds')

mzBackgroundPeakFilenameList <- list(
    c('2021-11-14/TiO2NP25+EtOH+NaCl-MS2/0_O19/1/111.1000.LIFT/Spectrum_processed-peaks.csv'),
    c('2021-11-08/TiO2NP25+EtOH+NaCl-MS2/0_C14/1/164.0000.LIFT/Spectrum_processed-peaks.csv'),
    c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/173.0000.LIFT/Spectrum_processed-peaks.csv'),
    c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/339.2000.LIFT/Spectrum_processed-peaks.csv'),
    c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/463.2000.LIFT/Spectrum_processed-peaks.csv'),
    c('2021-11-08/TiO2NP25+EtOH+NaCl-MS2/0_C15/1/603.5000.LIFT/Spectrum_processed-peaks.csv')
    # c('2021-11-14/TiO2NP25+EtOH+NaCl-MS2/0_O19/1/111.1000.LIFT/Spectrum_processed-peaks.csv'),
    # c('2021-11-14/TiO2NP25+EtOH+NaCl-MS2/0_O19/1/123.1000.LIFT/Spectrum_processed-peaks.csv'),
    # c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/135.1000.LIFT/Spectrum_processed-peaks.csv'),
    # c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/173.0000.LIFT/Spectrum_processed-peaks.csv'),
    # c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/177.0000.LIFT/Spectrum_processed-peaks.csv'),
    # c('2021-11-08/TiO2NP25+EtOH+NaCl-MS2/0_C14/1/195.0000.LIFT/Spectrum_processed-peaks.csv'),
    # c('2021-11-24/TiO2NP25+EtOH+NaCl-MS2/0_D5/1/217.0000.LIFT/Spectrum_processed-peaks.csv'),
    # c('2021-11-08/TiO2NP25+EtOH+NaCl-MS2/0_C17/1/605.5000.LIFT/Spectrum_processed-peaks.csv')
)
mzBackgroundPeakFilenameList <- addParentDirectory(mzBackgroundPeakFilenameList, mzDataDir)

mzPlotFilename <- 'TiO2NP25+PE(18-18)-MS2.png'
# mzPlotFilename <- 'TiO2NP25+PG(18-18)-MS2.png'
mzPlotFilename <- addParentDirectory(mzPlotFilename, mzPlotDir)

mzSpectrumNameList <- list(
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'
)

mzRange <- c(20, 1000)

mzPeakResolution <- 5e-5  # Average resolution of the detector (MS2 at m/z = 200)
mzPeakTolerance <- 0.1
mzRelativePeakTolerance <- FALSE

mzPlotDimensionX <- 2
mzPlotDimensionY <- 3  # PE(18-18)
# mzPlotDimensionY <- 4  # PG(18-18)
mzPlotSubDimensionX <- 1
mzPlotSubDimensionY <- 1
mzPlotDirection <- 'v'

mzPlotType <- c(
    rep('h', 6)  # PE(18-18)
    # rep('h', 8)  # PG(18-18)
)

# mzPlotLabelAxisX <- expression(italic('m') * '/' * italic('z'))
mzPlotLabelAxisX <- ''
# mzPlotLabelAxisY <- 'Relative Intensity (%)'
mzPlotLabelAxisY <- ''

mzPlotRangeListX <- list(
    c(20, 105),  # PE(18-18), m/z 111
    c(20, 160),  # PE(18-18), m/z 164
    c(20, 165),  # PE(18-18), m/z 173
    c(20, 335),  # PE(18-18), m/z 339
    c(20, 455),  # PE(18-18), m/z 463
    c(20, 600)  # PE(18-18), m/z 603
    # c(20, 105),  # PG(18-18), m/z 111
    # c(20, 115),  # PG(18-18), m/z 123
    # c(20, 130),  # PG(18-18), m/z 135
    # c(20, 170),  # PG(18-18), m/z 173
    # c(20, 170),  # PG(18-18), m/z 177
    # c(20, 190),  # PG(18-18), m/z 195
    # c(20, 210),  # PG(18-18), m/z 217
    # c(20, 600)  # PG(18-18),m/z 605
)
mzPlotRangeY <- c(0, 1.2)

mzPlotAxisTickX <- NULL
mzPlotAxisTickY <- seq(0, 1, 0.25)

mzPlotFontSize <- 50
mzPlotFontSizeTitle <- 20
mzPlotFontSizeAxisLabel <- 30
mzPlotFontSizeAnnotation <- 6

mzPlotLineColorList <- c(rep(list('#5F0000'), 8))
mzPlotAnnotationColor <- '#000000'

mzPlotAnnotationText <- '○'
mzPlotAnnotationThreshold <- (mzPlotRangeY[2] - mzPlotRangeY[1]) * 0.01
mzPlotAnnotationPos <- function(X)
{
    offset <- (mzPlotRangeY[2] - mzPlotRangeY[1]) * 0.05
    max <- mzPlotRangeY[2] - offset
    pos <- X + offset
    pos[pos > max] <- max
    return(pos)
}

mzPlotDevice <- png
mzPlotResolution <- c(1950, 1200)  # PE(18-18)
# mzPlotResolution <- c(1950, 1600)  # PG(18-18)


### MAIN ENTRY ###
# Collect the spectrum data
mzDataList <- lapply(mzSpectrumFilenameList, 
                     function(Y){ lapply(Y, function(X){ readSpectrum(X) }) })
mzPeakList <- lapply(mzPeakFileNameList, 
                     function(Y){ lapply(Y, 
                        function(X){ read.csv(X, as.is = TRUE)[,c(1,4)] }) })
mzBackgroundList <- lapply(mzBackgroundPeakFilenameList, 
                           function(Y){ lapply(Y, 
                               function(X){ read.csv(X, as.is = TRUE) }) })

# Filter the m/z values
mzDataList <- lapply(mzDataList,
                     min = mzRange[1],
                     max = mzRange[2],
                     function(data, min, max) { lapply(data,
                             function(X){ X[X[,1] >= min & X[,1] <= max,] })})

# Normalize the intensities (max peak normalization)
mzDataList <- lapply(mzDataList, 
                     function(data){ lapply(data, 
                        function(X)
                        {
                            X[,2] <- X[,2] / max(X[,2], na.rm = T)
                            return(X)
                        })})
mzPeakList <- lapply(mzPeakList, 
                     function(data){ lapply(data, 
                        function(X)
                        {
                            X[,2] <- X[,2] / max(X[,2], na.rm = T)
                            return(X)
                        })})

plotList <- list()
for (i in seq(1, length(mzDataList)))
{
    
    # Construct the data for plotting figures    
    plotData <- data.frame()
    plotAnnotationData <- data.frame(matrix(nrow = 0, ncol = 3))
    for (j in seq(1, length(mzDataList[[i]])))
    {
        plotData <- rbind(plotData,
                          cbind(mzDataList[[i]][[j]], j))
        if (j <= length(mzPeakList[[i]]) && j <= length(mzBackgroundList[[i]]))
        {
            # Calculate the position of annotation text
            mzIndexes <- mapUniqueMZ(mzBackgroundList[[i]][[j]][,1], 
                                     mzPeakList[[i]][[j]][,1], 
                                     tolerance = mzPeakTolerance,
                                     relativeTolerance = mzRelativePeakTolerance)
            
            # Exclude weak peaks
            mzIndexes <- mzIndexes[mzPeakList[[i]][[j]][mzIndexes, 2] >=
                                   mzPlotAnnotationThreshold]
            # Exclude peaks of general isotopic patterns (2H, 13C, etc.)
            # Always annotate the first peak
            firstPeak <- mzIndexes == 1
            mzIndexes <- mzIndexes[!firstPeak]
            mzIndexes <- mzIndexes[mzPeakList[[i]][[j]][mzIndexes - 1, 1] + 1.05 < 
                                   mzPeakList[[i]][[j]][mzIndexes, 1] | 
                                   mzPeakList[[i]][[j]][mzIndexes - 1, 2]  <= 
                                   mzPeakList[[i]][[j]][mzIndexes, 2]]
            mzIndexes <- c(rep(1, length(firstPeak[firstPeak])), mzIndexes)
            if (length(mzIndexes) == 0)
                next()
            
            mzIndexes <- mapUniqueMZ(mzPeakList[[i]][[j]][mzIndexes, 1], 
                                     mzDataList[[i]][[j]][,1],
                                     tolerance = mzPeakResolution, 
                                     relativeTolerance = FALSE)
            plotAnnotationData <- rbind(plotAnnotationData,
                                        cbind(mzDataList[[i]][[j]][mzIndexes, 1],
                                              mzPlotAnnotationPos(
                                                  mzDataList[[i]][[j]][mzIndexes, 2]), 
                                              j))
        }
    }
    colnames(plotData) <- c('mz', 'intensity', 'id')
    colnames(plotAnnotationData) <- c('x', 'y', 'id')
    
    # Plot figures
    p <- ggplot(plotData, aes(x = mz, y = intensity, group = id))
    if (mzPlotType[i] == 'h')
        p <- p + geom_histogram(aes(color = as.factor(id)), stat = 'identity', 
                                show.legend = FALSE)
    else
        p <- p + geom_line(aes(color = as.factor(id)), show.legend = FALSE)
    if (!is.null(mzPlotAxisTickX))
        p <- p + scale_x_continuous(breaks = mzPlotAxisTickX)
    plotScale <- c()
    if (!all(is.na(mzPlotRangeListX[[i]])))
        plotScale <- 'free_y'
    if (!all(is.na(mzPlotRangeY)))
        plotScale <- c(plotScale, 'free_x')
    if (length(plotScale) == 0)
        plotScale <- 'free'
    if (length(plotScale) > 1)
        plotScale <- 'fixed'
    p <- p + 
        scale_y_continuous(breaks = mzPlotAxisTickY,
                           labels = function(X)
                           { formatC(X * 100, digits = 0, format = 'f') }) +
        scale_color_manual(values = mzPlotLineColorList[[i]]) + 
        facet_wrap(~id, nrow = mzPlotSubDimensionY, ncol = mzPlotSubDimensionX,
                   dir = mzPlotDirection,
                   scales = plotScale, shrink = TRUE)
    if (any(is.na(mzPlotRangeListX[[i]])))
        mzPlotRangeListX[[i]] <- NULL
    if (any(is.na(mzPlotRangeY)))
        mzPlotRangeY <- NULL
    p <- p + 
        coord_cartesian(xlim = mzPlotRangeListX[[i]], 
                        ylim = mzPlotRangeY, expand = FALSE)
    if (nrow(plotAnnotationData) > 0)
        p <- p + geom_text(plotAnnotationData, 
                           mapping = aes(x = x, y = y, group = id,
                                         label = mzPlotAnnotationText),
                           color = mzPlotAnnotationColor,
                           size = mzPlotFontSizeAnnotation,
                           show.legend = FALSE)
    p <- p + theme_bw() + 
        theme(text = element_text(size = mzPlotFontSize),
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              axis.title.y = element_text(margin = margin(r = 0.5, unit = 'line')),
              axis.title.x = element_text(margin = margin(t = 0.5, unit = 'line')),
              axis.text = element_text(size = mzPlotFontSizeAxisLabel),
              axis.ticks.length = unit(mzPlotFontSizeAxisLabel / 2, 'point'),
              panel.spacing.y = unit(0, 'line'),
              # plot.margin = margin(t = 2, r = 4, b = 2, l = 2, unit = 'line')
              plot.margin = margin(t = 0, r = 2, b = -2, l = -4, unit = 'line')
        )
    p <- p + 
        xlab(mzPlotLabelAxisX) + ylab(mzPlotLabelAxisY)
    
    # Add title for subplots
    plotText <- data.frame(x = -Inf, y = Inf,
                           hjust = -1.0, vjust = 1.5,
                           label = mzSpectrumNameList[[i]],
                           id = unique(plotData$id))
    p <- p + geom_text(plotText, 
                       mapping = aes(x = x, y = y, 
                                     hjust = hjust, vjust = vjust,
                                     label = label),
                       size = mzPlotFontSizeTitle)
    
    plotList <- c(plotList, list(p))
}

# Combine multiple plots into a large one
mplot <- arrangeGrob(grobs = plotList,
                     nrow = mzPlotDimensionY, ncol = mzPlotDimensionX)

# Save the plot
ggsave(mzPlotFilename, mplot, device = mzPlotDevice,
       width = mzPlotResolution[1], height = mzPlotResolution[2], 
       limitsize = FALSE)
