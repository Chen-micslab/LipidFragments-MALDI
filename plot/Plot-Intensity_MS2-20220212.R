# R script for plotting 1-D mass spectra for publication

library('ggplot2')
library('gridExtra')
source('lib/PathUtils.R')
source('lib/MSFileIO.R')
source('lib/MSUtils.R')


mzDataDir <- '../data/spectra/MALDI-TOF'
mzPlotDir <- '../data/plot'

mzSpectrumFilenameList <- list(
    c('2022-02-11/TiO2NP25+NaCl+PC(18(1)-18(1))-MS2/0_M2/1/104.1000.LIFT/Spectrum_processed.rds',
      '2021-11-25/TiO2NP25+SerumExtract-100%EtOH-20210908-MS2/0_G7/1/104.1000.LIFT/Spectrum_processed.rds'),
    c('2022-02-12/TiO2NP25+NaCl+PC(18(1)-18(1))-MS2/0_B7/1/147.0000.LIFT/Spectrum_processed.rds',
      '2021-12-21/TiO2NP25+SerumExtract-100%EtOH-20210908-MS2/0_E22/1/147.0000.LIFT/Spectrum_processed.rds'),
    c('2022-02-10/TiO2NP25+NaCl+PC(18(1)-18(1))-MS2/0_K4/1/184.1000.LIFT/Spectrum_processed.rds',
      '2021-11-25/TiO2NP25+SerumExtract-100%EtOH-20210908-MS2/0_G10/1/184.1000.LIFT/Spectrum_processed.rds'),
    c('2022-02-10/TiO2NP25+NaCl+PC(18(1)-18(1))-MS2/0_K5/1/198.1000.LIFT/Spectrum_processed.rds',
      '2021-11-15/TiO2NP25+SerumExtract-100%EtOH-20210908-MS2/0_P14/1/198.1000.LIFT/Spectrum_processed.rds')
)        
mzSpectrumFilenameList <- addParentDirectory(mzSpectrumFilenameList, mzDataDir)

mzPeakFileNameList <- changeFilenameSuffix(mzSpectrumFilenameList,
                                           '-peaks.csv', '.rds')

mzBackgroundPeakFilenameList <- list(
    c('2022-02-11/TiO2NP25+EtOH+NaCl-MS2/0_O1/1/104.1000.LIFT/Spectrum_processed-peaks.csv',
      '2021-11-25/TiO2NP25+SPCS-100%EtOH-20210908-MS2/0_F7/1/104.1000.LIFT/Spectrum_processed-peaks.csv'),
    c('2022-02-12/TiO2NP25+EtOH+NaCl-MS2/0_A7/1/147.0000.LIFT/Spectrum_processed-peaks.csv',
      '2021-12-21/TiO2NP25+SPCS-100%EtOH-20210908-MS2/0_D22/1/147.0000.LIFT/Spectrum_processed-peaks.csv'),
    c('2022-02-10/TiO2NP25+EtOH+NaCl-MS2/0_J4/1/184.1000.LIFT/Spectrum_processed-peaks.csv',
      '2021-11-25/TiO2NP25+SPCS-100%EtOH-20210908-MS2/0_F10/1/184.1000.LIFT/Spectrum_processed-peaks.csv'),
    c('2022-02-10/TiO2NP25+EtOH+NaCl-MS2/0_J5/1/198.1000.LIFT/Spectrum_processed-peaks.csv',
      '2021-11-15/TiO2NP25+SPCS-100%EtOH-20210908-MS2/0_N14/1/198.1000.LIFT/Spectrum_processed-peaks.csv')
)
mzBackgroundPeakFilenameList <- addParentDirectory(mzBackgroundPeakFilenameList, mzDataDir)

mzPlotFilename <- 'TiO2NP25+PC(18-18)+SE1-MS2.png'
mzPlotFilename <- addParentDirectory(mzPlotFilename, mzPlotDir)

mzSpectrumNameList <- list(
    c('A', ''), c('B', ''), c('C', ''), c('D', '')
)

mzRange <- c(20, 1000)

mzPeakResolution <- 5e-5  # Average resolution of the detector (MS2 at m/z = 200)
mzPeakTolerance <- 0.15
mzRelativePeakTolerance <- FALSE

mzPlotDimensionX <- 4
mzPlotDimensionY <- 1
mzPlotSubDimensionX <- 1
mzPlotSubDimensionY <- 2
mzPlotDirection <- 'v'

mzPlotType <- c(
    rep('h', 4)
)

mzPlotLabelAxisX <- ''
# mzPlotLabelAxisX <- expression(italic('m') * '/' * italic('z'))
mzPlotLabelAxisY <- ''
# mzPlotLabelAxisY <- 'Relative Intensity (%)'

mzPlotRangeListX <- list(
    c(20, 100),  # m/z 104
    c(20, 140),  # m/z 147
    c(20, 180),  # m/z 184
    c(20, 195)  # m/z 198
)
mzPlotRangeY <- c(0, 1.1)
# mzPlotRangeY <- c(0, 0.2)

mzPlotAxisTickX <- NULL

mzPlotLineSize <- 1
mzPlotFontSize <- 100
mzPlotFontSizeTitle <- 40
mzPlotFontSizeAxisLabel <- 60
mzPlotFontSizeAnnotation <- 12

mzPlotLineColorList <- c(rep(list(c('#00A0FF', '#5F0000')), 4))
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
mzPlotResolution <- c(4800, 1900)


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
    plotScale <- c()
    if (!all(is.na(mzPlotRangeListX[[i]])))
        plotScale <- 'free_y'
    if (!all(is.na(mzPlotRangeY)))
        plotScale <- c(plotScale, 'free_x')
    if (length(plotScale) == 0)
        plotScale <- 'free'
    if (length(plotScale) > 1)
        plotScale <- 'fixed'
    p <- ggplot(plotData, aes(x = mz, y = intensity, group = id))
    if (mzPlotType[i] == 'h')
        p <- p + geom_histogram(aes(color = as.factor(id)), stat = 'identity', 
                                size = mzPlotLineSize, 
                                show.legend = FALSE)
    else
        p <- p + geom_line(aes(color = as.factor(id)), size = mzPlotLineSize, 
                           show.legend = FALSE)
    if (!is.null(mzPlotAxisTickX))
        p <- p + scale_x_continuous(breaks = mzPlotAxisTickX)
    p <- p + 
        scale_y_continuous(labels = function(X)
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
              plot.margin = margin(t = 0, r = 4, b = 0, l = -4, unit = 'line')
        )
    p <- p + 
        xlab(mzPlotLabelAxisX) + ylab(mzPlotLabelAxisY)
    
    # Add title for subplots
    plotText <- data.frame(x = -Inf, y = Inf,
                           hjust = -0.8, vjust = 2.5,
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
