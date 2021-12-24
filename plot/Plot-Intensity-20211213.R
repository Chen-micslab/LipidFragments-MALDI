# R script for plotting 1-D mass spectra for publication

library('ggplot2')
source('lib/PathUtils.R')
source('lib/MSFileIO.R')
source('lib/MSUtils.R')


mzDataDir <- '../data/spectra/MALDI-TOF'
mzStatDir <- '../data/summary/MALDI-TOF'
mzPlotDir <- '../data/plot'

mzSpectrumFilename <- c(
    '2021-11-09/TiO2NP25+NaCl+PC(18(0)-18(2))/0_D20/1/Spectrum_processed.rds',
    '2021-12-13/TiO2NP25+SerumExtract-100%EtOH-20210908/0_F22/1/Spectrum_processed.rds'
)
mzSpectrumFilename <- addParentDirectory(mzSpectrumFilename, mzDataDir)

mzPeakFileName <- changeFilenameSuffix(mzSpectrumFilename,
                                       '-peaks.csv', '.rds')

mzAnnotationPeakFilename <- c(
    'TiO2NP25+NaCl+PC(18:0-18:2)_noBG-selected.csv',
    'TiO2NP25+SerumExtract-100%EtOH-20210908-selected.csv'
)
mzAnnotationPeakFilename <- addParentDirectory(mzAnnotationPeakFilename, mzStatDir)

mzPlotFilename <- 'TiO2NP25+PC(18-18)+SE1.png'
# mzPlotFilename <- 'TiO2NP25+SE1-mz125.png'
# mzPlotFilename <- 'TiO2NP25+SE1-mz166.png'
mzPlotFilename <- addParentDirectory(mzPlotFilename, mzPlotDir)

mzSpectrumName <- c(
    'A', 'B'
    # '', ''
)

mzRange <- c(50, 202)

mzPeakResolution <- 2e-5  # Average resolution of the detector (at m/z = 500)
mzPeakTolerance <- 0.05
mzRelativePeakTolerance <- FALSE

mzPlotDimensionX <- 1
mzPlotDimensionY <- 2

mzPlotLabelAxisX <- expression(italic('m') * '/' * italic('z'))
mzPlotLabelAxisY <- 'Normalized Intensity (%)'

mzPlotRangeX <- c(69, 202)
# mzPlotRangeX <- c(123.8, 126.3)
# mzPlotRangeX <- c(165, 167.5)
mzPlotRangeY <- c(0, 1.11)
# mzPlotRangeY <- c(0, 0.01)

mzPlotAxisTickX <- seq(50, 200, 10)
# mzPlotAxisTickX <- seq(124, 126, 1)
# mzPlotAxisTickX <- seq(165, 168, 1)

mzPlotFontSize <- 50
mzPlotFontSizeTitle <- 20
mzPlotFontSizeAxisLabel <- 30
mzPlotFontSizeAnnotation <- 7
# mzPlotFontSizeAnnotation <- 16

mzPlotLineColor <- c('#00A0FF', 
                     '#5F0000')
mzPlotAnnotationColor <- '#000000'

mzPlotAnnotationText <- '▽'
# mzPlotAnnotationText <- ''
mzPlotAnnotationUseX <- TRUE
mzPlotAnnotationThreshold <- (max(mzPlotRangeY) - min(mzPlotRangeY)) * 0.01
mzPlotAnnotationPos <- function(X)
{
    offset <- (mzPlotRangeY[2] - mzPlotRangeY[1]) * 0.05
    max <- mzPlotRangeY[2] - offset
    pos <- X + offset
    pos[pos > max] <- max
    return(pos)
}

mzPlotDevice <- png
mzPlotResolution <- c(2000, 1000)
# mzPlotResolution <- c(800, 400)


### MAIN ENTRY ###
# Collect the spectrum data
mzDataList <- lapply(mzSpectrumFilename, function(X){ readSpectrum(X)} )
mzPeakList <- lapply(mzPeakFileName, 
                     function(X){ read.csv(X, as.is = TRUE)[,c(1,4)] })
mzAnnotationList <- lapply(mzAnnotationPeakFilename, 
                           function(X){ read.csv(X, as.is = TRUE) })

# Filter the m/z values
mzDataList <- lapply(mzDataList,
                     min = mzRange[1],
                     max = mzRange[2],
                     function(X, min, max){ X[X[,1] >= min & X[,1] <= max,] })

# Normalize the intensities (max peak normalization)
mzDataList <- lapply(mzDataList, function(X)
                     {
                         X[,2] <- X[,2] / max(X, na.rm = T)
                         return(X)
                     })
mzPeakList <- lapply(mzPeakList, function(X)
                     {
                         X[,2] <- X[,2] / max(X, na.rm = T)
                         return(X)
                     })

# Construct the data for plotting figures
plotData <- data.frame()
plotAnnotationData <- data.frame(matrix(nrow = 0, ncol = 3))
for (i in seq(1, length(mzDataList)))
{
    plotData <- rbind(plotData,
                      cbind(mzDataList[[i]], i))
    if (i <= length(mzPeakList) && i <= length(mzAnnotationList))
    {
        # Calculate the position of annotation text
        mzIndexes <- mapUniqueMZ(mzAnnotationList[[i]][,1], 
                                 mzPeakList[[i]][,1], 
                                 tolerance = mzPeakTolerance,
                                 relativeTolerance = mzRelativePeakTolerance)
        
        # Exclude weak peaks
        mzIndexes <- mzIndexes[mzPeakList[[i]][mzIndexes, 2] >=
                               mzPlotAnnotationThreshold]
        # Exclude peaks of general isotopic patterns (2H, 13C, etc.)
        # Always annotate the first peak
        firstPeak <- mzIndexes == 1
        mzIndexes <- mzIndexes[!firstPeak]
        mzIndexes <- mzIndexes[mzPeakList[[i]][mzIndexes - 1, 1] + 1.05 < 
                               mzPeakList[[i]][mzIndexes, 1] | 
                               mzPeakList[[i]][mzIndexes - 1, 2]  <= 
                               mzPeakList[[i]][mzIndexes, 2]]
        mzIndexes <- c(rep(1, length(firstPeak[firstPeak])), mzIndexes)
        if (length(mzIndexes) == 0)
            next()
        
        mzIndexes <- mapUniqueMZ(mzPeakList[[i]][mzIndexes, 1], 
                                 mzDataList[[i]][,1],
                                 tolerance = mzPeakResolution, 
                                 relativeTolerance = FALSE)
        if (length(mzIndexes) > 0)
            plotAnnotationData <- 
                rbind(plotAnnotationData,
                      cbind(mzDataList[[i]][mzIndexes, 1],
                            mzPlotAnnotationPos(mzDataList[[i]][mzIndexes, 2]), 
                            i))
    }
}
if (mzPlotAnnotationUseX)
{
    plotAnnotationData <- cbind(plotAnnotationData, 
                                paste0(formatC(round(plotAnnotationData[,1], 3),
                                               digits = 3, format = 'f'),
                                       '\n', mzPlotAnnotationText))
} else 
{
    plotAnnotationData <- cbind(plotAnnotationData, 
                                mzPlotAnnotationText)
}
colnames(plotData) <- c('mz', 'intensity', 'id')
colnames(plotAnnotationData) <- c('x', 'y', 'id', 'label')

# Plot figures
p <- ggplot(plotData, aes(x = mz, y = intensity, group = id)) + 
    geom_line(aes(color = as.factor(id)), show.legend = FALSE) + 
    # scale_y_continuous(labels = function(X)
    #                    { return(format(X, scientific = TRUE)) }) +
    scale_x_continuous(breaks = mzPlotAxisTickX) + 
    scale_y_continuous(labels = function(X)
                       { formatC(X * 100, digits = 0, format = 'f') }) +
    scale_color_manual(values = mzPlotLineColor) + 
    facet_wrap(~id, nrow = mzPlotDimensionY, ncol = mzPlotDimensionX,
               scales = 'free_y', shrink = TRUE) + 
    coord_cartesian(xlim = mzPlotRangeX, ylim = mzPlotRangeY, expand = FALSE)

if (nrow(plotAnnotationData) > 0)
    p <- p + geom_text(plotAnnotationData, 
                       mapping = aes(x = x, y = y, group = id,
                                     label = label),
                       color = mzPlotAnnotationColor,
                       size = mzPlotFontSizeAnnotation,
                       lineheight = .75,
                       show.legend = FALSE)
p <- p + theme_bw() + 
    theme(text = element_text(size = mzPlotFontSize),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.title.y = element_text(margin = margin(r = 0.5, unit = 'line')),
          axis.title.x = element_text(margin = margin(t = 0.5, unit = 'line')),
          axis.text = element_text(size = mzPlotFontSizeAxisLabel),
          axis.ticks.length = unit(mzPlotFontSizeAxisLabel / 2, 'point'),
          panel.spacing.y = unit(2, 'line'),
          plot.margin = margin(t = 2, r = 4, b = 2, l = 2, unit = 'line')
         )
p <- p + 
    xlab(mzPlotLabelAxisX) + ylab(mzPlotLabelAxisY)

# Add title for subplots
plotTitlePosX <- mzPlotRangeX[1] * 0.05 + mzPlotRangeX[2] * 0.95
plotTitlePosY <- mzPlotRangeY[1] * 0.1 + mzPlotRangeY[2] * 0.9
plotText <- data.frame(x = plotTitlePosX, y = plotTitlePosY,
                       label = mzSpectrumName,
                       id = unique(plotData$id))
p <- p + geom_text(plotText, mapping = aes(x = x, y = y, label = label),
                   size = mzPlotFontSizeTitle)

# Save the plot
ggsave(mzPlotFilename, p, device = mzPlotDevice,
       width = mzPlotResolution[1], height = mzPlotResolution[2], 
       limitsize = FALSE)
