# R script for plotting 1-D mass spectra for publication

library('ggplot2')
source('lib/PathUtils.R')
source('lib/MSFileIO.R')
source('lib/MSUtils.R')


mzDataDir <- '../data/spectra/MALDI-TOF'
mzStatDir <- '../data/summary/MALDI-TOF'
mzPlotDir <- '../data/plot'

mzSpectrumFilename <- c(
    '2022-02-11/TiO2NP25+NaCl+PC(18(1)-18(1))/0_N3/1/Spectrum_processed.rds',
    '2021-11-07/TiO2NP25+NaCl+PE(18-18)/0_E13/1/Spectrum_processed.rds',
    '2021-11-09/TiO2NP25+NaCl+PG(18(0)-18(1))/0_E20/1/Spectrum_processed.rds'
)
mzSpectrumFilename <- addParentDirectory(mzSpectrumFilename, mzDataDir)

mzPeakFileName <- changeFilenameSuffix(mzSpectrumFilename,
                                       '-peaks.csv', '.rds')

mzAnnotationPeakFilename <- c(
    'TiO2NP25+NaCl+PC(18(9Z)-18(9Z)).csv',
    'TiO2NP25+NaCl+PE(18-18)_noBG.csv',
    'TiO2NP25+NaCl+PG(18:0-18:1)_noBG.csv'
)
mzAnnotationPeakFilename <- addParentDirectory(mzAnnotationPeakFilename, mzStatDir)

mzPlotFilename <- 'TiO2NP25+PC+PE+PG.png'
# mzPlotFilename <- 'TiO2NP25+PC_mz86.png'
# mzPlotFilename <- 'TiO2NP25+PC_mz184.png'
# mzPlotFilename <- 'TiO2NP25+PC_mz805mz815.png'
mzPlotFilename <- addParentDirectory(mzPlotFilename, mzPlotDir)

mzSpectrumName <- c(
    'A', 'B', 'C'
    # ''
)

mzAnnotationPeakList <- list(
    # PC(18:1/18:1)
    c(55.0542, 
      # 72.0808, # Too crowded; annotated manually
      86.0964, 104.1070, 
      # 124.9998, # Too crowded; annotated manually
      146.9818, 166.0628, 184.0733, 198.0890, 504.3449,
      808.5827),
    # PE(18:1/18:1)
    c(55.0542, 
      # 69.0699, # Too crowded; annotated manually
      83.0855, 97.1012, 
      # 109.1012, # Too crowded; annotated manually
      164.0083, 
      # 185.9903, # Too crowded; annotated manually
      339.2894, 463.2196, 603.5347, 
      766.5357, 788.5177),
    # PG(18:0/18:1)
    c(57.0699, 71.0855, 
      # 83.0855, # Too crowded; annotated manually
      97.1012, 
      # 109.1012, # Too crowded; annotated manually
      173.0210, 195.0029, 216.9848, 341.3050, 605.5503,
      799.5460, 821.5279)
)

mzRange <- c(50, 1000)

mzPeakResolution <- 2e-5  # Average resolution of the detector (at m/z = 500)
mzPeakTolerance <- 0.05
mzRelativePeakTolerance <- FALSE

mzPlotDimensionX <- 1
mzPlotDimensionY <- 3

mzPlotLabelAxisX <- expression(italic('m') * '/' * italic('z'))
mzPlotLabelAxisY <- 'Relative Intensity (%)'

mzPlotRangeX <- c(40, 840)
# mzPlotRangeX <- c(85.5, 92.5)
# mzPlotRangeX <- c(183.5, 190.5)
# mzPlotRangeX <- c(805, 815)
mzPlotRangeY <- c(0, 1.12)
# mzPlotRangeY <- c(0, 0.02)
# mzPlotRangeY <- c(0, 0.02)
# mzPlotRangeY <- c(0, 0.003)

mzPlotAxisTickX <- seq(100, 800, 100)
# mzPlotAxisTickX <- seq(86, 92, 3)
# mzPlotAxisTickX <- seq(184, 190, 3)
# mzPlotAxisTickX <- seq(805, 815, 5)
mzPlotAxisTickY <- ggplot2::waiver()
# mzPlotAxisTickY <- seq(0, 0.02, 0.01)
# mzPlotAxisTickY <- seq(0, 0.003, 0.001)

mzPlotLineSize <- 1.5
mzPlotFontSize <- 100
mzPlotFontSizeTitle <- 40
mzPlotFontSizeAxisLabel <- 60
mzPlotFontSizeAnnotation <- 14
# mzPlotFontSizeAnnotation <- 20

mzPlotLineColor <- rep('#5F0000', 3)
mzPlotAnnotationColor <- '#000000'

mzPlotAnnotationText <- '▽'
# mzPlotAnnotationText <- ''
mzPlotAnnotationUseX <- TRUE
# mzPlotAnnotationUseX <- FALSE
mzPlotAnnotationThreshold <- (max(mzPlotRangeY) - min(mzPlotRangeY)) * 0.005
mzPlotAnnotationPos <- function(X)
{
    offset <- (mzPlotRangeY[2] - mzPlotRangeY[1]) * 0.06
    max <- mzPlotRangeY[2] - offset / 2
    pos <- X + offset
    pos[pos > max] <- max
    return(pos)
}

mzPlotDevice <- png
mzPlotResolution <- c(4000, 2400)
# mzPlotResolution <- c(800, 800)
# mzPlotResolution <- c(640, 530)
# mzPlotResolution <- c(900, 600)


### MAIN ENTRY ###
# Collect the spectrum data
mzDataList <- lapply(mzSpectrumFilename, function(X){ readSpectrum(X)} )
mzPeakList <- lapply(mzPeakFileName, 
                     function(X){ read.csv(X, as.is = TRUE)[,c(1,4)] })
mzAnnotationList <- lapply(mzAnnotationPeakFilename, 
                           function(X){ read.csv(X, as.is = TRUE) })

# Filter the m/z values
mzDataList <- lapply(mzDataList,
                     min = mzRange[1], max = mzRange[2],
                     function(X, min, max){ X[X[,1] >= min & X[,1] <= max,] })
mzPeakList <- lapply(mzPeakList,
                     min = mzRange[1], max = mzRange[2],
                     function(X, min, max){ X[X[,1] >= min & X[,1] <= max,] })
mzAnnotationList <- 
    lapply(seq(1, length(mzAnnotationList)),
           function(i)
           {
               return(mzAnnotationList[[i]]
                          [mapUniqueMZ(mzAnnotationPeakList[[i]],
                                       mzAnnotationList[[i]][,1],
                                       mzPeakTolerance,
                                       mzRelativePeakTolerance),])
           })

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
plotAnnotationData <- data.frame()
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
        mzIndexes <- mzIndexes[mzPeakList[[i]][mzIndexes - 1, 1] < 
                               mzPeakList[[i]][mzIndexes, 1] | 
                               mzPeakList[[i]][mzIndexes - 1, 1] + 1.05 < 
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
        plotAnnotationData <- rbind(plotAnnotationData,
                                    cbind(mzDataList[[i]][mzIndexes, 1],
                                          mzPlotAnnotationPos(
                                              mzDataList[[i]][mzIndexes, 2]), 
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
    geom_line(aes(color = as.factor(id)), size = mzPlotLineSize, 
              show.legend = FALSE) + 
    # scale_y_continuous(labels = function(X)
    #                    { return(format(X, scientific = TRUE)) }) +
    scale_x_continuous(breaks = mzPlotAxisTickX) + 
    scale_y_continuous(breaks = mzPlotAxisTickY,
                       labels = function(X)
                       { formatC(X * 100, digits = 1, format = 'f') }) +
    scale_color_manual(values = mzPlotLineColor) + 
    facet_wrap(~id, nrow = mzPlotDimensionY, ncol = mzPlotDimensionX,
               scales = 'free_y', shrink = TRUE) + 
    coord_cartesian(xlim = mzPlotRangeX, ylim = mzPlotRangeY, expand = FALSE)
p <- p + geom_text(plotAnnotationData, 
                   mapping = aes(x = x, y = y, group = id, label = label),
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
plotTitlePosX <- min(mzPlotRangeX) * 0.05 + max(mzPlotRangeX) * 0.98
plotTitlePosY <- min(mzPlotRangeY) * 0.1 + max(mzPlotRangeY) * 0.9
plotText <- data.frame(x = plotTitlePosX, y = plotTitlePosY,
                       label = mzSpectrumName,
                       id = unique(plotData$id))
p <- p + geom_text(plotText, mapping = aes(x = x, y = y, label = label),
                   size = mzPlotFontSizeTitle)

# Save the plot
ggsave(mzPlotFilename, p, device = mzPlotDevice,
       width = mzPlotResolution[1], height = mzPlotResolution[2], 
       limitsize = FALSE)
