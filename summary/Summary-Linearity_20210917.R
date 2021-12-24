# R script for evaluating linearity of response of targeted metabolites
# Using functions provided by ggplot2 package

library('ggplot2')
source('lib/MSUtils.R')


mzStatDir <- '../data/summary/HPLC-MS'
mzPlotDir <- '../data/plot'

mzPeakFileName <- c(
    'SPCS+d9-Choline+HomoPhe.csv',
    'SPCS+d9-Choline+HomoPhe+Choline(5uM).csv',
    'SPCS+d9-Choline+HomoPhe+Choline(10uM).csv',
    'SPCS+d9-Choline+HomoPhe+Choline(20uM).csv',
    'SPCS+d9-Choline+HomoPhe+Choline(30uM).csv',
    'SPCS+d9-Choline+HomoPhe+Choline(40uM).csv'
)

mzPlotSummaryFileName <- c(
    'SPCS+Choline_RefPeakIntensity_mz104.png'
)

mzStatMZ <- c(
    104.1070 # Choline, [C5H14NO]+
)

mzBackgroundMZ <- c(
    # 80.9484 # NaCl, [NaCl + Na]+
)
mzBackgroundWeight <- c(
    # 0.31995776 # {37Cl} / {35Cl}
)

mzPeakTolerance <- 0.1
mzPeakRelativeTolerance <- FALSE

mzStatVariableList <- list(
    c(0, 5, 10, 20, 30, 40) # Choline, umol/L
)

# mzStatIntensityIndex <- 2 # Raw intensities
mzStatIntensityIndex <- 6 # Intensities normalized by m/z=113 (d9-Choline)
# mzStatIntensityCVIndex <- 3 # Raw Coeff. of Variance
mzStatIntensityCVIndex <- 7 # Coeff. of Variance of normalized intensities

mzPlotAxisNameX <- 'Added Concentration (mmol/L)'
mzPlotAxisNameY <- 'Normalized Intensity (m/z=104 vs 113)'
mzPlotTitle <- c('Response of Choline (m/z = 104)')
mzPlotSummaryResolution <- c(600, 600)
    
mzPeakFileName <- paste0(mzStatDir, '/', mzPeakFileName)
mzPlotSummaryFileName <- paste0(mzPlotDir, '/', mzPlotSummaryFileName)


### MAIN ENTRY ###

# Read peaks' intensities and their CVs from data files
for (i in seq(1, length(mzStatMZ)))
{
    intensities <- c()
    intensityCV <- c()
    
    for (mzPeakFile in mzPeakFileName)
    {
        mzPeaks <- read.csv(mzPeakFile, as.is = TRUE)
        
        # Select specified m/z and their intensities
        mzIndexes <- mapUniqueMZ(mzStatMZ[i], mzPeaks[,1], 
                                 tolerance = mzPeakTolerance,
                                 relativeTolerance = mzPeakRelativeTolerance,
                                 noMatchAsNA = TRUE)
        bgIndexes <- mapUniqueMZ(mzBackgroundMZ[i], mzPeaks[,1], 
                                 tolerance = mzPeakTolerance, 
                                 relativeTolerance = mzPeakRelativeTolerance,
                                 noMatchAsNA = TRUE)
        if (is.na(mzIndexes))
        {
            intensities <- c(intensities, 0)
            intensityCV <- c(intensityCV, NA)
        }
        else
        {
            intensity <- mzPeaks[mzIndexes, mzStatIntensityIndex]
            if (all(is.na(bgIndexes)))
                bgIntensity <- 0
            else
                bgIntensity <- mzPeaks[bgIndexes, mzStatIntensityIndex] *
                               mzBackgroundWeight[i]
            intensities <- c(intensities, 
                             intensity - bgIntensity)
            intensityCV <- c(intensityCV, 
                             mzPeaks[mzIndexes, mzStatIntensityCVIndex])
        }
    }
    
    # Linear regression
    statData <- data.frame(x = mzStatVariableList[[i]], y = intensities)
    statLinearity <- lm(y ~ x, statData)
    
    # Plot statistics
    if (length(mzPlotSummaryFileName) < i)
        next()
    
    p <- ggplot(as.data.frame(cbind(x = mzStatVariableList[[i]], 
                                    y = intensities)),
                aes(x = x, y = y)) + 
         theme_bw() + 
         theme(plot.title = element_text(size = 30, hjust = 0.5),
               axis.title.x = element_text(size = 15),
               axis.title.y = element_text(size = 15)) + 
         geom_point() +
         geom_errorbar(aes(x = x, 
                           ymin = intensities - intensities * intensityCV,
                           ymax = intensities + intensities * intensityCV),
                       width = (max(mzStatVariableList[[i]]) - 
                                min(mzStatVariableList[[i]])) / 20)
    # Plot the fitted curve
    p <- p +
        geom_smooth(method = 'lm', se = FALSE,
                    colour = 'blue', alpha = 0.6, linetype = 'dotted') +
        geom_text(aes(label = paste(
            'y = ', round(statLinearity$coefficients[2], 3),
            ' x + ', round(statLinearity$coefficients[1], 3),
            '\n',
            'R^2 =',
            round(summary(statLinearity)[['r.squared']], 3)),
            x = max(statData$x) * 0.3 + min(statData$x) * 0.7,
            y = max(intensities) * 0.9 + min(intensities) * 0.1),
            size = 10,
            show.legend = FALSE)
    
    p <- p + 
         xlab(mzPlotAxisNameX) + ylab(mzPlotAxisNameY) + ggtitle(mzPlotTitle[i])
    ggsave(filename = mzPlotSummaryFileName[i], plot = p, device = png,
           width = mzPlotSummaryResolution[1],
           height = mzPlotSummaryResolution[2], 
           dpi = 200, limitsize = FALSE)
}