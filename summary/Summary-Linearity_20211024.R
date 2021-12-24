# R script for evaluating linearity of response of targeted metabolites
# Using functions provided by ggplot2 package

library('ggplot2')
source('lib/MSUtils.R')


mzStatDir <- '../data/summary/MALDI-TOF'
mzPlotDir <- '../data/plot'

mzPeakFileName <- c(
    '25DHB+d3-Carnitine+SPCS.csv',
    '25DHB+d3-Carnitine+Carnitine(5uM).csv',
    '25DHB+d3-Carnitine+Carnitine(10uM).csv',
    '25DHB+d3-Carnitine+Carnitine(20uM).csv',
    '25DHB+d3-Carnitine+Carnitine(30uM).csv'
)

mzPlotSummaryFileName <- c(
    # '25DHB+SPCS+Carnitine_RefPeakIntensity_mz162.png'
    '25DHB+SPCS+Carnitine_RefPeakIntensity_mz184.png'
)

mzStatMZ <- c(
    # 162.1125 # Carnitine, [C7H15NO3 + H]+
    184.0944 # Carnitine, [C7H15NO3 + Na]+
)

mzBackgroundMZ <- c(
    # 80.9484 # NaCl, [NaCl + Na]+
)
mzBackgroundWeight <- c(
    # 0.31995776 # {37Cl} / {35Cl}
)

mzPeakTolerance <- 1e-4

mzStatVariableList <- list(
    c(0, 5, 10, 20, 30) # Carnitine, umol/L
)

# mzStatIntensityIndex <- 2 # Raw intensities
# mzStatIntensityIndex <- 6 # Intensities normalized by m/z=165 (d3-Carnitine)
mzStatIntensityIndex <- 8 # Intensities normalized by m/z=187 (d3-Carnitine)
# mzStatIntensityCVIndex <- 3 # Raw Coeff. of Variance
# mzStatIntensityCVIndex <- 7 # Coeff. of Variance of normalized intensities
mzStatIntensityCVIndex <- 9 # Coeff. of Variance of normalized intensities

mzPlotAxisNameX <- 'Added Concentration (mmol/L)'
# mzPlotAxisNameY <- 'Normalized Intensity (m/z=162 vs 165)'
mzPlotAxisNameY <- 'Normalized Intensity (m/z=184 vs 187)'
# mzPlotTitle <- c('Response of Carnitine (m/z = 162)')
mzPlotTitle <- c('Response of Carnitine (m/z = 184)')
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
                                 noMatchAsNA = TRUE)
        bgIndexes <- mapUniqueMZ(mzBackgroundMZ[i], mzPeaks[,1], 
                                 tolerance = mzPeakTolerance, 
                                 noMatchAsNA = TRUE)
        if (is.na(mzIndexes))
        {
            intensities <- c(intensities, 0)
            intensityCV <- c(intensityCV, 0)
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
                             formatC(round(
                                     summary(statLinearity)[['r.squared']], 3),
                                     digits = 3, format = 'f')),
                       x = quantile(statData$x)[2],
                       y = max(intensities)),
                   size = 10,
                   show.legend = FALSE)
    
    p <- p + 
         xlab(mzPlotAxisNameX) + ylab(mzPlotAxisNameY) + ggtitle(mzPlotTitle[i])
    ggsave(filename = mzPlotSummaryFileName[i], plot = p, device = png,
           width = mzPlotSummaryResolution[1],
           height = mzPlotSummaryResolution[2], 
           dpi = 200, limitsize = FALSE)
}