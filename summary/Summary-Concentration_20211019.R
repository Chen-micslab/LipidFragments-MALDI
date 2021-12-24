# R script for calculating properties (e.g. analyte concentration) based on 
# original data (e.g. peak intensities in spectra) and 
# established transforming functions (e.g. calibration curves)

source('lib/PathUtils.R')
source('lib/MSUtils.R')


mzStatDir <- '../data/summary/MALDI-TOF/'

mzDataFileName <- c(
    '25DHB+d3-Carnitine+SerumExtract-100%EtOH-20210908.csv',
    '25DHB+d3-Carnitine+SerumExtract-50%CHCl3-50%MeOH-20210908.csv'
)
mzDataFileName <- addParentDirectory(mzDataFileName, mzStatDir)

mzStatFilename <- 'Summary-Concentration+25DHB+d3-Carnitine+SerumExtract-20210908.csv'
mzStatFilename <- addParentDirectory(mzStatFilename, mzStatDir)

mzSampleName <- sapply(mzDataFileName, 
                       function(filename)
                       {
                           pos <- which(strsplit(filename, '')[[1]] == '/')
                           if (length(pos) > 0)
                               filename <- substr(filename, 
                                                  pos[length(pos)] + 1, 
                                                  nchar(filename))
                           pos <- which(strsplit(filename, '')[[1]] == '.')
                           if (length(pos) > 0)
                               filename <- substr(filename, 
                                                  1, pos[length(pos)] - 1)
                           return(filename)
                       })

mzRowID <- c(
    rep(184.0944, 2) # Carnitine, [C7H15NO3 + Na]+
)
mzRowIDTolerance <- 2e-4

mzColumnindex <- c(
    8, 9  # Int(m/z=184)
)
mzColumnName <- c(
    'Conc(Carnitine,uM)', 'SD(Carnitine,uM)'
)

mzTransformFunction <- c(
    # # Calibration of Carnitine
    function(X) { (X - 0) / 0.056 * 4 },
    # Std. deviation of Choline
    function(X) { X / 0.056 * 4 }
)


### MAIN ENTRY ###
# Read values of specified variables (columns) in given data files
variableList <- list()
for (dataFile in mzDataFileName)
{
    data <- read.csv(dataFile, as.is = TRUE)
    
    variables <- c()
    for (i in seq(1, length(mzRowID)))
    {
        selectedData <- data[data[,1] == mzRowID[i], mzColumnindex[i]]
        if (length(selectedData) <= 0)
        {
            # No exact match; try to find a sample by value match
            sampleIndexes <- mapUniqueMZ(mzRowID[i], data[,1], 
                                         mzRowIDTolerance)
            if (length(sampleIndexes) > 0)
                selectedData <- data[sampleIndexes, mzColumnindex[i]]
            else
                selectedData <- NA
        }
        variables <- c(variables, selectedData)
    }
    variableList <- c(variableList, list(variables))
}

# Swap the dimension (row <==> column)
variableList <- lapply(seq(1, length(mzColumnindex)), 
                       function(i) { sapply(variableList, function(X) {X[i]}) })

# Ugly Hack: convert Coeff. of Variance (CV) to Standard Deviation
for (i in seq(1, as.integer(length(variableList) / 2)))
    variableList[[i * 2]] <- variableList[[i * 2]] * variableList[[i * 2 - 1]]

# Apply transform functions to variables
for (i in seq(1, length(variableList)))
{
    if (i <= length(mzTransformFunction) && 
        !is.null(mzTransformFunction[i]))
        variableList[[i]] <- sapply(variableList[[i]], mzTransformFunction[[i]])
}

# Write results
outputData <- cbind(mzSampleName, as.data.frame(variableList))
colnames(outputData) <- c('Sample', mzColumnName)
write.csv(outputData, mzStatFilename, row.names = FALSE)