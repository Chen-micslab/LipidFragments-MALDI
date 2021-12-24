# R script for calculating properties (e.g. analyte concentration) based on 
# original data (e.g. peak intensities in spectra) and 
# established transforming functions (e.g. calibration curves)

source('lib/PathUtils.R')
source('lib/MSUtils.R')


mzStatDir <- '../data/summary/HPLC-MS/'

mzDataFileName <- c(
    'SerumExtract-100%EtOH-20210908+d3-Carnitine.csv',
    'SerumExtract-50%CHCl3-50%MeOH-20210908+d3-Carnitine.csv'
)
mzDataFileName <- addParentDirectory(mzDataFileName, mzStatDir)

mzStatFilename <- 'Summary-Concentration-SerumExtract-20210908+d3-Carnitine.csv'
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
    rep(162.1125, 2) # Carnitine, [C7H15NO3 + H]+
)
mzRowIDTolerance <- 1e-3

mzColumnindex <- c(
    6, 7  # Int(m/z=162)
)
mzColumnName <- c(
    'Conc(Carnitine,uM)', 'SD(Carnitine,uM)'
)

mzTransformFunction <- c(
    # # Calibration of Carnitine
    function(X) { (X + 0.071) / 0.06 * 4 },
    # Std. deviation of Carnitine
    function(X) { X / 0.06 * 4 }
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