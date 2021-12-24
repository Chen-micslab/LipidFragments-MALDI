# R script for dealing with file/directory path

insertFilenameSuffix <- function(filename, suffix)
{
    if (is.list(filename))
        return(lapply(filename, suffix = suffix, insertFilenameSuffix))
    
    pos <- which(strsplit(filename, '')[[1]] == '.')
    if (length(pos) > 0)
    {
        pos <- pos[length(pos)]
        return(paste0(substr(filename, 1, pos - 1),
                      suffix,
                      substr(filename, pos, nchar(filename))))
    } else
        return(paste0(filename, suffix))
}

changeFilenameSuffix <- function(filename, newSuffix, oldSuffix = '')
{
    if (is.list(filename))
        return(lapply(filename, 
                      newSuffix = newSuffix, 
                      oldSuffix = oldSuffix,
                      changeFilenameSuffix))
    
    return(sapply(filename, 
                  oldSuffix = oldSuffix, 
                  newSuffix = newSuffix,
                  function(X, oldSuffix, newSuffix)
                  {
                      if (nchar(X) == 0)
                          return('')
                      
                      if (nchar(oldSuffix) == 0)
                      {
                          pos <- grep('.', strsplit(X, '')[[1]][length(X):1], 
                                      fixed = TRUE)
                          if (length(pos) > 0)
                              return(paste0(substr(X, 1, pos[1] - 1), 
                                            newSuffix))
                          else
                              return(paste0(X, newSuffix))
                      }
                      else
                      {
                          pos <- regexpr(oldSuffix, X, fixed = TRUE)[1]
                          if (pos > 0)
                              return(paste0(substr(X, 1 , pos - 1), 
                                            newSuffix))
                          else
                              return(X)
                      }
                  }))
}

addParentDirectory <- function(filename, directory, skipEmpty = TRUE)
{
    if (is.list(filename))
        return(lapply(filename, 
                      directory = directory, 
                      skipEmpty = skipEmpty,
                      addParentDirectory))
    
    if (length(filename) == 0)
        return(NULL)
    
    lastChar <- substr(directory, nchar(directory), nchar(directory))
    if (lastChar != '/' && lastChar != '\\')
        directory <- paste0(directory, '/')
    return(sapply(filename, 
                  USE.NAMES = FALSE,
                  skipEmpty = skipEmpty,
                  function(X, skipEmpty)
                  {
                      if (is.na(X))
                          return(NA)
                      if (skipEmpty && nchar(X) == 0)
                          return('')
                      return(paste0(directory, X))
                  }))
}
