# All functions exist in this file. If they are too many, they can be separeted in
# several functions_XXX.R files.

# Clear enviroment
rm(list = grep("^global.variables|^project.variables", ls(), value = TRUE, invert = TRUE))

# Functions scripts to load
functions.subfiles <- c("functions_build.R")#,
                        #"functions_analyze.R")


invisible(lapply(functions.subfiles, source))
