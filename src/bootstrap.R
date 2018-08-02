# Bootstraps the whole project for testing purposies
rm(list=ls())

# Return the memory to the OS
gc(verbose = FALSE,
   reset = TRUE)

setwd("/home/theodor/Documents/ProteoSign-2.0/")

source("src/main.R")