# Script for data import of actin and TOCA fluorescence in identified structures over time

#-------------------------------------------------------------------------------
# DEPENDENCIES:

# If you do not have TidyR package installed, execute this line (once only):
# install.packages('tidyr', dependencies=TRUE, repos='http://cran.rstudio.com/')
library(tidyr)

#-------------------------------------------------------------------------------
# Set directory location, find data files:

getwd()
setwd("/Users/Lab/Documents/Postdoc/Hanae/CrosscorrelationActinToca/DummyData")
list.files()

#-------------------------------------------------------------------------------

# Example import table:
curr.table <- "Results_ACTIN.txt"
  
# import table, skipping title line:
mytable <- read.table(curr.table, skip = 1)

# read title line
topline = readLines(curr.table)[1]
topline.names = unlist(strsplit(topline, "\t"))

# some cleaning up with column names: 
ncols <- ncol(mytable)
nvars <- length(topline.names)
ncols  # 22
nvars  # 7
colnames(mytable)[(ncols-nvars+1):ncols] <- unlist(topline.names)

# Create ROI column?
# (notes - this is a quick and  dirty fix; for futureproofing, make a regex
# read with 'grep' through columns to find the ROI info there instead)

mytable["ROI"] <- rep(NA, nrow(mytable))

n.slices <- max(mytable["Slice"])
n.spots <- nrow(mytable)/n.slices

mytable["ROI"] <- rep(1:n.spots, n.slices)

mytable[, c(10, 19:ncol(mytable))]  
  ## ATTENTION! The names are shifted relative to Fiji ROI names wherever an ROI was deleted

actin <- mytable


# Repeat the above with actin.bg, TOCA, TOCA.bg

#-------------------------------------------------------------------------------
# rewrite above import as function:
rm(mytable)

ImportData <- function(filename) {
  curr.table <- filename
  mytable <- read.table(curr.table, skip = 1)
  topline = readLines(curr.table)[1]
  topline.names = unlist(strsplit(topline, "\t"))
  ncols <- ncol(mytable)
  nvars <- length(topline.names)
  colnames(mytable)[(ncols-nvars+1):ncols] <- unlist(topline.names)
  mytable["ROI"] <- rep(NA, nrow(mytable))
  n.slices <- max(mytable["Slice"])
  n.spots <- nrow(mytable)/n.slices
  mytable["ROI"] <- rep(1:n.spots, n.slices)
  return(mytable)
}
 
# check identity: 
identical(ImportData("Results_ACTIN.txt"), actin)
 
#------------------------------------------------------------------------------- 
# Import each data table:

list.files()  # shows files in home folder
actin <- ImportData("Results_ACTIN.txt")  
actin.bg <- ImportData("Results_ACTIN_rollingball.txt")  
toca <- ImportData("Results_TOCA.txt")  
toca.bg <- ImportData("Results_TOCA_rollingball.txt")  
  
# This is the 'tidy' data format...  
# for plugging it into my scripts, create time series per column

help(spread)
# minimal example
miniactin <- actin[1:58, c(17:23)]; miniactin  
spread(miniactin, ROI, Mean)  # not working as desired. Another dirty fix for now...

miniactin

getwd()

# AIM: Wrangle data into such a format so that: 
# - each column = 1 structure
# - each row = 1 timepoint
# - one such table for actin measurements, one for toca

mytable <- actin

wrangletable <- data.frame(matrix(NA, ncol = n.spots, nrow = n.slices))
colnames(wrangletable) = paste0("ROI", 1:n.spots)
rownames(wrangletable) = 1:n.slices
wrangletable 
  
for(i in 1:n.spots) {
  curr.rows <- which(mytable$ROI == i)  
  wrangletable[, i] <- mytable[curr.rows, ]
  wrangletable[, i] <- mytable[curr.rows, ]$Mean
}

return(wrangletable)

# Rewrite this as function to apply to all datasets:

Wrangle <- function(input.table) {
  mytable <- input.table
  wrangletable <- data.frame(matrix(NA, ncol = n.spots, nrow = n.slices))
  colnames(wrangletable) = paste0("ROI", 1:n.spots)
  rownames(wrangletable) = 1:n.slices
  for(i in 1:n.spots) {
    curr.rows <- which(mytable$ROI == i)  
    wrangletable[, i] <- mytable[curr.rows, ]
    wrangletable[, i] <- mytable[curr.rows, ]$Mean
  }
  return(wrangletable)
}
 
identical(Wrangle(actin), wrangletable) # Yes! Good

w.actin <- Wrangle(actin)
w.toca <- Wrangle(toca)
w.actin.bg <- Wrangle(actin.bg)
w.toca.bg <- Wrangle(toca.bg)


matplot(w.toca, type = "l", col = "blue", lty = 1,
        xlab = "Timepoint",
        ylab = "Fluorescence [au]",
        ylim = range(w.actin, w.toca))
matplot(w.actin, type = "l", col= "red", lty = 1, add = TRUE)

write.table(w.actin, "w_actin.txt")
write.table(w.actin.bg, "w_actin-bg.txt")
write.table(w.toca, "w_toca.txt")
write.table(w.toca.bg, "w_toca-bg.txt")







