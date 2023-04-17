library(data.table)
library(ggplot2)
library(readxl)
require(tidyverse)
library(dplyr)

rm(list = ls())    # remove any variables in R's memory 

#-------------------------------------------------------------------------------
# open file
#-------------------------------------------------------------------------------
hathibaseline <- fread("C:/Users/Hyolim/OneDrive - London School of Hygiene and Tropical Medicine/SNU/HATHI/hathi_baseline.csv")
hathibaseline <- fread("D:/OneDrive - London School of Hygiene and Tropical Medicine/SNU/HATHI/hathi_baseline.csv")

cols_to_convert <- 7:163

for (colname in names(hathibaseline)[cols_to_convert]) {
  hathibaseline[, (colname) := as.integer(.SD[[colname]] == "Yes")]
}

fwrite(hathibaseline, "D:/OneDrive - London School of Hygiene and Tropical Medicine/SNU/HATHI/hathibinary.csv")

col_to_plot <- c(7:163)

colnames(hathibaseline) <- paste(1:163)

hathibaseline[is.na(hathibaseline)] <- 0

hathibaseline <- data.frame(hathibaseline)

quantile <- apply(hathibaseline[,87:92], 2, mean)
meandf <- data.frame(variable = names(quantile), value = quantile)


