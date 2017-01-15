library(stats)  # To avoid non-masking of stats filter function.  Need to move to package workflow.
library(fitdistrplus)
library(rprojroot)
library(dplyr)
library(ggplot2)
library(lme4)
library(magrittr)
library(numDeriv)

# Used mainly in render functions
library(RColorBrewer)   # For brewer.pal function
library(scales)         # For hue_pal function
library(plot3D)         # For mesh function
library(gtable)         # For gtable_add_cols
library(gridExtra)      # For grid.arrange

mwROOT <- is_rstudio_project$make_fix_file()

# Weird need to put the sourcing in a function.  If not, only the first source command works.
# Also, have to actually call the .First function as it won't be called otherwise.  Maybe because
#   .First doesn't work as expected as a project .Rprofile?
.First <- function() {
  source(mwROOT("includes","mwMod.R"))
  source(mwROOT("includes","mwIPM.R"))
  cat("\nWelcome at", date(), "\n")
}

.First()