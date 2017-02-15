library(stats)  # To avoid non-masking of stats filter function.  Need to move to package workflow.
library(fitdistrplus)
library(rprojroot)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(AICcmodavg)
library(magrittr)
library(rootSolve)
library(numDeriv)

# Used mainly in render functions and figure scripts
library(RColorBrewer)   # For brewer.pal function
library(scales)         # For hue_pal function
library(plot3D)         # For mesh function
library(gtable)         # For gtable_add_cols
library(grid)           # For unit.pmax
library(gridExtra)      # For grid.arrange
library(latex2exp)      # For latex expressions in plots
library(doParallel)     # For bootstrapping in Figures 4 and 5

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