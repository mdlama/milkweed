#### -- Packrat Autoloader (version 0.4.8-1) -- ####
source("packrat/init.R")
#### -- End Packrat Autoloader -- ####

# Load libraries ----
require(fitdistrplus)
require(tidyr)
require(dplyr)
require(ggplot2)
require(lme4)
require(AICcmodavg)
require(magrittr)
require(rootSolve)
require(numDeriv)

# Used mainly in render functions and figure scripts
require(RColorBrewer)   # For brewer.pal function
require(scales)         # For hue_pal function
require(plot3D)         # For mesh function
require(gtable)         # For gtable_add_cols
require(grid)           # For unit.pmax
require(gridExtra)      # For grid.arrange
require(latex2exp)      # For latex expressions in plots
require(doParallel)     # For bootstrapping in Figures 4 and 5

if (require(rprojroot)) {
  mwROOT <- is_rstudio_project$make_fix_file()
}

# Weird need to put the sourcing in a function.  If not, only the first source command works.
# Also, have to actually call the .First function as it won't be called otherwise.  Maybe because
#   .First doesn't work as expected as a project .Rprofile?
.First <- function() {
  source(mwROOT("includes","mwMod.R"))
  source(mwROOT("includes","mwIPM.R"))
  cat("\nWelcome at", date(), "\n")
}
.First()