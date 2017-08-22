# Load libraries ----
require(fitdistrplus, quietly = TRUE, warn.conflicts = FALSE)
require(tidyr, quietly = TRUE, warn.conflicts = FALSE)
require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
require(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
require(lme4, quietly = TRUE, warn.conflicts = FALSE)
require(AICcmodavg, quietly = TRUE, warn.conflicts = FALSE)
require(magrittr, quietly = TRUE, warn.conflicts = FALSE)
require(rootSolve, quietly = TRUE, warn.conflicts = FALSE)
require(numDeriv, quietly = TRUE, warn.conflicts = FALSE)

# Used mainly in render functions and figure scripts
require(RColorBrewer, quietly = TRUE, warn.conflicts = FALSE)   # For brewer.pal function
require(scales, quietly = TRUE, warn.conflicts = FALSE)         # For hue_pal function
require(plot3D, quietly = TRUE, warn.conflicts = FALSE)         # For mesh function
require(gtable, quietly = TRUE, warn.conflicts = FALSE)         # For gtable_add_cols
require(grid, quietly = TRUE, warn.conflicts = FALSE)           # For unit.pmax
require(gridExtra, quietly = TRUE, warn.conflicts = FALSE)      # For grid.arrange
require(latex2exp, quietly = TRUE, warn.conflicts = FALSE)      # For latex expressions in plots
require(doParallel, quietly = TRUE, warn.conflicts = FALSE)     # For bootstrapping in Figures 5 and 6

if (require(rprojroot, quietly = TRUE, warn.conflicts = FALSE)) {
  mwROOT <- is_rstudio_project$make_fix_file()
  source(mwROOT("R","mw_mod.R"))
  source(mwROOT("R","mw_ipm.R"))
  cat("\nWelcome at", date(), "\n")
}
