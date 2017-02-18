This repository contains the code for the above publication as an RStudio Project.

# Installation

 1. Download and install [R](https://www.r-project.org/) and [RStudio](https://www.rstudio.com/products/rstudio/download/).
 1. Download the code as a [ZIP archive](https://github.com/mdlama/milkweed/archive/master.zip) or in a terminal with the command `git clone https://github.com/mdlama/milkweed ./milkweed` 
 1. Navigate to the resulting directory on your computer.  *If you downloaded from the ZIP archive, you might have to first unzip the `milkweed-master.zip` file if it doesn't happen automatically.*
 1. Open the file `pubdat.Rproj`.  RStudio should open the project and begin installing all necessary R packages.  You should see the following
 
 ```
 Packrat is not installed in the local library -- attempting to bootstrap an installation...
 ```
 
 Note that these packages are being installed locally to this project and will *not* be available in R or RStudio outside of this project.
 1. Please wait while the R packages are installed.  It will be complete after you see `Packrat bootstrap successfully completed. Restarting R and entering packrat mode...`.  R will then restart, and a few libraries will be loaded.
 1. After the previous step is complete, just this once, type `source(".Rprofile")` in the R console.  This will make sure all necessary libraries are loaded.  *Note that this command needs to be run only after the very first time this project is opened.*
 1. You should now be able to run all commands and scripts located in the `scripts/figs` and `scripts/model_selection` directories.
 
# How to create the IPM
You can create the IPM with the command
 
```
ipm <- mwIPM()
```
 
The IPM will be built, and resulting computations saved in the `data/calculated` directory to avoid unnecessary computations in the future.  If at any time you would like to recompute the IPM, you can pass `compute=TRUE` as an argument in a list:
 
```
ipm <- mwIPM(list(compute = TRUE))
```
 
This will by default not overwrite any previously saved computations in the `data/calculated` directory.  To save the results of the computations, you can pass `saveresults=TRUE` as an additional argument in a list:
 
```
ipm <- mwIPM(list(compute = TRUE, saveresults = TRUE))
```

# Figure scripts

[Figure 1](scripts/figs/Figure1/Figure1.html) <br \>
[Figure 2](scripts/figs/Figure2/Figure2.html) <br \>
[Figure 3](scripts/figs/Figure3/Figure3.html) <br \>
[Figure 4](scripts/figs/Figure4/Figure4.html) <br \>
[Figure 5](scripts/figs/Figure5/Figure5.html) <br \>
[Appendix: Figure 2](scripts/figs/AppendixFigure2/AppendixFigure2.html) <br \>
[Appendix: Figure 3](scripts/figs/AppendixFigure3/AppendixFigure3.html) <br \>
[Appendix: Figure 4](scripts/figs/AppendixFigure4/AppendixFigure4.html)
