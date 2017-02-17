# Herbivory decreases population growth of common milkweed (Asclepias syriaca) through negative effects on clonal reproduction

This repository contains the code for the above publication as an RStudio Project.

## Installation

 1. Download and install [R](https://www.r-project.org/) and [RStudio](https://www.rstudio.com/products/rstudio/download/).
 2. Download the code as a [ZIP archive](https://github.com/mdlama/milkweed/archive/master.zip) or in a terminal with the command `git clone https://github.com/mdlama/milkweed ./milkweed` 
 3. Navigate to the resulting directory on your computer.  *If you downloaded from the ZIP archive, you might have to first unzip the `milkweed-master.zip` file if it doesn't happen automatically.*
 4. Open the file `pubdat.Rproj`.  RStudio should open the project and begin installing all necessary R packages.  You should see the following
 
 ```
 Packrat is not installed in the local library -- attempting to bootstrap an installation...
 ```
 
 Note that these packages are being installed locally to this project and will *not* be available in R or RStudio outside of this project.
 5. Please wait while the R packages are installed.  It will be complete after you see `Packrat bootstrap successfully completed. Restarting R and entering packrat mode...`.  R will then restart, and a few libraries will be loaded.
 6. After Step 5 is complete, just this once, type `source(".Rprofile")` in the R console.  This will make sure all necessary libraries are loaded.  *Note that this command needs to be run only after the very first time this project is opened.*
 7. You should now be able to run all commands and scripts located in the `scripts` directory.  
 
 ## How to create the IPM
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
