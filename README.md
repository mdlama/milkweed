This repository contains the code for the above publication as an RStudio Project.

# Installation

 - Download and install [R](https://www.r-project.org/) and [RStudio](https://www.rstudio.com/products/rstudio/download/).
 - Download the code as a [ZIP archive](https://github.com/mdlama/milkweed/archive/master.zip) or in a terminal with the command `git clone https://github.com/mdlama/milkweed ./milkweed` 
 - Navigate to the resulting directory on your computer.  *If you downloaded from the ZIP archive, you might have to first unzip the `milkweed-master.zip` file if it doesn't happen automatically.*
 - Open the file `milkweed.Rproj`.  RStudio should open the project and begin installing all necessary R packages.  You should see the following
 
 ```
 Packrat is not installed in the local library -- attempting to bootstrap an installation...
 ```
 
 Note that these packages are being installed locally to this project and will *not* be available in R or RStudio outside of this project.
 - Please wait while the R packages are installed.  It will be complete after you see `Packrat bootstrap successfully completed. Restarting R and entering packrat mode...`.  R will then restart, and a few libraries will be loaded.
 - After the previous step is complete, just this once, type `source(".Rprofile")` in the R console.  This will make sure all necessary libraries are loaded.  *Note that this command needs to be run only after the very first time this project is opened.*
 - You should now be able to run all commands and scripts located in the `scripts/figs` and `scripts/model_selection` directories.
 
# How to create the IPM
You can create the IPM with the command
 
```r
ipm <- mwIPM()
```
 
The IPM will be built, and resulting computations saved in the `data/calculated` directory to avoid unnecessary computations in the future.  If at any time you would like to recompute the IPM, you can pass `compute=TRUE` as an argument in a list:
 
```r
ipm <- mwIPM(list(compute = TRUE))
```
 
This will by default not overwrite any previously saved computations in the `data/calculated` directory.  To save the results of the computations, you can pass `saveresults=TRUE` as an additional argument in a list:
 
```r
ipm <- mwIPM(list(compute = TRUE, saveresults = TRUE))
```

# Figure scripts

 - [Figure 2: Flowering and growth](scripts/figs/Figure2_FloweringAndGrowth/Figure2_FloweringAndGrowth.html)
 - [Figure 3: Herbivory slope parameters across sites](scripts/figs/Figure3_HerbSlopePars/Figure3_HerbSlopePars.html)
 - [Figure 4: Clonal growth versus herbivory](scripts/figs/Figure4_ClonalGrowthVsHerb/Figure4_ClonalGrowthVsHerb.html)
 - [Figure 5: Population growth versus herbivory](scripts/figs/Figure5_PopGrowthVsHerb/Figure5_PopGrowthVsHerb.html)
 - [Figure 6: Elasticity analysis](scripts/figs/Figure6_ElasticityAnalysis/Figure6_ElasticityAnalysis.html)
 - [Appendix Figure 2: Recruitment and herbivory distributions](scripts/figs/AppendixFigure2_RecruitAndHerbDist/AppendixFigure2_RecruitAndHerbDist.html)
 - [Appendix Figure 3: Herbivory versus apical height](scripts/figs/AppendixFigure3_HerbVsHeight/AppendixFigure3_HerbVsHeight.html)
 - [Appendix Figure 4: Sprouts per stem versus apical height](scripts/figs/AppendixFigure4_SproutsPerStemVsHeight/AppendixFigure4_SproutsPerStemVsHeight.html)

# Model selection scripts

 - [Herbivory Distributions](scripts/model_selection/FitHerbivoryDistribution.html)
 - [Sprout Distributions](scripts/model_selection/FitBudlingDistribution.html)
 - [Seedling Distribution](scripts/model_selection/FitSeedlingDistribution.html)
 - [Sprouts-per-stem Distribution](scripts/model_selection/FitBudlingsPerStem.html)
