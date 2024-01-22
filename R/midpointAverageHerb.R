# midpoint average for herbivory
setwd("C:/Users/Rebecca/Downloads") # working directory may be different for each computer
mdata <- read.csv("milkweed_demography_2023_full_wcoords_Oct12_HM.csv") # to import the whole table
# put data into data frame
mdf <- data.frame(herb0 = mdata$herb0, herb1 = mdata$herb1, herb2= mdata$herb2, herb3= mdata$herb3, herb4=mdata$herb4, herb5=mdata$herb5, herb6=mdata$herb6 )
# remove NA's
mdf <- mdf[!is.na(mdf$herb0),]
# weights is the midpoint of percent tissue removed for each section
# weights <- c(0, (0+5)/2, (5+25)/2, (25+50)/2, (50+75)/2, (75+99)/2, 100)
mdf$"avgHerb(midpoint average)" <- (mdf$herb0 * 0 +
                  mdf$herb1 * (0+.05)/2 +
                  mdf$herb2 * (.05+.25)/2 +
                  mdf$herb3 * (.25+.50)/2 +
                  mdf$herb4 * (.50+.75)/2 +
                  mdf$herb5 * (.75+.99)/2 +
                  mdf$herb6 * 1
                ) / ( mdf$herb0  + mdf$herb1 + mdf$herb2 +
                    mdf$herb3  + mdf$herb4  + mdf$herb5  + mdf$herb6) 
