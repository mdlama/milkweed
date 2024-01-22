################# Assign Genet IDs ############
##### adapted from clonal_sampling_simulation 4/7/23
#### Hannah Machiorlete 10/24/23

library(ggplot2)
library(tidyr) 
library(dplyr)

# import population demographic data
  setwd("C:/Users/Rebecca/Downloads") # working directory may be different for each computer
  demo_data <- read.csv("milkweed_demography_2023_full_wcoords_Oct12_HM.csv") # to import the whole table

# alternative way to get only the columns x, y, id
  #library("data.table")
  #demo_data <- fread("milkweed_demography_2023_full_wcoords_Oct12_HM.csv", select = c("plantID", "x", "y"))


  
  # data frame with plant ID and locations
  base_pop <- data.frame(site=demo_data$site_full, x=demo_data$x, y=demo_data$y, id=demo_data$plantID, umbels=demo_data$N_umbels)
  # remove duplicates
  base_pop <- distinct(base_pop, id, .keep_all = TRUE) # cut down to 1219 unique plant IDs (1234 entries total)
  # remove NA
  base_pop <- base_pop[!is.na(base_pop$site) & !is.na(base_pop$x) & !is.na(base_pop$y) & !is.na(base_pop$id) & !is.na(base_pop$umbels) ,] # down to 1207 plants with no NA entries

  blandy <- base_pop[base_pop$site=="blandy meadow",]           # should be rows 1-582    (563 plants)
  turner <- base_pop[base_pop$site=="turner pond",]             # should be rows 583-782  (197 plants)
  research <- base_pop[base_pop$site=="research village",]      # should be rows 783-985  (202 plants)
  sky <- base_pop[base_pop$site=="sky meadows",]                # should be rows 986-1235 (245 plants)
  
  popSize_blandy <- nrow(blandy)
  popSize_turner <- nrow(turner)
  popSize_research <- nrow(research)
  popSize_sky <- nrow(sky)
  
  # lower number of genets --> higher clonality
  numGenets_blandy <- .25 * popSize_blandy
  numGenets_turner <- .25 * popSize_turner
  numGenets_research <- .25 * popSize_research
  numGenets_sky <- .25 * popSize_sky
  
  
  ### consider these function parameters -- change these, the csv file names, and the graph labels 
  numGenets <- numGenets_blandy
  popSize <- popSize_blandy
  pop <- blandy
  
  ###### Genetic ID based on Performance ######
  trait_pop <- pop
  trait_pop$genetID <- NA
  # assign ramets with h_apical is NA genetic id 0
  trait_pop <- trait_pop %>% 
    mutate(genetID = ifelse(is.na(umbels),0,genetID))
  # assign ramets with h_apical < median genetic id 1
  trait_pop <- trait_pop %>% 
    mutate(genetID = ifelse(umbels < median(trait_pop$umbels),1,genetID))
  # assign ramets with h_apical > median genetic id 2
  trait_pop <- trait_pop %>% 
    mutate(genetID = ifelse(umbels >= median(trait_pop$umbels),2,genetID))
  
  ################ spatial pop #################### 
  ################ distance matrix ################ 
  xy_pop <- data.frame(x=pop$x, y=pop$y)                    # x, y data frame for calculating distance matrices
  
  d_pop <- as.matrix(dist(xy_pop, method = "euclidean"))    # distance matrix
  d_pop <- as.data.frame(d_pop)                             # convert to df
  colnames(d_pop) <- pop$id                                 # change column names to plant IDs
  d_pop$genetID <- pop$id                                   # add ID column
  dl_pop <- gather(d_pop, "id", distance, 1:popSize)        # long form
  
  ################ choose moms ###################
  moms <- sample(pop$id, numGenets, replace = F)          # choose moms
  pop$clone_m <- pop$id %in% moms
  orphanage <- dl_pop[dl_pop$id %in% moms == FALSE &      # plants except moms
                                dl_pop$distance != 0 &             # no self comparisons
                               !is.na(dl_pop$distance) &           # need a distance to assign genet
                                dl_pop$genetID %in% moms == TRUE,] # in relation to just the moms
  
  ################ create three pops ###################
  ###### phalanx ###### 
  phal_pop <- pop
  phal_pop$genetID <- NA
  phal_pop$distance <- NA
  
  phal_pop$id <- as.character(phal_pop$id)
  orphanage$id <- as.character(orphanage$id)
  
  ## fill in genetID and distance for moms (genetID = id, distance = 0)
  phal_pop <- phal_pop %>% 
    mutate(genetID = ifelse(clone_m == TRUE,id,genetID))
  phal_pop <- phal_pop %>% 
    mutate(distance = ifelse(clone_m == TRUE,0,distance))
  
  ## make data the same 
  orphanage$genetID <- as.numeric(orphanage$genetID)
  phal_pop$genetID <- as.numeric(phal_pop$genetID)
  
  ## choose plants until all plants have been assigned a genet
  while (sum(is.na(phal_pop$distance)) > 0) {
    newMoms <- orphanage %>% # choose miminum distance ramet
      group_by(genetID) %>% slice_min(n=1, distance) 
    newMoms <- newMoms %>% # only choose one ramet 
      group_by(id) %>% slice_min(n=1, distance) 
    phal_pop <- rows_update(phal_pop, select(newMoms, c("genetID", "id", "distance")), by = "id") # update the rows with the new genetIDs
    
    orphanage <- orphanage[!(orphanage$id %in% newMoms$id),] # remove the ramets that now have genetIDs 
  }
  phal2_pop <- phal_pop

  hist(table(as.character(phal_pop$genetID)), main="Frequency of Genet Sizes (# Ramets) at Blandy", xlab="Number of Ramets in Genet")

  ##### vary phalanx genet size ######
  test <- as.data.frame(table(phal2_pop$genetID))
  biggerGenets <- round(numGenets * 0.2)
  
  make_bigger <- sample(test$Var1, biggerGenets, replace = F)          # randomly choose some genets to make bigger
  for(w in 1:length(make_bigger)) {
    ramets <- phal2_pop$id[phal2_pop$genetID != make_bigger[w] & 
                                phal2_pop$clone_m == FALSE]                # choose all non-mom, non-genet ramets
    orphanage2 <- dl_pop[dl_pop$id %in% ramets == TRUE &                   # only ramets not in the genet already
                              dl_pop$genetID %in% make_bigger[w] == TRUE,] # genets to grow only
    new_ramets <- orphanage2 %>%                                           # choose 1-5 miminum distance ramets
      group_by(genetID) %>% slice_min(n= sample(1:5, 1), distance)
    phal2_pop <- rows_update(phal2_pop, select(new_ramets, c("id", "genetID", "distance")), by = "id") # update the rows with the new genetIDs
    
  }
  ### we could change distribution we sample from (e.g. poisson with different mean)
  ## could add rpois(1,4) or sample(1:5, 1)
  ## a few criteria: genets can't go extinct -- no self comparisons
  
  hist(table(phal_pop$genetID), main="Frequency of Genet Sizes (# Ramets) at Blandy", xlab="Number of Ramets in Genet")
  hist(table(phal2_pop$genetID), main="Frequency of (Varied) Genet Sizes (# Ramets) at Blandy", xlab="Number of Ramets in Genet")
  
  p1 <- as.data.frame(table(phal_pop$genetID))
  p2 <- as.data.frame(table(phal2_pop$genetID))
  mean(table(phal_pop$genetID))
  mean(table(phal2_pop$genetID))
  mean(p1$Freq)
  mean(p2$Freq)
  var(p1$Freq)
  var(p2$Freq)

  phal_pop <- phal_pop[order(phal_pop$genetID),]
  phal2_pop <- phal2_pop[order(phal2_pop$genetID),]
  
  write.csv(phal_pop, file="plalanx_population_blandy.csv")
  write.csv(phal2_pop, file="plalanx_population_VariedGenetSizes_blandy.csv")
  
  phal_pop$genetID <- as.factor(phal_pop$genetID)
  phalPlot =  
    ggplot(phal_pop, aes(x = x, y = y, color = genetID))+
    ggtitle("Phalanx Population")+
    geom_point(alpha = .7, size = 3)+
    theme_minimal()+
    theme(legend.position = "topleft")
    #scale_color_manual(values=c("#9B1200","#DA3907", "#FB8222", "#FDAF35", 
                                # "#F2CA3A", "#C0F435", "#7CD250", "#29AF7F", 
                                # "#1E9B8A", "#2F6F8E", "#3C508B", "#433E85", 
                                # "#742081", "#872781", "#C13A75", "#F4685C"))
  phalPlot
  
  phal2_pop$genetID <- as.factor(phal2_pop$genetID)
  phalPlot2 =  
    ggplot(phal2_pop, aes(x = x, y = y, color = genetID))+
    ggtitle("Phalanx Population, with Varied Genet Sizes")+
    geom_point(alpha = .7, size = 3)+
    theme_minimal()+
    theme(legend.position = "topleft")
    # scale_colour_gradient2(low="#210A10",mid="#DD4762", high="#FFFE9F")
  #   scale_color_manual(values=c("#9B1200","#DA3907", "#FB8222", "#FDAF35", 
  #                               "#F2CA3A", "#C0F435", "#7CD250", "#29AF7F", 
  #                               "#1E9B8A", "#2F6F8E", "#3C508B", "#433E85", 
  #                               "#742081", "#872781", "#C13A75", "#F4685C"))
  phalPlot2
  
  
  ###### intermediate ######
  #set.seed(199)
  #set.seed(runs[i])
  newclones <- sample(phal_pop$id[phal_pop$clone_m == FALSE], 
                      (popSize - numGenets) / 2, replace = F) 
  
  orphanage_int <- dl_pop[dl_pop$id %in% newclones == TRUE & # new genet plants
                               dl_pop$distance != 0 & # no self comparisons
                               dl_pop$genetID %in% moms == TRUE,] # moms only 
  
  # randomly assign genet ID to newclones 
  #set.seed(runs[i])
  #set.seed(199)
  new_moms_int <- orphanage_int %>%
    group_by(id) %>%
    slice_sample(n=1)
  
  # update rows, make consistent variable type
  new_moms_int$genetID <- as.character(new_moms_int$genetID)
  phal_pop$genetID <- as.character(phal_pop$genetID)
  
  int_pop <- rows_update(phal_pop, 
                         select(new_moms_int, c("id", "genetID", "distance")), 
                         by = "id")
  
  int_pop$genetID <- as.factor(int_pop$genetID)
  int_pop_plot =  
    ggplot(int_pop, aes(x = x, y = y, color = genetID))+
    ggtitle("Intermediate Population")+
    geom_point(alpha = .7, size = 3)+
    theme_minimal()+
    theme(legend.position = "topleft")
    #scale_colour_gradient2(low="#210A10",mid="#DD4762", high="#FFFE9F")
    # scale_color_manual(values=c("#9B1200","#DA3907", "#FB8222", "#FDAF35", 
    #                             "#F2CA3A", "#C0F435", "#7CD250", "#29AF7F", 
    #                             "#1E9B8A", "#2F6F8E", "#3C508B", "#433E85", 
    #                             "#742081", "#872781", "#C13A75", "#F4685C"))
  int_pop_plot
  write.csv(int_pop, file="intermediate_population_blandy.csv")
  
  
  ###### guerrilla ######
  orphanage_guer <- dl_pop[dl_pop$id %in% moms == FALSE & # no moms
                                dl_pop$distance != 0 & # no self comparisons
                                dl_pop$genetID %in% moms == TRUE,] # moms only 
  
  # randomly choose genets
  #set.seed(runs[i])
  #set.seed(199)
  new_moms_guer <- orphanage_guer %>%
    group_by(id) %>%
    slice_sample(n=1)
  
  # update rows, make consistent variable type
  new_moms_guer$genetID <- as.character(new_moms_guer$genetID)
  guer_pop <- rows_update(phal_pop, 
                          select(new_moms_guer, c("id", "genetID", "distance")), 
                          by = "id")
  
  guer_pop$genetID <- as.factor(guer_pop$genetID)
  guer_pop_plot =  
    ggplot(guer_pop, aes(x = x, y = y, color = genetID))+
    ggtitle("Guerrilla Population")+
    geom_point(alpha = .7, size = 3)+
    theme_minimal()+
    theme(legend.position = "topleft")
    #scale_colour_gradient2(low="#210A10",mid="#DD4762", high="#FFFE9F")
    # scale_color_manual(values=c("#9B1200","#DA3907", "#FB8222", "#FDAF35", 
    #                             "#F2CA3A", "#C0F435", "#7CD250", "#29AF7F", 
    #                             "#1E9B8A", "#2F6F8E", "#3C508B", "#433E85", 
    #                             "#742081", "#872781", "#C13A75", "#F4685C"))
  guer_pop_plot
  write.csv(guer_pop, file="guerrilla_population_blandy.csv")