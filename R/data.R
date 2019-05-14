#' Stem Data!
#'
#' A dataset containing demographic information for milkweed plants from 7 sites across VA and PA from 2013-2017.
#'
#' @format A data frame with 5829 rows and 16 variables:
#' \describe{
#'   \item{year}{collection year}
#'   \item{site}{3-digit site code: BLD = Blandy Experimental Farm (actually two sites, 1 indicating the area that is unburned, 2 the area that is burned)
#'                                  SKY = Sky Meadows State Park
#'                                  PWR = Presquile National Wildlife Refuge
#'                                  YTB = Yorktown National Battlefield
#'                                  GRN = Greensprings Historic Site
#'                                  GET = Gettysburg National Battlefield}
#'   \item{transect}{transect number, unique to each transect}
#'   \item{plantID}{plant identification number used to track growth and survival, unique to each plant within a transect}
#'   \item{h_apical}{plant height measured in cm at the beginning of the growing season (June)}
#'   \item{h_apical.next}{plant height measured in cm at the end of the growing season (the following September)}
#'   \item{stem_width}{width taken at the base of the stem measured in mm}
#'   \item{alive_june}{factor indicating whether the plant was alive in June}
#'   \item{surv}{factor indicating whether the plant survived to reproduce}
#'   \item{fec.flower}{factor indicating whether the plant produced flowers}
#'   \item{N_pods}{the number of seed pods produced}
#'   \item{seedling}{factor indicating seedling status}
#'   \item{munched}{factor indicating if the plant experienced herbivory}
#'   \item{per.leaves.dam}{percent of leaves that had at least 1 cm^2 of tissue removed by herbivores; herbivory measurement used in 2013 and then transformed to herb_avg to match subsequent years}
#'   \item{herb_avg}{average herbivory score experienced: 0 = intact
#'                                                        1 = 1-5% removed
#'                                                        2 = 6-24% removed
#'                                                        3 = 25-50% removed
#'                                                        4 = 51-75% removed
#'                                                        5 = 76-99% removed
#'                                                        6 = only petiole remaining}
#'   ...
#' }
"stemdata"

#' Seed data!
#'
"seeddata"

#'
#' demdata
#'
#' #' A dataset containing demographic information for milkweed plants from 7 sites across VA and PA from 2013-2017... truncated to 2015
#'
#' @format A data frame with 5829 rows and 16 variables:
#' \describe{
#'   \item{year}{collection year}
#'   \item{site}{3-digit site code: BLD = Blandy Experimental Farm (actually two sites, 1 indicating the area that is unburned, 2 the area that is burned)
#'                                  SKY = Sky Meadows State Park
#'                                  PWR = Presquile National Wildlife Refuge
#'                                  YTB = Yorktown National Battlefield
#'                                  GRN = Greensprings Historic Site
#'                                  GET = Gettysburg National Battlefield}
#'   \item{transect}{transect number, unique to each transect}
#'   \item{plantID}{plant identification number used to track growth and survival, unique to each plant within a transect}
#'   \item{h_apical}{plant height measured in cm at the beginning of the growing season (June)}
#'   \item{h_apical.next}{plant height measured in cm at the end of the growing season (the following September)}
#'   \item{stem_width}{width taken at the base of the stem measured in mm}
#'   \item{alive_june}{factor indicating whether the plant was alive in June}
#'   \item{surv}{factor indicating whether the plant survived to reproduce}
#'   \item{fec.flower}{factor indicating whether the plant produced flowers}
#'   \item{N_pods}{the number of seed pods produced}
#'   \item{seedling}{factor indicating seedling status}
#'   \item{munched}{factor indicating if the plant experienced herbivory}
#'   \item{per.leaves.dam}{percent of leaves that had at least 1 cm^2 of tissue removed by herbivores; herbivory measurement used in 2013 and then transformed to herb_avg to match subsequent years}
#'   \item{herb_avg}{average herbivory score experienced: 0 = intact
#'                                                        1 = 1-5% removed
#'                                                        2 = 6-24% removed
#'                                                        3 = 25-50% removed
#'                                                        4 = 51-75% removed
#'                                                        5 = 76-99% removed
#'                                                        6 = only petiole remaining}
#'   ...
#' }
"dem.data"

#' traitdata
#'
#' A dataset containing leaf chemistry data collected from milkweed plants in the field using spectroradiometry.
#'
#' @format A data frame with 5829 rows and 16 variables:
#' \describe{
#'   \item{year}{collection year}
#'   \item{site}{3-digit site code: BLD = Blandy Experimental Farm (actually two sites, 1 indicating the area that is unburned, 2 the area that is burned)
#'                                  SKY = Sky Meadows State Park
#'                                  PWR = Presquile National Wildlife Refuge
#'                                  YTB = Yorktown National Battlefield
#'                                  GRN = Greensprings Historic Site
#'                                  GET = Gettysburg National Battlefield}
#'   \item{transect}{transect number, unique to each transect}
#'   \item{plantID}{plant identification number used to track growth and survival, unique to each plant within a transect}
#'   \item{h_apical}{plant height measured in cm at the beginning of the growing season (June)}
#'   \item{h_apical.next}{plant height measured in cm at the end of the growing season (the following September)}
#'   \item{stem_width}{width taken at the base of the stem measured in mm}
#'   \item{alive_june}{factor indicating whether the plant was alive in June}
#'   \item{surv}{factor indicating whether the plant survived to reproduce}
#'   \item{fec.flower}{factor indicating whether the plant produced flowers}
#'   \item{N_pods}{the number of seed pods produced}
#'   \item{seedling}{factor indicating seedling status}
#'   \item{munched}{factor indicating if the plant experienced herbivory}
#'   \item{per.leaves.dam}{percent of leaves that had at least 1 cm^2 of tissue removed by herbivores; herbivory measurement used in 2013 and then transformed to herb_avg to match subsequent years}
#'   \item{herb_avg}{average herbivory score experienced: 0 = intact
#'                                                        1 = 1-5% removed
#'                                                        2 = 6-24% removed
#'                                                        3 = 25-50% removed
#'                                                        4 = 51-75% removed
#'                                                        5 = 76-99% removed
#'                                                        6 = only petiole remaining}
#'   ...
#' }
"trait.data"

#' fulldata
#'
#' A dataset containing leaf chemistry data collected from milkweed plants in the field using spectroradiometry, combined with their respective demographic data.
#'
#' @format A data frame with 5829 rows and 16 variables:
#' \describe{
#'   \item{year}{collection year}
#'   \item{site}{3-digit site code: BLD = Blandy Experimental Farm (actually two sites, 1 indicating the area that is unburned, 2 the area that is burned)
#'                                  SKY = Sky Meadows State Park
#'                                  PWR = Presquile National Wildlife Refuge
#'                                  YTB = Yorktown National Battlefield
#'                                  GRN = Greensprings Historic Site
#'                                  GET = Gettysburg National Battlefield}
#'   \item{transect}{transect number, unique to each transect}
#'   \item{plantID}{plant identification number used to track growth and survival, unique to each plant within a transect}
#'   \item{h_apical}{plant height measured in cm at the beginning of the growing season (June)}
#'   \item{h_apical.next}{plant height measured in cm at the end of the growing season (the following September)}
#'   \item{stem_width}{width taken at the base of the stem measured in mm}
#'   \item{alive_june}{factor indicating whether the plant was alive in June}
#'   \item{surv}{factor indicating whether the plant survived to reproduce}
#'   \item{fec.flower}{factor indicating whether the plant produced flowers}
#'   \item{N_pods}{the number of seed pods produced}
#'   \item{seedling}{factor indicating seedling status}
#'   \item{munched}{factor indicating if the plant experienced herbivory}
#'   \item{per.leaves.dam}{percent of leaves that had at least 1 cm^2 of tissue removed by herbivores; herbivory measurement used in 2013 and then transformed to herb_avg to match subsequent years}
#'   \item{herb_avg}{average herbivory score experienced: 0 = intact
#'                                                        1 = 1-5% removed
#'                                                        2 = 6-24% removed
#'                                                        3 = 25-50% removed
#'                                                        4 = 51-75% removed
#'                                                        5 = 76-99% removed
#'                                                        6 = only petiole remaining}
#'   ...
#' }
"full.data"

#' mapdata
#'
#' A dataset containing x and y location data for a sample of transects in the study.
#'
#' @format A data frame with 5829 rows and 16 variables:
#' \describe{
#'   \item{year}{collection year}
#'   \item{site}{3-digit site code: BLD = Blandy Experimental Farm (actually two sites, 1 indicating the area that is unburned, 2 the area that is burned)
#'                                  SKY = Sky Meadows State Park
#'                                  PWR = Presquile National Wildlife Refuge
#'                                  YTB = Yorktown National Battlefield
#'                                  GRN = Greensprings Historic Site
#'                                  GET = Gettysburg National Battlefield}
#'   \item{transect}{transect number, unique to each transect}
#'   \item{plantID}{plant identification number used to track growth and survival, unique to each plant within a transect}
#'   \item{h_apical}{plant height measured in cm at the beginning of the growing season (June)}
#'   \item{h_apical.next}{plant height measured in cm at the end of the growing season (the following September)}
#'   \item{stem_width}{width taken at the base of the stem measured in mm}
#'   \item{alive_june}{factor indicating whether the plant was alive in June}
#'   \item{surv}{factor indicating whether the plant survived to reproduce}
#'   \item{fec.flower}{factor indicating whether the plant produced flowers}
#'   \item{N_pods}{the number of seed pods produced}
#'   \item{seedling}{factor indicating seedling status}
#'   \item{munched}{factor indicating if the plant experienced herbivory}
#'   \item{per.leaves.dam}{percent of leaves that had at least 1 cm^2 of tissue removed by herbivores; herbivory measurement used in 2013 and then transformed to herb_avg to match subsequent years}
#'   \item{herb_avg}{average herbivory score experienced: 0 = intact
#'                                                        1 = 1-5% removed
#'                                                        2 = 6-24% removed
#'                                                        3 = 25-50% removed
#'                                                        4 = 51-75% removed
#'                                                        5 = 76-99% removed
#'                                                        6 = only petiole remaining}
#'   ...
#' }
"map.data"

#' buddata
#'
#' A dataset containing budlings in 2015 produced per stem in 2014, as well as 2014 demographic data and 2015 traitdata averaged on the 1m transect block level.
#'
#' @format A data frame with 5829 rows and 16 variables:
#' \describe{
#'   \item{year}{collection year}
#'   \item{site}{3-digit site code: BLD = Blandy Experimental Farm (actually two sites, 1 indicating the area that is unburned, 2 the area that is burned)
#'                                  SKY = Sky Meadows State Park
#'                                  PWR = Presquile National Wildlife Refuge
#'                                  YTB = Yorktown National Battlefield
#'                                  GRN = Greensprings Historic Site
#'                                  GET = Gettysburg National Battlefield}
#'   \item{transect}{transect number, unique to each transect}
#'   \item{plantID}{plant identification number used to track growth and survival, unique to each plant within a transect}
#'   \item{h_apical}{plant height measured in cm at the beginning of the growing season (June)}
#'   \item{h_apical.next}{plant height measured in cm at the end of the growing season (the following September)}
#'   \item{stem_width}{width taken at the base of the stem measured in mm}
#'   \item{alive_june}{factor indicating whether the plant was alive in June}
#'   \item{surv}{factor indicating whether the plant survived to reproduce}
#'   \item{fec.flower}{factor indicating whether the plant produced flowers}
#'   \item{N_pods}{the number of seed pods produced}
#'   \item{seedling}{factor indicating seedling status}
#'   \item{munched}{factor indicating if the plant experienced herbivory}
#'   \item{per.leaves.dam}{percent of leaves that had at least 1 cm^2 of tissue removed by herbivores; herbivory measurement used in 2013 and then transformed to herb_avg to match subsequent years}
#'   \item{herb_avg}{average herbivory score experienced: 0 = intact
#'                                                        1 = 1-5% removed
#'                                                        2 = 6-24% removed
#'                                                        3 = 25-50% removed
#'                                                        4 = 51-75% removed
#'                                                        5 = 76-99% removed
#'                                                        6 = only petiole remaining}
#'   ...
#' }
"bud.data"
