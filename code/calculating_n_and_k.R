setwd("C:/Users/kewo0002/confidence_in_interactions/data/Salix_example")

#######################################################
#### Create list of all combinations of galler-salix and 
#### galler-parasitoid and how often they interact and cooccur
#######################################################

library("plyr")
library("igraph")

#### FUNCTIONS ####

# Function to calculate the number of times each galler-salix combination
# cooccurs and interacts
calc_cooccurence <- function(galler_salix,dataset){
  # the first variable in the galler_salix dataframe is the salix
  salix <- galler_salix[1]
  # the second is the galler
  galler <- galler_salix[2]
  
  # get only those rows where those two species both appear (and therefore interact)
  inter_subs <- subset(dataset,Rgaller==galler & Salix == salix)
  
  # print(inter_subs)
  
  # the number of sites they interact with is the number of unique sites in that list
  interact <- length(unique(inter_subs$site_ID))
  
  if(interact>1){print(head(inter_subs))}
  
  # run through each site and calculate co-occurrences
  cooccurence <- ddply(dataset, .(site_ID), count_cooccur, salix=salix, galler=galler)
  # if(sum(cooccurence$V1) == 0 ){print("no cooccurence")}
  
  # return the salix, galler, number of co-occurrences and number of interactions
  return(c(salix,galler,sum(cooccurence$V1), interact))
}


# function to count co-occurrences (for galler-salix)
count_cooccur <- function(site_network, salix, galler){
  #print(paste(salix,galler))
  
  # for each network (passed to the function), if both salix and galler appear
  # then it's a cooccurrence, return one. Otherwise return 0 (no co-occurrence)
  if(salix %in% site_network$Salix & galler %in% site_network$Rgaller){  
    #print(site_network)
    #print(paste("cooccurrence,", salix, galler, " site_ID = ", site_network[1,"site_ID"]))
    return(1)
  }else{return(0)}
}


# calculate cooccurrences for parasitoid-galler (data is in a matrix so it has to be done differently)
matr_calc_cooccur <- function(paras_galler, dataset){
  # set the galler and parasitoid we are interested in
  para <- paras_galler[2]
  galler <- paras_galler[1]
  
  # take just the data which has that galler
  intlist <- subset(dataset, Rgaller==galler)
  
  # take just the data where the galler interacts with the parasitoid of interest
  intlistboth <- intlist[which(intlist[,as.character(para)]>0),]

  # count how many sites that's at so we can get a number of interactions
  int.count <- length(unique(intlistboth$site_ID))
  
  # Now run through all the sites and check for cooccurrences
  site_list <- unique(dataset$site_ID)
  cooccur <- sapply(site_list, matr_count_cooccur, paras_galler = paras_galler, dataset = dataset)
  
  # cooccur is a list of 1s and 0s, so sum gives the total number of cooccurrences
  return(c(galler, para, sum(cooccur), int.count))
}

# count cooccurences for galler-parasitoid data (in a matrix so has to be done differently)
matr_count_cooccur <- function(site,paras_galler,dataset)
{
  # set our galler and parasitoid
  para <- paras_galler[2]
  galler <- paras_galler[1]
  
  # just take the data from that site
  site.matrix <- dataset[which(dataset$site_ID == site),]
  
  # if both the galler and the parasitoid appear, then they cooccur, return 1.
  if(galler %in% site.matrix$Rgaller & (sum(site.matrix[,as.character(para)]) > 0)){
    return(1)
  }else{
    # otherwise they don't cooccur, return 0
    return(0)
  }
}
#### END OF FUNCTIONS ####


#### GALLER AND SALIX ####

# read in full matrix of all interactions
interactions <- read.csv("\\\\storage.slu.se/Home$/kewo0002/My Documents/Projects/Gallers/original.data.from.paper/Salix_webs.csv",stringsAsFactors=FALSE)

# for just the galler-salix interactions, we only need certain columns
galler_salix <- interactions[,c("X","X.5","X.9","X.14","X.18")]
# the columns should be named
colnames(galler_salix) <- c("Year","site","Salix","Rgaller","ngalls")
# the first 11 rows are empty, we should cut them out
galler_salix_clean <- galler_salix[11:nrow(galler_salix),]

# create a data frame with each possible combination of galler and salix
galler_salix_combos <- expand.grid(unique(galler_salix_clean$Salix),unique(galler_salix_clean$Rgaller))
# give them the same names as the data
colnames(galler_salix_combos) <- c("Salix","Rgaller")

# some sites were visited more than once, so create a new variable that combines site and year
#(can't use date because some sites were sampled over multiple days)
galler_salix_clean$site_time <- as.factor(paste0(galler_salix_clean$site,galler_salix_clean$Year))


#TEST
# there should be 641 unique site-year combinations
if(length(unique(galler_salix_clean$site_time))==641){print("pass")}else{print(paste("Something wrong, should have 641 sites visits, actually have ", length(unique(galler_salix_clean$site_time))))}

# give each site-year combination a digit ID (just makes it easier to work with)
galler_salix_clean$site_ID <- mapvalues(galler_salix_clean$site_time, from = levels(galler_salix_clean$site_time), to = 1:length(levels(galler_salix_clean$site_time)))


# for each combination of galler and salix (running down the rows in the data frame of combinations)
# calculate cooccurrences and interactions
output <- apply(galler_salix_combos, 1, calc_cooccurence, dataset=galler_salix_clean)


#trial runs (on a subset of the data)
#shortlist <- galler_salix_clean[which(galler_salix_clean$site %in% site_list[1:20]),]
#output <- apply(galler_salix_combos[1,], 1, calc_cooccurence, dataset=shortlist)


# transpose the output 
toutput <- as.data.frame(t(output))
# give it column names
colnames(toutput) <- c("Salix","Rgaller","cooccur","interact")
# make the counts numeric
toutput$cooccur <- as.numeric(levels(toutput$cooccur))[toutput$cooccur]
toutput$interact <- as.numeric(levels(toutput$interact))[toutput$interact]


# work out what range of values we've got
range(toutput$cooccur)
range(toutput$interact)

# add numeric IDs to make plotting easier
toutput$salix_ID <- mapvalues(toutput$Salix, from = levels(toutput$Salix), to = 1:length(levels(toutput$Salix)))
toutput$galler_ID <- mapvalues(toutput$Rgaller, from = levels(toutput$Rgaller), to = 1:length(levels(toutput$Rgaller)))

# get an idea of what the data looks like
#plot(toutput$cooccur,toutput$interact)

write.csv(toutput,file="cooccur_interact_galler_salix.csv",quote=F)



#### GALLER AND PARASITOID ####

# create a data frame of galler-parasitoid interactions
int.matrix <- read.csv("\\\\storage.slu.se/Home$/kewo0002/My Documents/Projects/Gallers/data.dataframes/df_interaction_matrix.csv",stringsAsFactors=FALSE)
# set the nas to 0
int.matrix[is.na(int.matrix)] <- 0
# create a variable for each combination of site and time
int.matrix$site_time <- as.factor(paste0(int.matrix$site,int.matrix$Year.of.coll))
# and then give it a number to make it easier to work with
int.matrix$site_ID <- mapvalues(int.matrix$site_time, from = levels(int.matrix$site_time), to = 1:length(levels(int.matrix$site_time)))


#TEST
# there should be 641 unique site-year combinations
if(length(unique(int.matrix$site_time))==641){print("pass")}else{print(paste("Something wrong, should have 641 sites visits, actually have ", length(unique(int.matrix$site_time))))}

# create a data frame with each possible combination of galler and paras
paras_galler_combos<- expand.grid(unique(int.matrix$Rgaller),unique(colnames(int.matrix)[4:(ncol(int.matrix)-2)]))
# give them the same names as the data
colnames(paras_galler_combos) <- c("Rpara","Rgaller")

# I broke this into bits because it took so long to run. 
# But otherwise use output <- apply(....)
output1 <- apply(paras_galler_combos[1:600,],1,matr_calc_cooccur, dataset = int.matrix)
output2 <- apply(paras_galler_combos[601:1209,],1,matr_calc_cooccur, dataset = int.matrix)
output3 <- apply(paras_galler_combos[1210:1800,],1,matr_calc_cooccur, dataset = int.matrix)
output4 <- apply(paras_galler_combos[1801:2400,],1,matr_calc_cooccur, dataset = int.matrix)
output5 <- apply(paras_galler_combos[2401:3000,],1,matr_calc_cooccur, dataset = int.matrix)
output6 <- apply(paras_galler_combos[3001:3600,],1,matr_calc_cooccur, dataset = int.matrix)
output7 <- apply(paras_galler_combos[3601:4200,],1,matr_calc_cooccur, dataset = int.matrix)
output8 <- apply(paras_galler_combos[4201:4800,],1,matr_calc_cooccur, dataset = int.matrix)
output9 <- apply(paras_galler_combos[4801:5400,],1,matr_calc_cooccur, dataset = int.matrix)
output10 <- apply(paras_galler_combos[5401:6000,],1,matr_calc_cooccur, dataset = int.matrix)
output11 <- apply(paras_galler_combos[6001:6600,],1,matr_calc_cooccur, dataset = int.matrix)
output12 <- apply(paras_galler_combos[6601:7200,],1,matr_calc_cooccur, dataset = int.matrix)
output13 <- apply(paras_galler_combos[7201:7800,],1,matr_calc_cooccur, dataset = int.matrix)
output14 <- apply(paras_galler_combos[7801:8400,],1,matr_calc_cooccur, dataset = int.matrix)
output15 <- apply(paras_galler_combos[8401:9000,],1,matr_calc_cooccur, dataset = int.matrix)
output16 <- apply(paras_galler_combos[9001:9600,],1,matr_calc_cooccur, dataset = int.matrix)
output17 <- apply(paras_galler_combos[9601:12096,],1,matr_calc_cooccur, dataset = int.matrix)

# stick all the bits together
output <- cbind(output1,output2,output3,output4,output5,output6,output7,output8,output9,output10,output11,output12,output13,output14,output15,output16,output17)

# this takes ages to run, but does it in one go
#output <- apply(paras_galler_combos,1,matr_calc_cooccur, dataset = int.matrix)
 
#trial runs (on a subset of the data)
#shortlist <- int.matrix[which(int.matrix$site %in% unique(int.matrix$site)[1:20]),]
#output <- apply(paras_galler_combos, 1, matr_calc_cooccur, dataset=int.matrix)


# transpose the output 
toutput <- as.data.frame(t(output))
# give it column names
colnames(toutput) <- c("Rgaller","Rpara","cooccur","interact")
# make the counts numeric
toutput$cooccur <- as.numeric(levels(toutput$cooccur))[toutput$cooccur]
toutput$interact <- as.numeric(levels(toutput$interact))[toutput$interact]


# work out what range of values we've got
range(toutput$cooccur)
range(toutput$interact)

# add a numeric ID to make plotting easier
toutput$para_ID <- mapvalues(toutput$Rpara, from = levels(toutput$Rpara), to = 1:length(levels(toutput$Rpara)))
toutput$galler_ID <- mapvalues(toutput$Rgaller, from = levels(toutput$Rgaller), to = 1:length(levels(toutput$Rgaller)))

#take a look at the data
#plot(toutput$cooccur,toutput$interact)

write.csv(toutput,file="cooccur_interact_galler_parasit.csv",quote=F)




#####################
# Everything below this was playing around trying to get it to work
# and is mostly a mess. Probably best to ignore :)
#####################

##Check 
mcgaller <- "Hplicl"
mcsalix <- "lapponum"
  
most_common_combo <- subset(galler_salix_clean, Rgaller==mcgaller & Salix == mcsalix)
interact <- length(unique(most_common_combo$site_ID))
cooccur <- by(galler_salix_clean, site_ID, function(x) if(mcgaller %in% x$Rgaller & mcsalix %in% x$Salix){return(1)}else{return(0)} )




#Trials
for(site in unique(shortlist$site_ID)){
  print(site)
  site_network <- subset(shortlist,site_ID == site)
  occurences <- count_cooccur(site_network, salix = "elaeagnos", galler = "Oelaea")}
#  print(c(cooccur, interact_site))
#output <- list(cooccur,interact_site)
#print(cooccur)
#return(c(cooccur)



################# Trials
#functions
count_interact <- function(row.to.check,salix,galler){
  if(row.to.check["Salix"] == salix && row.to.check["Rgaller"] == galler) check <- 1
  else check <- 0
  check
}
count_cooccur <- function(site_network, salix, galler){
  print(site_network)
  cooccur <- 0
  interact_site <- 0
  if(salix %in% site_network$Salix & galler %in% site_network$Rgaller){
    #print(paste0("cooccurrence,", salix, galler))
    cooccur <- 1
  }
  
 # if(sum(apply(site_network,1,count_interact, salix = salix, galler = galler)) > 0)
 # { interact_site <- 1}
#  print(c(cooccur, interact_site))
 # output <- list(cooccur,interact_site)
  print(cooccur)
  return(cooccur)
}


interactions <- read.csv("../original.data.from.paper/Salix_webs.csv")

galler_salix <- interactions[,c("X.5","X.9","X.14","X.18")]
colnames(galler_salix) <- c("site","Salix","Rgaller","ngalls")
galler_salix_clean <- galler_salix[11:nrow(galler_salix),]

site_list <- unique(galler_salix_clean$site)
galler_salix_clean$site_ID <- mapvalues(galler_salix_clean$site, from = levels(galler_salix_clean$site), to = 1:length(levels(galler_salix_clean$site)))

shortlist <- galler_salix_clean[which(galler_salix_clean$site %in% site_list[1:20]),]


inter_list <- data.frame()
row <- 1
gallerlist <- unique(dataset$Rgaller)
salixlist <- unique(dataset$Salix)

calc_cooccurence <- function(dataset,galler,salix){
  cooccurence <- ddply(dataset, .(site_ID), count_cooccur, salix=salix, galler=galler)
  cooccurence
}

whateverthisis <- sapply(gallerlist, function(galler1) sapply(salixlist, function(salix1) calc_cooccurence, dataset=shortlist, salix = salix1, galler=galler1))

calc.n.k <- function(dataset){
  for(galler in unique(dataset$Rgaller)){
    for(salix in unique(dataset$Salix)){      
      inter_list[row,"galler"] <- galler
      inter_list[row,"salix"] <- salix
      cooccurence <- ddply(dataset, .(site_ID), count_cooccur, salix=salix, galler=galler)
      
      #cooccur <- 0
      #interact_site <- 0
      #for(site in unique(dataset$site)){   
            # use sapply to put the data back to a matrix
      cooccur_matrix<- t(sapply(cooccurence, I))
      # make a data frame
      cooccur_df <- as.data.frame(cooccur_matrix)  
      
      
        #interact_count <- 0
        #site_network <- dataset[which(dataset$site == site),]
        
       # if(salix %in% site_network$Salix & galler %in% site_network$Rgaller){
       #   cooccur <- cooccur + 1 
       # }
        
        #if(sum(apply(site_network,1,count_interact, salix = salix, galler = galler)) > 0)
        #{ interact_site <- interact_site + 1}
        
        #for(rownum in 1:nrow(site_network))
        #{
         # if(site_network[rownum,"Salix"] == salix & site_network[rownum,"Rgaller"] == galler){
          #  interact_count <- interact_count + 1
            
          #}
           
        
        #}
        #if(interact_count > 0){ interact_site <- interact_site + 1}
       
      }
      inter_list[row,"n"] <- cooccur
      inter_list[row,"k"] <- interact_site
      row <- row + 1  
    }
    
    
  }
  inter_list
}

trial_inter_list <- calc.n.k(shortlist)

write.csv(trial_inter_list,file="galler_salix_nk_20sites.csv",quote=FALSE,)
