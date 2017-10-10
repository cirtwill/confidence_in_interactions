setwd("C:/Users/kewo0002/confidence_in_interactions/data/Salix_example") 

###############################################################
#### Build the galler-parasitoid web from Balbour 2016
#### to be used as the posterior for Kopelke's data
###############################################################

# read in the data
prior_web_data=read.csv('tree_level_interaxn_all_plants_traits_size.csv')
# just keep what we need
prior_web_data <- prior_web_data[,8:26]

# Build a web with interactions=1 for any parasitoid on any galler
prior_web=data.frame()

# each column contains all interactions for a particular parasitoid-galler combination
for(column in 1:ncol(prior_web_data)){
  # the column name gives the species which makes the gall, followed by the species emerging from the gall
  # they need to be separated 
  galler_emerged=strsplit(colnames(prior_web_data)[column], "_")
  galler = galler_emerged[[1]][1]
  emerged = galler_emerged[[1]][2]
  
  # the strength of the interaction is the sum of the column (number of galls across all branches sampled)
  linkstrength = sum(prior_web_data[,column])
  # put that into a dataframe
  prior_web[as.character(galler),as.character(emerged)] <- linkstrength
} 

# Make the naS into 0
prior_web[is.na(prior_web)] <- 0

# some of the emergers were the gallers, remove those columns if we're only
# interested in galler-parasitoid interactions
prior_web_para_only <- prior_web[,c(2:6,8:9,12)]
write.csv(prior_web_para_only, file="prior_web_para_only.csv", quote=F)
