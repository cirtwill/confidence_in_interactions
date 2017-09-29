source('handy_R_functions.R')

# Read in and process the raw data
prior_web_data=read.csv('../data/Salix_example/tree_level_interaxn_all_plants_traits_size.csv')
# Remove extra data on parasitoids, etc.
prior_web_data[,31:57]<-NULL
prior_web_data[,4:26]<-NULL
prior_web_data[,2]<-NULL
# Build a web with interactions=1 for any galler observed on any genotype
prior_web=matrix(nrow=26,ncol=4)
for(r in 1:length(levels(prior_web_data$Genotype))){
  gen=levels(prior_web_data$Genotype)[r]
  subset=prior_web_data[which(prior_web_data$Genotype==gen),]
  for(c in 1:4){
    if(sum(subset[,2+c])>0){
      prior_web[r,c]=1
    } else {
      prior_web[r,c]=0
    }
  }        
}


# Now get the degree distributions from the web
deg_dist_Salix=rowSums(prior_web)/ncol(prior_web)
# Remove one genotype that never interacted
deg_dist_Salix<-deg_dist_Salix[which(rowSums(prior_web)>0)]
deg_dist_galler=colSums(prior_web)/nrow(prior_web)
# Interaction probabilities are the product of plant and galler probabilities
int_probs=as.numeric(deg_dist_galler%*%t(deg_dist_Salix))


fulldata=read.table('../data/Salix_example/cooccur_interact_galler_salix.csv',header=TRUE,sep=',',row.names=1)
fulldata$post_mean=0
fulldata$post_SD=0
for(r in 1:nrow(fulldata)){
      n=fulldata$cooccur[r]
      k=fulldata$interact[r]
      post_params=calculate_distribution(calculate_parameters(int_probs,n,k))
      fulldata[r,'post_mean']=post_params[[1]]
      fulldata[r,'post_SD']=post_params[[2]]
      print(r) # There are about 5000 rows.
}

write.table(fulldata,file='../data/Salix_example/Salix_Galler/posterior_probabilities.tsv',sep='\t')

