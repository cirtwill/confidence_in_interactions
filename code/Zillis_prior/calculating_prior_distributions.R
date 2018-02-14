source('../handy_R_functions.R')
# This prior isn't working with fitdistr, somehow. :(

# To create the data file:
for(subweb in  c('Zillis','Zillertal')){
  SG_dataset=read.csv('../../data/Salix_example/cooccur_interact_galler_salix.csv',row.names=1)
  SG_web=read.csv(paste0('../../data/Salix_example/Zillis/',subweb,'_SG_prior.csv',sep=''),row.names=1)
  deg_dist_Salix=rowSums(SG_web)/ncol(SG_web)
  deg_dist_galler=colSums(SG_web)/nrow(SG_web)
  # Interaction probabilities are the product of plant and galler probabilities
  sg_int_probs=as.numeric(deg_dist_galler%*%t(deg_dist_Salix))
  SG_dataset$post.mean<-1
  SG_dataset$post.sd<-0
  for(r in 1:nrow(SG_dataset)){
    if(SG_dataset$interact[r]==0){
      pars=calculate_parameters(sg_int_probs,SG_dataset$cooccur[r],0)
      dist=calculate_distribution(pars)
      SG_dataset[r,7:8]=dist
    }
  }
  write.table(SG_dataset,file=paste0('../../data/Salix_example/Zillis/Salix_Galler/posterior_probabilities_',subweb,'.tsv',sep=''))
  print(subweb)

  GP_dataset=read.csv('../../data/Salix_example/cooccur_interact_galler_parasit.csv',row.names=1)
  GP_web=read.csv(paste0('../../data/Salix_example/Zillis/',subweb,'_GP_prior.csv',sep=''),row.names=1)
  deg_dist_galler=rowSums(GP_web)/ncol(GP_web)
  deg_dist_paras=colSums(GP_web)/nrow(GP_web)
  # Interaction probabilities are the product of plant and galler probabilities
  gp_int_probs=as.numeric(deg_dist_paras%*%t(deg_dist_galler))
  GP_dataset$post.mean<-1
  GP_dataset$post.sd<-0
  for(r in 1:nrow(GP_dataset)){
    if(GP_dataset$interact[r]==0){
      pars=calculate_parameters(gp_int_probs,GP_dataset$cooccur[r],0)
      dist=calculate_distribution(pars)
      GP_dataset[r,7:8]=dist
    }
  }
  write.table(GP_dataset,file=paste0('../../data/Salix_example/Zillis/Galler_Parasitoid/posterior_probabilities_',subweb,'.tsv',sep=''))

}

