import os
import sys
import math
import random
from decimal import *
import numpy as np
import scipy as sp

# # To create the data file:
# SG_dataset=read.csv('../../data/Salix_example/cooccur_interact_galler_salix.csv',row.names=1)
# SG_web=read.csv('figure_generation/binary_prior_web_SG.csv',row.names=1)
# deg_dist_Salix=rowSums(SG_web)/ncol(SG_web)
# deg_dist_galler=colSums(SG_web)/nrow(SG_web)
# # Interaction probabilities are the product of plant and galler probabilities
# sg_int_probs=as.numeric(deg_dist_galler%*%t(deg_dist_Salix))
# SG_dataset$post.mean<-1
# SG_dataset$post.sd<-0
# for(r in 1:nrow(SG_dataset)){
#   if(SG_dataset$interact[r]==0){
#     pars=calculate_parameters(sg_int_probs,SG_dataset$cooccur[r],0)
#     dist=calculate_distribution(pars)
#     SG_dataset[r,7:8]=dist
#   }
# }
# write.table(SG_dataset,file='../../data/Salix_example/Salix_Galler/posterior_probabilities.tsv')

# GP_dataset=read.csv('../../data/Salix_example/cooccur_interact_galler_parasit.csv',row.names=1)
# GP_web=read.csv('figure_generation/binary_prior_web_para_only.csv',row.names=1)
# deg_dist_galler=rowSums(GP_web)/ncol(GP_web)
# deg_dist_paras=colSums(GP_web)/nrow(GP_web)
# # Interaction probabilities are the product of plant and galler probabilities
# gp_int_probs=as.numeric(deg_dist_paras%*%t(deg_dist_galler))
# GP_dataset$post.mean<-1
# GP_dataset$post.sd<-0
# for(r in 1:nrow(GP_dataset)){
#   if(GP_dataset$interact[r]==0){
#     pars=calculate_parameters(gp_int_probs,GP_dataset$cooccur[r],0)
#     dist=calculate_distribution(pars)
#     GP_dataset[r,7:8]=dist
#   }
# }
# write.table(GP_dataset,file='../../data/Salix_example/Galler_Parasitoid/posterior_probabilities.tsv')


#   numpy.random.binomial(n,p,size) draws [size] samples from a binomial distribution
#   % Each observation is a Bernoulli, or a binomial with n=1. Use p=mean(posterior).

def read_data(datafile,site):
  ints={'SG':set(),'GP':set()}
  Salix_sp=set()
  galler_sp=set()
  para_sp=set()
  f=open(datafile,'r')
  for line in f:
    newline=line.split('\r')[0]
    if site=='Zillertal':
      if newline.split(',')[0]!='REARING_NUMBER':
        Salix_ID=newline.split(',')[10]
        galler=newline.split(',')[11]
        para=newline.split(',')[12] 
        n_galls=int(newline.split(',')[13])
        n_galls_parasit=int(newline.split(',')[14])
        Salix=newline.split(',')[15]
    else:
      if newline.split(',')[1]!='REARING_NUMBER':
        Salix=newline.split(',')[10]
        galler=newline.split(',')[11]
        para=newline.split(',')[12] 
        n_galls=int(newline.split(',')[13])
        n_galls_parasit=int(newline.split(',')[14])
    if newline.split(',')[0]!='REARING_NUMBER' and newline.split(',')[1]!='REARING_NUMBER':
      Salix_sp.add(Salix)
      galler_sp.add(galler)
      ints['SG'].add((Salix,galler))
      # Some galls were never parasitized
      if para!='none':
        para_sp.add(para)
        ints['GP'].add((galler, para))
  f.close()
  return ints, Salix_sp, galler_sp, para_sp

def write_prior_web(ints,resources,consumers,flavour,site):
  if flavour=='SG':
    f=open('../../data/Salix_example/Zillis/'+site+'_SG_prior.csv','w')
  else:
    f=open('../../data/Salix_example/Zillis/'+site+'_GP_prior.csv','w')
  f.write(','+','.join(sorted(consumers))+'\n')
  for res in sorted(resources):
    f.write(res)
    for con in sorted(consumers):
      if (res,con) in ints:
        f.write(',1')
      else:
        f.write(',0')
    f.write('\n')

  f.close()

def main():

  for site in ['Zillis','Zillertal']:
    print site
    datafile='../../data/Salix_example/Zillis/'+site+'_web.csv'

    ints, Salix_sp, galler_sp, para_sp=read_data(datafile,site)

    write_prior_web(ints['SG'],Salix_sp,galler_sp,'SG',site)
    write_prior_web(ints['GP'],galler_sp,para_sp,'GP',site)
  
if __name__ == '__main__':
  main()

