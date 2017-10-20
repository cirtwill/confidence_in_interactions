import os
import sys
import math
import random
from decimal import *
import numpy as np
import scipy as sp

# # To create the data file:
# SG_dataset=read.csv('../data/Salix_example/cooccur_interact_galler_salix.csv',row.names=1)
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
# write.table(SG_dataset,file='../data/Salix_example/Salix_Galler/posterior_probabilities.tsv')

# GP_dataset=read.csv('../data/Salix_example/cooccur_interact_galler_parasit.csv',row.names=1)
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
# write.table(GP_dataset,file='../data/Salix_example/Galler_Parasitoid/posterior_probabilities.tsv')


#   numpy.random.binomial(n,p,size) draws [size] samples from a binomial distribution
#   % Each observation is a Bernoulli, or a binomial with n=1. Use p=mean(posterior).

def read_in_data(datafile):
  pdict={}
  f=open(datafile,'r')
  for line in f:
    if line.split()[0] not in ['"Salix"','"Rgaller"']:
      if len(line.split())==9:
        plant=line.split()[1][1:-1]
        galler=line.split()[2][1:-1]
        occurs=int(line.split()[3])
        interacts=int(line.split()[4])
        plantID=int(line.split()[5])
        gallerID=int(line.split()[6])
        postmean=float(line.split()[7])
        postSD=float(line.split()[8])
        sppair=(plant,galler)
        pdict[sppair]=((postmean,postSD))
      elif len(line.split())==12:
        plant='.'.join([line.split()[1][1:],line.split()[2],line.split()[3],'?']) # 4 is always '?'
        galler=line.split()[5][1:-1]
        occurs=int(line.split()[6])
        interacts=int(line.split()[7])
        plantID=int(line.split()[8])
        gallerID=int(line.split()[9])
        postmean=float(line.split()[10])
        postSD=float(line.split()[11])
        sppair=(plant,galler)
        pdict[sppair]=((postmean,postSD))
      elif len(line.split())==11:
        plant='.'.join([line.split()[1][1:],line.split()[2],line.split()[3][:-1]]) 
        galler=line.split()[4][1:-1]
        occurs=int(line.split()[5])
        interacts=int(line.split()[6])
        plantID=int(line.split()[7])
        gallerID=int(line.split()[8])
        postmean=float(line.split()[9])
        postSD=float(line.split()[10])
        sppair=(plant,galler)
        pdict[sppair]=((postmean,postSD))
      elif len(line.split())==10:
        plant='.'.join([line.split()[1][1:],line.split()[2][:-1]]) 
        galler=line.split()[3][1:-1]
        occurs=int(line.split()[4])
        interacts=int(line.split()[5])
        plantID=int(line.split()[6])
        gallerID=int(line.split()[7])
        postmean=float(line.split()[8])
        postSD=float(line.split()[9])
        sppair=(plant,galler)
        pdict[sppair]=((postmean,postSD))        
      else:
        print len(line.split())
        print line.split()
  f.close()

  return pdict


def posterior_sampling(pdict):
  random_ints={}
  for pair in pdict:
    samples=np.random.binomial(1,pdict[pair][0],100)
    random_ints[pair]=samples

  return random_ints


def write_posterior_webs(random_ints,postdir):
  for i in range(1,101):
    outfile=open(postdir+'/P'+str(i)+'.web','w')
    for (plant, galler) in random_ints:
      if random_ints[(plant,galler)][i-1]==1:
        outfile.write(plant+'\t'+galler+'\n')
    outfile.close()

def filter_posterior_webs(postdir,proportion,altdir):
  for postweb in os.listdir(postdir):
    linklist=[]
    f=open(postdir+postweb,'r')
    for line in f:
      linklist.append('&'.join([line.split()[0],line.split()[1]]))
    f.close()
    for j in range(1,101):
      newcount=int(proportion*len(linklist))
      newlist=np.random.choice(linklist,newcount,replace=False)
      g=open(altdir+str(proportion)+'/'+postweb.split('.web')[0]+'_'+str(j)+'.web','w')
      for item in newlist:
        g.write(item.split('&')[0]+'\t'+item.split('&')[1]+'\n')
      g.close()


def main():

  for nettype in ['SG','GP']:
    if nettype=='SG':
      datafile='../data/Salix_example/Salix_Galler/posterior_probabilities.tsv'
    else:
      datafile='../data/Salix_example/Galler_Parasitoid/posterior_probabilities.tsv'
    if nettype=="SG":
      postdir='../data/randomised_webs/posterior/'
      altdir='../data/randomised_webs/detection_filter/'
    else:
      postdir='../data/randomised_webs/gp_posterior/'
      altdir='../data/randomised_webs/gp_detection_filter/'

    pdict=read_in_data(datafile) # Get posterior dist for each interaction 
    random_ints=posterior_sampling(pdict) # Calculate a set of random trails for each int
    write_posterior_webs(random_ints,postdir) # Create a set of random posterior webs

    for proportion in [0.5,0.6,0.7,0.8,0.9,0.95,0.99]:
      print nettype, proportion
    # proportion=0.8 # Eventually I probably want to do this for different proportions of links detected
      filter_posterior_webs(postdir,proportion,altdir)


  
if __name__ == '__main__':
  main()

