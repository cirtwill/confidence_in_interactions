import os
import sys
import math
import random
from decimal import *
import numpy as np
import scipy as sp

#   numpy.random.binomial(n,p,size) draws [size] samples from a binomial distribution
#   % Each observation is a Bernoulli, or a binomial with n=1. Use p=mean(posterior).

def read_in_data(datafile):
  pdict={}
  f=open(datafile,'r')
  for line in f:
    if line.split()[0]!='"Salix"':
      if len(line.split())==7:
        plant=line.split()[1][1:-1]
        galler=line.split()[2][1:-1]
        occurs=int(line.split()[3])
        interacts=int(line.split()[4])
        postmean=float(line.split()[5])
        postSD=float(line.split()[6])
        sppair=(plant,galler)
        pdict[sppair]=((postmean,postSD))
      elif len(line.split())==10:
        plant='.'.join([line.split()[1][1:],line.split()[2],line.split()[3],'?']) # 4 is always '?'
        galler=line.split()[5][1:-1]
        occurs=int(line.split()[6])
        interacts=int(line.split()[7])
        postmean=float(line.split()[8])
        postSD=float(line.split()[9])
        sppair=(plant,galler)
        pdict[sppair]=((postmean,postSD))
      elif len(line.split())==9:
        plant='.'.join([line.split()[1][1:],line.split()[2],line.split()[3][:-1]]) 
        galler=line.split()[4][1:-1]
        occurs=int(line.split()[5])
        interacts=int(line.split()[6])
        postmean=float(line.split()[7])
        postSD=float(line.split()[8])
        sppair=(plant,galler)
        pdict[sppair]=((postmean,postSD))
      elif len(line.split())==8:
        plant='.'.join([line.split()[1][1:],line.split()[2][:-1]]) 
        galler=line.split()[3][1:-1]
        occurs=int(line.split()[4])
        interacts=int(line.split()[5])
        postmean=float(line.split()[6])
        postSD=float(line.split()[7])
        sppair=(plant,galler)
        pdict[sppair]=((postmean,postSD))        
      else:
        print len(line.split())
  f.close()

  return pdict


def posterior_sampling(pdict):
  random_ints={}
  for pair in pdict:
    samples=np.random.binomial(1,pdict[pair][0],10)
    random_ints[pair]=samples

  return random_ints


def write_posterior_webs(random_ints):
  for i in range(1,11):
    outfile=open('../data/randomised_webs/gp_posterior/P'+str(i)+'.web','w')
    for (plant, galler) in random_ints:
      if random_ints[(plant,galler)][i-1]==1:
        outfile.write(plant+'\t'+galler+'\n')
    outfile.close()

def filter_posterior_webs(postdir,proportion):
  for postweb in os.listdir(postdir):
    print postweb
    linklist=[]
    f=open(postdir+postweb,'r')
    for line in f:
      linklist.append('&'.join([line.split()[0],line.split()[1]]))
    f.close()
    for j in range(1,11):
      newcount=int(proportion*len(linklist))
      newlist=np.random.choice(linklist,newcount,replace=False)
      g=open('../data/randomised_webs/gp_detection_filter/'+str(proportion)+'/'+postweb.split('.web')[0]+'_'+str(j)+'.web','w')
      for item in newlist:
        g.write(item.split('&')[0]+'\t'+item.split('&')[1]+'\n')
      g.close()


def main():

  datafile='../data/Salix_example/Galler_Parasitoid/posterior_probabilities.tsv'
  pdict=read_in_data(datafile) # Get posterior dist for each interaction 
  random_ints=posterior_sampling(pdict) # Calculate a set of random trails for each int
  write_posterior_webs(random_ints) # Create a set of random posterior webs

  proportion=0.8 # Eventually I probably want to do this for different proportions of links detected
  postdir='../data/randomised_webs/gp_posterior/'
  filter_posterior_webs(postdir,proportion)


  
if __name__ == '__main__':
  main()

