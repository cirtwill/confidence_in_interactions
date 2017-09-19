# ######### R code to generate data files:
# # % Getting cdfs:
# pars=calculate_parameters(int_probs,0,0)
# n=seq(0,500,1)

# CDFs=matrix(nrow=4,ncol=502)
# colnames(CDFs)=c("Threshold",n)
# r=1
# for(threshold in c(0.5,0.25,0.1,0.05)){
#   # for(confidence in c(0.9,0.95,0.975)){
#     cdf=pbeta(threshold,shape1=pars[[1]],shape2=pars[[2]]+n)
#     # samples=length(which(cdf<confidence))
#     CDFs[r,1]=threshold
#     # CDFs[r,2]=confidence
#     # CDFs[r,3]=samples
#     CDFs[r,2:502]=cdf
#     print(r)
#     r=r+1
#   }
# # }


# samples=matrix(nrow=12,ncol=3)
# colnames(samples)=c("Threshold","Confidence","Samples")
# r=1
# for(threshold in c(0.5,0.25,0.1,0.05)){
#   for(confidence in c(0.9,0.95,0.975)){
#     cdf=pbeta(threshold,shape1=pars[[1]],shape2=pars[[2]]+n)
#     n_obs=length(which(cdf<confidence))
#     samples[r,1]=threshold
#     samples[r,2]=confidence
#     samples[r,3]=n_obs
#     print(r)
#     r=r+1
#   }
# }

# write.table(samples,file='../data/Salix_example/Salix_Galler/samples_for_threshold.tsv',sep='\t')
# write.table(CDFs,file='../data/Salix_example/Salix_Galler/samplefigure.tsv',sep='\t')

import sys
import math
import random
from decimal import *
import numpy as np

#Pygrace libraries
from PyGrace.grace import Grace
from PyGrace.colors import RandomColorScheme, MarkovChainColorScheme, ColorBrewerScheme
from PyGrace.dataset import SYMBOLS
from PyGrace.Extensions.panel import NetworkPanel,Panel,MultiPanelGrace
from PyGrace.drawing_objects import DrawText, DrawLine, DrawBox
from PyGrace.axis import LINEAR_SCALE, LOGARITHMIC_SCALE
from PyGrace.Extensions.network import Network

from PyGrace.Extensions.distribution import CDFGraph, PDFGraph
from PyGrace.Extensions.latex_string import LatexString, CONVERT
from PyGrace.Extensions.colorbar import SolidRectangle, ColorBar
from PyGrace.Styles.el import ElGraph, ElLinColorBar, ElLogColorBar

colors=ColorBrewerScheme('Spectral')  # The blue is very beautiful but maybe harder to see.
# colors.add_color(120,120,120,'grey')

def read_Rfiles(filename):
  datadict={}
  sampledict={}

  f=open(filename,'r')
  for line in f:
    if line.split()[0]!='"Threshold"':
      threshold=float(line.split()[1])
      ys=line.split()[2:]
      datadict[threshold]=ys
  f.close()

  g=open('../../data/Salix_example/Salix_Galler/samples_for_threshold.tsv')
  for line in g:
    if line.split()[0]!='"Threshold"':
      threshold=float(line.split()[1])
      confidence=float(line.split()[2])
      samples=int(line.split()[3])
      if confidence in sampledict:
        sampledict[confidence][threshold]=samples
      else:
        sampledict[confidence]={threshold:samples}
  g.close()
  return datadict,sampledict

def combiner(datadict):
  pointdict={}
  for threshold in datadict:
    pointlist=[]
    for n in range(0,501):
      pointlist.append((n,float(datadict[threshold][n])))
    pointdict[threshold]=pointlist
  return pointdict

def format_graph(graph):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1
  graph.world.xmin=0
  graph.world.xmax=200
  graph.world.ymin=-0
  graph.world.ymax=1

  graph.yaxis.tick.configure(major=.20,onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='normal',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=1)

  graph.xaxis.tick.configure(major=50,onoff='on',minor_ticks=0,major_size=.5,place='normal',minor_size=.5,major_linewidth=1,minor_linewidth=1)
  graph.xaxis.ticklabel.configure(char_size=.75,format='decimal',prec=0)

  graph.xaxis.label.configure(text="Number of samples",char_size=1,just=2,place='normal')
  graph.yaxis.label.configure(text="Cumulative density",char_size=1,just=2)
  graph.legend.configure(box_linestyle=0,fill=0,fill_pattern=0,char_size=.75,
    loc=(125,.75),loctype='world')
  # graph.add_drawing_object(DrawText,text="Threshold",x=150,y=.9,char_size=.75,just=2,loctype='world')

  return graph

def populate_graph(graph,pointdict):
  for threshold in pointdict:
    if threshold==.5:
      x=1
    elif threshold==.25:
      x=3
    elif threshold==.1:
      x=5
    else:
      x=2

    data=graph.add_dataset(pointdict[threshold])
    data.symbol.shape=0
    data.line.configure(linewidth=2,linestyle=x)
    data.legend="Threshold="+str(threshold)

  bar95=graph.add_dataset([(0,0.95),(1000,.95)])
  bar95.symbol.shape=0
  bar95.line.configure(linewidth=1,linestyle=1,color=2)

  return graph

def add_samplelines(graph,sampledict):
  dots95=[]
  dots90=[]
  dots975=[]

  for confidence in sampledict:
    for threshold in sampledict[confidence]:
      liner=graph.add_dataset([(sampledict[confidence][threshold],0.1),(sampledict[confidence][threshold],0)])
      liner.symbol.shape=0
      liner.line.configure(linewidth=.75,linestyle=1,color=1)

  for threshold in sampledict[0.9]:
    dots90.append((sampledict[0.9][threshold],0.9)) 
  for threshold in sampledict[0.95]:
    dots95.append((sampledict[0.95][threshold],0.95)) 
  for threshold in sampledict[0.975]:
    dots975.append((sampledict[0.975][threshold],0.975)) 

  do=graph.add_dataset(dots90)
  do.line.linestyle=0
  do.symbol.configure(shape=2,color=5,size=.75,fill_color=5)
  do.legend="P=0.90"

  dot=graph.add_dataset(dots95)
  dot.line.linestyle=0
  dot.symbol.configure(shape=1,color=2,size=.75,fill_color=2)
  dot.legend="P=0.95"

  dots=graph.add_dataset(dots975)
  dots.line.linestyle=0
  dots.symbol.configure(shape=3,color=11,size=.75,fill_color=11)
  dots.legend="P=0.975"

  print 'sqif'
  return graph

grace=Grace(colors=colors)

datadict,sampledict=read_Rfiles('../../data/Salix_example/Salix_Galler/samplefigure.tsv')

datasets=combiner(datadict)

graph=grace.add_graph()
graph=format_graph(graph)
graph=populate_graph(graph,datasets)
graph=add_samplelines(graph,sampledict)
graph.set_view(0.15,0.15,0.95,0.65)

grace.write_file('../../manuscript/figures/Salix_Galler_samples_and_cdfs.eps')
