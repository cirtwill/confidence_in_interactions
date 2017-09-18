######### R code to generate data files:
# % Getting cdfs:
# pars=calculate_parameters(int_probs,0,0)
# n=seq(0,100,1)

# CDFs=matrix(nrow=9,ncol=104)
# colnames(CDFs)=c("Threshold","Confidence","Samples",n)
# r=1
# for(threshold in c(0.5,0.1,0.01)){
#   for(confidence in c(0.9,0.95,0.975)){
#     cdf=pbeta(threshold,shape1=pars[[1]],shape2=pars[[2]]+n)
#     samples=length(which(cdf<confidence))
#     CDFs[r,1]=threshold
#     CDFs[r,2]=confidence
#     CDFs[r,3]=samples
#     CDFs[r,4:104]=cdf
#     print(r)
#     r=r+1
#   }
# }

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
      confidence=float(line.split()[2])
      samples=int(line.split()[3])
      ys=line.split()[4:]
      datadict[threshold]=ys
      if threshold in sampledict:
        sampledict[threshold][confidence]=samples
      else:
        sampledict[threshold]={confidence:samples}
  f.close()
  print sampledict
  return datadict,sampledict

def combiner(datadict):
  pointdict={}
  for threshold in datadict:
    pointlist=[]
    for n in range(0,101):
      pointlist.append((float(datadict[threshold][n]),float(datadict[threshold][n])))
    pointdict[threshold]=pointlist
  return pointdict

def format_graph(graph):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1
  graph.world.xmin=0
  graph.world.xmax=.7
  graph.world.ymin=-0
  graph.world.ymax=500

  graph.yaxis.tick.configure(major=50,onoff='off',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.ticklabel.configure(char_size=0,format='decimal',prec=0)

  graph.xaxis.tick.configure(major=.1,onoff='on',minor_ticks=0,major_size=.5,place='normal',minor_size=.5,major_linewidth=1,minor_linewidth=1)
  graph.xaxis.ticklabel.configure(char_size=.75,format='decimal',prec=1)

  graph.xaxis.label.configure(text="Probability of interaction",char_size=1,just=2,place='normal')
  graph.yaxis.label.configure(text="Probability density",char_size=1,just=2)
  graph.legend.configure(box_linestyle=0,fill=0,fill_pattern=0,char_size=.75,
    loc=(0.6,450),loctype='world')
  graph.add_drawing_object(DrawText,text="N",x=0.655,y=450,char_size=.75,just=2,loctype='world')

  return graph

def populate_graph(graph,pointdict):
  x=1
  for threshold in pointdict:
    data=graph.add_dataset(pointdict[N])
    data.symbol.shape=0
    data.line.configure(linewidth=1,linestyle=1,color=x)
    data.fill.configure(pattern=1,color=x,type=1)
    data.legend=str(N)
    if N>10 and N<20:
      x+=2
    else:
      x+=1

  return graph


grace=Grace(colors=colors)

datadict,sampledict=read_Rfiles('../../data/Salix_example/Salix_Galler/samplefigure.tsv')

datasets=combiner(datadict)

graph=grace.add_graph()
graph=format_graph(graph)
graph=populate_graph(graph,datasets)
graph.set_view(0.15,0.15,0.95,0.65)

grace.write_file('../../manuscript/figures/Salix_Galler_pdfs_increasing_N.eps')
