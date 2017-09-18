######### R code to generate data files:
# % Getting distributions:
# xdata=matrix(ncol=1001,nrow=201)
# ydata=matrix(ncol=1001,nrow=201)
# for(n in 0:200){
#   dist=calculate_distribution(calculate_parameters(int_probs,2*n,0))
#   x=seq(-10,10,length=1000)*sqrt(dist[[2]])+dist[[1]]
#   hx=dnorm(x)
#   xdata[n+1,]=c(2*n,x)
#   ydata[n+1,]=c(2*n,hx)
# }

# % Data are along rows for pythonic convenience
# colnames(xdata)=c("N",seq(1,1000))
# colnames(ydata)=c("N",seq(1,1000))
# write.table(xdata,file='../data/Salix_example/Salix_Galler/distfigure_xvals.tsv',sep='\t')
# write.table(ydata,file='../data/Salix_example/Salix_Galler/distfigure_yvals.tsv',sep='\t')

# Something's funny about these data.

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

  f=open(filename,'r')
  for line in f:
    if line.split()[0]!='"N"':
      N=int(line.split()[1])
      vals=line.split()[2:]
      datadict[N]=vals
  f.close()
  return datadict

def combiner(xdata,ydata):
  pointdict={}
  for N in xdata:
    pointlist=[]
    for i in range(0,len(xdata[N])):
      pointlist.append((float(xdata[N][i]),float(ydata[N][i])))
    pointdict[N]=pointlist
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
  for N in [0,2,4,6,8,10,14,20,30,40,50]:
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

xdata=read_Rfiles('../../data/Salix_example/Salix_Galler/distfigure_xvals.tsv')
ydata=read_Rfiles('../../data/Salix_example/Salix_Galler/distfigure_yvals.tsv')

datasets=combiner(xdata,ydata)

graph=grace.add_graph()
graph=format_graph(graph)
graph=populate_graph(graph,datasets)
graph.set_view(0.15,0.15,0.95,0.65)

grace.write_file('../../manuscript/figures/Salix_Galler_pdfs_increasing_N.eps')
