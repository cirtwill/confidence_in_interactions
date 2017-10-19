# ####### R code to generate data files:
# # Getting distributions:
# xdata=matrix(ncol=1001,nrow=12)
# ydata=matrix(ncol=1001,nrow=12)
# ns=as.vector(c(0,5,10,15,20,25,50,100,150,200,300,374))
# for(i in 1:length(ns)){
#   n=ns[i]
#   dist=calculate_distribution(calculate_parameters(sg_int_probs,n,0))
#   x=seq(-5,5,length=1000)*sqrt(dist[[2]])+dist[[1]]
#   hx=dnorm(x,dist[[1]],sqrt(dist[[2]]))
#   xdata[i,]=c(n,x)
#   ydata[i,]=c(n,hx)
# }

# # Data are along rows for pythonic convenience
# colnames(xdata)=c("N",seq(1,1000))
# colnames(ydata)=c("N",seq(1,1000))


# write.table(xdata,file='../../data/Salix_example/Salix_Galler/distfigure_xvals.tsv',sep='\t')
# write.table(ydata,file='../../data/Salix_example/Salix_Galler/distfigure_yvals.tsv',sep='\t')

# xdata=matrix(ncol=1001,nrow=12)
# ydata=matrix(ncol=1001,nrow=12)
# ns=as.vector(c(0,5,10,15,20,25,50,100,150,200,300,374))
# for(i in 1:length(ns)){
#   n=ns[i]
#   dist=calculate_distribution(calculate_parameters(gp_int_probs,n,0))
#   x=seq(-5,5,length=1000)*sqrt(dist[[2]])+dist[[1]]
#   hx=dnorm(x,dist[[1]],sqrt(dist[[2]]))
#   xdata[i,]=c(n,x)
#   ydata[i,]=c(n,hx)
# }

# # Data are along rows for pythonic convenience
# colnames(xdata)=c("N",seq(1,1000))
# colnames(ydata)=c("N",seq(1,1000))


# write.table(xdata,file='../../data/Salix_example/Galler_Parasitoid/distfigure_xvals.tsv',sep='\t')
# write.table(ydata,file='../../data/Salix_example/Galler_Parasitoid/distfigure_yvals.tsv',sep='\t')


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

def format_graph(graph,nettype):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1
  graph.world.xmin=0
  graph.world.xmax=.6
  graph.world.ymin=-0
  graph.world.ymax=40

  graph.yaxis.tick.configure(major=50,onoff='off',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.ticklabel.configure(char_size=0,format='decimal',prec=0)

  graph.xaxis.tick.configure(major=.1,onoff='on',minor_ticks=0,major_size=.5,place='normal',minor_size=.5,major_linewidth=1,minor_linewidth=1)
  graph.xaxis.ticklabel.configure(char_size=.75,format='decimal',prec=1)

  graph.xaxis.label.configure(text="Probability of interaction",char_size=1,just=2,place='normal')
  graph.yaxis.label.configure(text="Probability density",char_size=1,just=2)
  graph.legend.configure(box_linestyle=0,fill=0,fill_pattern=0,char_size=.75,
    loc=(0.475,37),loctype='world')
  if nettype=='SG':
    graph.add_drawing_object(DrawText,text="N",x=0.545,y=37.5,char_size=.75,just=2,loctype='world')
  graph.panel_label.configure(placement='iul',char_size=.75,dx=.02,dy=.01)


  return graph

def populate_graph(graph,pointdict,nettype):
  cols=[1,2,5,9,11,12]
  x=0
  for N in [0,5,10,20,50,100]:
    data=graph.add_dataset(pointdict[N])
    data.symbol.shape=0
    data.line.configure(linewidth=1,linestyle=1,color=cols[x])
    data.fill.configure(pattern=1,color=cols[x],type=1)
    if nettype=='SG': 
      data.legend=str(N)
    x+=1

  return graph


grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['A: Salix-Galler','B: Galler-Parasitoid'])
grace.set_label_scheme('dummy')

for nettype in ['SG','GP']:
  if nettype=='SG':
    xdata=read_Rfiles('../../data/Salix_example/Salix_Galler/distfigure_xvals.tsv')
    ydata=read_Rfiles('../../data/Salix_example/Salix_Galler/distfigure_yvals.tsv')
  else:
    xdata=read_Rfiles('../../data/Salix_example/Galler_Parasitoid/distfigure_xvals.tsv')
    ydata=read_Rfiles('../../data/Salix_example/Galler_Parasitoid/distfigure_yvals.tsv')

  datasets=combiner(xdata,ydata)

  graph=grace.add_graph(Panel)
  graph=format_graph(graph,nettype)
  graph=populate_graph(graph,datasets,nettype)
  graph.set_view(0.15,0.15,0.95,0.65)

grace.multi(rows=2,cols=1,vgap=.04)
grace.hide_redundant_labels()
grace.set_col_yaxislabel(rowspan=(None,None),col=0,label="Probability density",just=2,char_size=1,perpendicular_offset=.04)

grace.write_file('../../manuscript/figures/Salix_Galler_pdfs_increasing_N.eps')
