# # ######### R code to generate data files:
# galldata=read.csv('../../data/Salix_example/cooccur_interact_galler_salix.csv',header=TRUE,rownames=TRUE)
# n_pairs=nrow(galldata)
# unobserved=length(which(galldata$cooccur==0))
# observed=galldata[which(galldata$cooccur>0),]
# observed$int_freq=observed$interact/observed$cooccur
# obs_ints=length(which(observed$interact>0))

# observation_histdata=as.data.frame(tapply(observed$cooccur,observed$cooccur,length))
# interaction_histdata=as.data.frame(tapply(observed$interact,observed$interact,length))
# colnames(observation_histdata)=c("cooccurances")
# colnames(interaction_histdata)=c("interactions")

# write.table(observation_histdata,file='../../data/Salix_example/Salix_Galler/nonzero_cooccurances.tsv',sep='\t')
# write.table(interaction_histdata,file='../../data/Salix_example/Salix_Galler/observed_interactions.tsv',sep='\t')



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
  pointlist=[]

  f=open(filename,'r')
  for line in f:
    if line.split()[0] not in ['"cooccurances"','"interactions"']:
      x=int(line.split()[0][1:-1]) # Row names are number of observations.
      count=int(line.split()[1])
      pointlist.append((x,count))
  f.close()

  return pointlist

def read_csv(filename): # Just using the R input file here.
  dotlist=[]
  f=open(filename,'r')
  for line in f:
    newline=line.split('\r')[0].split(',')
    if newline[0]!='':
      cooccurances=int(newline[-2])
      interactions=int(newline[-1])
      dotlist.append((cooccurances,interactions))
  f.close()

  return dotlist

def format_graph(graph,graphtype):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  if graphtype!='dotplot':
    graph.world.xmin=-0.5
    graph.world.xmax=70
    graph.world.ymin=0.5
    graph.world.ymax=900

    graph.yaxis.tick.configure(onoff='on',minor_ticks=1,major_size=.5,minor_size=.3,place='normal',major_linewidth=1,minor_linewidth=1)
    graph.yaxis.scale='Logarithmic'
    graph.xaxis.tick.configure(major=10,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
  else:
    graph.world.xmin=-0.15
    graph.world.ymin=-0.15
    graph.world.xmax=70
    graph.world.ymax=70
    graph.yaxis.tick.configure(major=20,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    graph.xaxis.tick.configure(major=20,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)

  graph.xaxis.ticklabel.configure(char_size=.75,format='decimal',prec=0)
  graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=0)

  if graphtype=='interacts':
    xtext="Observed interactions"
  else:
    xtext="Observed co-occurances"
  if graphtype=='dotplot':
    ytext="Observed interactions"
  else:
    ytext='Species pairs'
  graph.xaxis.label.configure(text=xtext,char_size=1,just=2,place='normal')
  graph.yaxis.label.configure(text=ytext,char_size=1,just=2)
  graph.legend.configure(box_linestyle=0,fill=0,fill_pattern=0,char_size=.75,
    loc=(100,.75),loctype='world')
  # graph.add_drawing_object(DrawText,text="Threshold",x=150,y=.9,char_size=.75,just=2,loctype='world')
  graph.panel_label.configure(loc='iur',char_size=.75,dx=.02,dy=.02)

  if graphtype=='occurs':
    graph.add_drawing_object(DrawText,text="Salix-Galler",x=35,y=1500,char_size=1.5,just=2,loctype='world')

  return graph

def populate_graph(graph,occurs,graphtype):
  if graphtype=='occurs':
    col=11
  elif graphtype=='interacts':
    col=2
  if graphtype!='dotplot':
    occ=graph.add_dataset(occurs,type='bar')
    occ.symbol.configure(fill_color=col,color=1,size=.2)
    occ.line.linestyle=0
  else:
    dots=graph.add_dataset(occurs,type='xy')
    dots.symbol.configure(fill_color=1,color=1,size=.4)
    dots.line.linestyle=0

    lines=graph.add_dataset([(0,0),(1000,1000)])
    lines.symbol.shape=0
    lines.line.configure(color=2,linestyle=3,linewidth=2)

  return graph

grace=MultiPanelGrace(colors=colors)

occurs=read_Rfiles('../../data/Salix_example/Salix_Galler/nonzero_cooccurances.tsv')
interacts=read_Rfiles('../../data/Salix_example/Salix_Galler/observed_interactions.tsv')
dots=read_csv('../../data/Salix_example/cooccur_interact_galler_salix.csv')

for graphtype in ['occurs','interacts','dotplot']:
  if graphtype=='occurs':
    dataset=occurs
  elif dataset=='interacts':
    dataset=interacts
  else:
    dataset=dots
  graph=grace.add_graph(Panel)
  graph=format_graph(graph,graphtype)
  graph=populate_graph(graph,dataset,graphtype)
# graph.set_view(0.15,0.15,0.95,0.65)

grace.multi(rows=3,cols=1,vgap=.07,hgap=.04)
# graph.set_view(0.15,0.15,0.95,0.65)
grace.hide_redundant_labels()
# grace.set_col_yaxislabel(rowspan=(None,None),col=0,label="Species pairs",just=2,char_size=1)
# grace.graphs[0].set_view(0.15,0.73,0.55,0.95)
# grace.graphs[1].set_view(0.15,0.44,0.55,0.66)
# grace.graphs[2].set_view(0.15,0.15,0.55,0.37)
# for graph in grace.graphs:
#   print graph.get_view()
grace.write_file('../../manuscript/figures/Salix_Galler_histogram.eps')
