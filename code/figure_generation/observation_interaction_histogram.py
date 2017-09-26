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


def format_graph(graph,graphtype):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1
  graph.world.xmin=-0.5
  graph.world.xmax=70
  graph.world.ymin=0.5
  graph.world.ymax=1001

  graph.yaxis.tick.configure(onoff='on',minor_ticks=1,major_size=.5,minor_size=.3,place='normal',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=0)
  graph.yaxis.scale='Logarithmic'

  graph.xaxis.tick.configure(major=10,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
  graph.xaxis.ticklabel.configure(char_size=.75,format='decimal',prec=0)

  if graphtype=='occurs':
    xtext="Observed co-occurances"
  else:
    xtext="Observed interactions"
  graph.xaxis.label.configure(text=xtext,char_size=1,just=2,place='normal')
  graph.yaxis.label.configure(text="Species pairs",char_size=1,just=2)
  graph.legend.configure(box_linestyle=0,fill=0,fill_pattern=0,char_size=.75,
    loc=(100,.75),loctype='world')
  # graph.add_drawing_object(DrawText,text="Threshold",x=150,y=.9,char_size=.75,just=2,loctype='world')
  graph.panel_label.configure(loc='iur',char_size=.75,dx=.02,dy=.02)
  return graph

def populate_graph(graph,occurs,graphtype):
  if graphtype=='occurs':
    col=11
    star=graph.add_dataset([(0.1,860)],type='xy')
    star.symbol.configure(shape=4,size=.65)

  else:
    col=2
  occ=graph.add_dataset(occurs,type='bar')
  occ.symbol.configure(fill_color=col,color=1,size=.5)
  occ.line.linestyle=0

  return graph

grace=MultiPanelGrace(colors=colors)

occurs=read_Rfiles('../../data/Salix_example/Salix_Galler/nonzero_cooccurances.tsv')
interacts=read_Rfiles('../../data/Salix_example/Salix_Galler/observed_interactions.tsv')

for graphtype in ['occurs','interacts']:
  if graphtype=='occurs':
    dataset=occurs
  else:
    dataset=interacts
  graph=grace.add_graph(Panel)
  graph=format_graph(graph,graphtype)
  graph=populate_graph(graph,dataset,graphtype)
# graph.set_view(0.15,0.15,0.95,0.65)

grace.multi(rows=2,cols=1,vgap=.07)
# graph.set_view(0.15,0.15,0.95,0.65)
grace.hide_redundant_labels()
grace.set_col_yaxislabel(rowspan=(None,None),col=0,label="Species pairs",just=2,char_size=1)
grace.graphs[0].set_view(0.15,0.585,0.95,0.95)
grace.graphs[1].set_view(0.15,0.15,0.95,0.515)
for graph in grace.graphs:
  print graph.get_view()
grace.write_file('../../manuscript/figures/Salix_Galler_histogram.eps')
