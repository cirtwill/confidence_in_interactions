# ######## R code to generate data files:
# % Getting cdfs:
# pars=calculate_parameters(sg_int_probs,0,0)
# n=seq(0,2000,1)

# CDFs=matrix(nrow=3,ncol=2002)
# colnames(CDFs)=c("Threshold",n)
# r=1
# for(threshold in c(0.1,0.05,0.01)){
#     cdf=pbeta(threshold,shape1=pars[[1]],shape2=pars[[2]]+n)
#     CDFs[r,1]=threshold
#     CDFs[r,2:2002]=cdf
#     print(r)
#     r=r+1
# }


# samples=matrix(nrow=9,ncol=3)
# colnames(samples)=c("Threshold","Confidence","Samples")
# r=1
# for(threshold in c(0.1,0.05,0.01)){
#   for(confidence in c(0.9,0.95,0.975)){
#     n_obs=samples_for_threshold(threshold,confidence,pars)
#     samples[r,1]=threshold
#     samples[r,2]=confidence
#     samples[r,3]=n_obs
#     print(r)
#     r=r+1
#   }
# }

# write.table(samples,file='../../data/Salix_example/Zillis/Salix_Galler/samples_for_threshold.tsv',sep='\t')
# write.table(CDFs,file='../../data/Salix_example/Zillis/Salix_Galler/samplefigure.tsv',sep='\t')
# # write.table(samples,file='../data/Salix_example/Salix_Galler/samples_for_threshold.tsv',sep='\t')
# # write.table(CDFs,file='../data/Salix_example/Salix_Galler/samplefigure.tsv',sep='\t')


# # Galler-parasitoid stuff
# ######### R code to generate data files:
# # % Getting cdfs:
# gp_pars=calculate_parameters(gp_int_probs,0,0)
# n=seq(0,500,1)

# gp_CDFs=matrix(nrow=3,ncol=502)
# colnames(gp_CDFs)=c("Threshold",n)
# r=1
# for(threshold in c(0.1,0.05,0.01)){
#     cdf=pbeta(threshold,shape1=gp_pars[[1]],shape2=gp_pars[[2]]+n)
#     gp_CDFs[r,1]=threshold
#     gp_CDFs[r,2:502]=cdf
#     print(r)
#     r=r+1
#   }



# gp_samples=matrix(nrow=9,ncol=3)
# colnames(gp_samples)=c("Threshold","Confidence","Samples")
# r=1
# for(threshold in c(0.1,0.05,0.01)){
#   for(confidence in c(0.9,0.95,0.975)){
#     n_obs=samples_for_threshold(threshold,confidence,gp_pars)
#     gp_samples[r,1]=threshold
#     gp_samples[r,2]=confidence
#     gp_samples[r,3]=n_obs
#     print(r)
#     r=r+1
#   }
# }

# write.table(gp_samples,file='../../data/Salix_example/Zillis/Galler_Parasitoid/samples_for_threshold.tsv',sep='\t')
# write.table(gp_CDFs,file='../../data/Salix_example/Zillis/Galler_Parasitoid/samplefigure.tsv',sep='\t')
# # write.table(gp_samples,file='../data/Salix_example/Galler_Parasitoid/samples_for_threshold.tsv',sep='\t')
# # write.table(gp_CDFs,file='../data/Salix_example/Galler_Parasitoid/samplefigure.tsv',sep='\t')

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

def read_Rfiles(filename,nettype):
  datadict={}
  sampledict={}

  f=open(filename,'r')
  for line in f:
    if line.split()[0]!='"Threshold"':
      threshold=float(line.split()[1])
      ys=line.split()[2:]
      datadict[threshold]=ys
  f.close()

  if nettype=='SG':
    g=open('../../data/Salix_example/Zillis/Salix_Galler/samples_for_threshold.tsv','r')
    # g=open('../../data/Salix_example/Salix_Galler/samples_for_threshold.tsv','r')
  else:
    g=open('../../data/Salix_example/Zillis/Galler_Parasitoid/samples_for_threshold.tsv','r')    
    # g=open('../../data/Salix_example/Galler_Parasitoid/samples_for_threshold.tsv','r')    
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
    for n in range(0,len(datadict[threshold])):
      pointlist.append((n,float(datadict[threshold][n])))
    pointdict[threshold]=pointlist
  return pointdict

def format_graph(graph,nettype):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1
  graph.world.xmin=0
  if nettype=='SG':
    graph.world.xmax=1300
    major=250
    graph.legend.configure(box_linestyle=0,fill=0,fill_pattern=0,char_size=.75,
      loc=(750,.65),loctype='world')
  else:
    graph.world.xmax=300
    major=50
    graph.legend.configure(box_linestyle=0,fill=0,fill_pattern=0,char_size=.75,
      loc=(155,.65),loctype='world')
  graph.world.ymin=-0
  graph.world.ymax=1
  # graph.world.xmax=750
  # major=250

  graph.yaxis.tick.configure(major=.20,onoff='on',minor_ticks=1,major_size=.5,minor_size=.3,place='normal',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=1)

  graph.xaxis.tick.configure(major=major,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
  graph.xaxis.ticklabel.configure(char_size=.75,format='decimal',prec=0)

  graph.xaxis.label.configure(text="Number of samples",char_size=1,just=2,place='normal')
  graph.yaxis.label.configure(text="Cumulative density",char_size=1,just=2)
    # loc=(400,.65),loctype='world')
  # graph.add_drawing_object(DrawText,text="Threshold",x=150,y=.9,char_size=.75,just=2,loctype='world')
  # graph.panel_label.configure(placement='iur',char_size=.75,dx=.02,dy=.03)

  return graph

def populate_graph(graph,pointdict,nettype):
  for threshold in sorted(pointdict):
    print threshold
    if threshold==.1:
      x=1
    elif threshold==.05:
      x=3
    elif threshold==.01:
      x=2

    print len(pointdict[threshold])
    data=graph.add_dataset(pointdict[threshold])
    data.symbol.shape=0
    data.line.configure(linewidth=2,linestyle=x)
    # if nettype=='SG':
    data.legend="Threshold="+str(threshold)

  bar95=graph.add_dataset([(0,0.95),(2000,.95)])
  bar95.symbol.shape=0
  bar95.line.configure(linewidth=1,linestyle=1,color=2)

  # if nettype=="SG":
    # graph.add_drawing_object(DrawText,text="Salix-Galler",char_size=1,just=0,x=0,y=1.05,loctype='world')
  # else:
    # graph.add_drawing_object(DrawText,text="Galler-Parasitoid",char_size=1,just=2,x=40,y=1.05,loctype='world')
    # graph.add_drawing_object(DrawText,text="Galler-Parasitoid",char_size=1,just=0,x=0,y=1.05,loctype='world')

  return graph

def add_samplelines(graph,sampledict,nettype):
  dots95=[]
  dots90=[]
  dots975=[]

  for confidence in sampledict:
    for threshold in sampledict[confidence]:
      liner=graph.add_dataset([(sampledict[confidence][threshold],0.1),(sampledict[confidence][threshold],0)])
      liner.symbol.shape=0
      if confidence==0.9:
        col=5
      elif confidence==0.95:
        col=2
      else:
        col=11
      liner.line.configure(linewidth=2.75,linestyle=1,color=col)

  for threshold in sampledict[0.9]:
    dots90.append((sampledict[0.9][threshold],0.9)) 
  for threshold in sampledict[0.95]:
    dots95.append((sampledict[0.95][threshold],0.95)) 
  for threshold in sampledict[0.975]:
    dots975.append((sampledict[0.975][threshold],0.975)) 

  do=graph.add_dataset(dots90)
  do.line.linestyle=0
  do.symbol.configure(shape=2,color=5,size=.75,fill_color=5)

  dot=graph.add_dataset(dots95)
  dot.line.linestyle=0
  dot.symbol.configure(shape=1,color=2,size=.75,fill_color=2)

  dots=graph.add_dataset(dots975)
  dots.line.linestyle=0
  dots.symbol.configure(shape=3,color=11,size=.75,fill_color=11)

  # if nettype=='SG':
  do.legend="P=0.90"
  dot.legend="P=0.95"
  dots.legend="P=0.975"

  print 'dots added'
  return graph



for nettype in ['SG','GP']:
  # grace=MultiPanelGrace(colors=colors)
  grace=Grace(colors=colors)
  if nettype=='SG':
    # datadict,sampledict=read_Rfiles('../../data/Salix_example/Salix_Galler/samplefigure.tsv',nettype)
    datadict,sampledict=read_Rfiles('../../data/Salix_example/Zillis/Salix_Galler/samplefigure.tsv',nettype)
  else:
    # datadict,sampledict=read_Rfiles('../../data/Salix_example/Galler_Parasitoid/samplefigure.tsv',nettype)    
    datadict,sampledict=read_Rfiles('../../data/Salix_example/Zillis/Galler_Parasitoid/samplefigure.tsv',nettype)    

  datasets=combiner(datadict)

  # graph=grace.add_graph(Panel)
  graph=grace.add_graph()
  graph=format_graph(graph,nettype)
  graph=populate_graph(graph,datasets,nettype)
  graph=add_samplelines(graph,sampledict,nettype)
  # graph.set_view(0.15,0.15,0.95,0.65)

# grace.automulti()
  # grace.multi(rows=2,cols=1,vgap=.08)
  # grace.hide_redundant_labels()
  # grace.set_col_yaxislabel(rowspan=(None,None),col=0,label="Cumulative density",just=2,char_size=1,perpendicular_offset=0.07)

  # grace.write_file('../../manuscript/figures/Salix_Galler_samples_and_cdfs.eps')
  grace.write_file('../../manuscript/figures/'+nettype+'_samples_and_cdfs_Zillis.eps')
