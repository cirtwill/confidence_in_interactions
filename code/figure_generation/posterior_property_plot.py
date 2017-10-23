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

colors=ColorBrewerScheme('Blues')  # The blue is very beautiful but maybe harder to see.
# colors.add_color(120,120,120,'grey')

def read_summary_file(filename):
  pointdict={}
  f=open(filename,'r')
  for line in f:
    if line.split()[0]!='Web':
      proportion=float(line.split()[1]) # Proportion of links observed
      obs_C=float(line.split()[2]) # C in the full posterior web
      sample_mean=float(line.split()[3]) # Mean C in the detection filter webs
      sample_SD=float(line.split()[4])
      sample_SE=sample_SD/10 # SD/sqrt(100 filter webs)
      # CI is +- 1.96SE
      if proportion not in pointdict:
        pointdict[proportion]=[(obs_C,sample_mean,sample_SD)]
      else:
        pointdict[proportion].append((obs_C,sample_mean,sample_SD))
  f.close()
  return pointdict

def format_graph(graph,prop,nettype):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  if prop=='C':
    if nettype=='SG':
      graph.world.xmin=0.52
      graph.world.ymin=0.20
      graph.world.xmax=0.575
      graph.world.ymax=0.60
    else:
      graph.world.xmin=0.18
      graph.world.ymin=0.075
      graph.world.xmax=0.20
      graph.world.ymax=0.20
    graph.yaxis.tick.configure(major=0.05,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    graph.xaxis.tick.configure(major=0.01,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    xtext="Posterior Connectance (C)"
    ytext="Detected C"
    graph.xaxis.ticklabel.configure(char_size=.75,format='decimal',prec=2)
    graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=2)
  elif prop=='LS':
    if nettype=='SG':
      graph.world.xmin=27
      graph.world.ymin=12
      graph.world.xmax=30
      graph.world.ymax=31
    else:
      graph.world.xmin=17
      graph.world.ymin=8
      graph.world.xmax=19
      graph.world.ymax=20

    graph.yaxis.tick.configure(major=2,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    graph.xaxis.tick.configure(major=1,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    xtext="Posterior Links per Resource (L/R)"
    ytext="Detected L/R"
    graph.xaxis.ticklabel.configure(char_size=.75,format='decimal',prec=0)
    graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=0)
  else:
    if nettype=='SG':
      graph.world.xmin=50
      graph.world.ymin=20
      graph.world.xmax=55
      graph.world.ymax=57.5
    else:
      graph.world.xmin=22.5
      graph.world.ymin=10
      graph.world.xmax=25
      graph.world.ymax=27.5      
    graph.yaxis.tick.configure(major=5,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    graph.xaxis.tick.configure(major=1,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    xtext="Posterior Links per Consumer (L/C)"
    ytext="Detected L/C"

    graph.xaxis.ticklabel.configure(char_size=.75,format='decimal',prec=0)
    graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=0)

  graph.xaxis.label.configure(text=xtext,char_size=1,just=2,place='normal')
  graph.yaxis.label.configure(text=ytext,char_size=1,just=2)
  graph.legend.configure(box_linestyle=0,fill=0,fill_pattern=0,char_size=.75,
    loc=(0.55,.25),loctype='world')
  # graph.add_drawing_object(DrawText,text="Threshold",x=150,y=.9,char_size=.75,just=2,loctype='world')
  graph.panel_label.configure(placement='iul',char_size=.75,dx=.02,dy=.01)

  # if graphtype=='occurs':
  #   graph.add_drawing_object(DrawText,text="Salix-Galler",x=35,y=1500,char_size=1.5,just=2,loctype='world')

  return graph

def populate_graph(graph,pointdict,prop,nettype):
  col=10
  for proportion in sorted(pointdict,reverse=True):
    dots=graph.add_dataset(pointdict[proportion],type='xydy')
    dots.line.linestyle=0
    # if prop=='C':
    #   dots.legend=str(proportion)+'% of links detected'
    dots.symbol.configure(shape=1,size=.5,fill_color=col)
    col-=1

  
  return graph

def main():

  grace=MultiPanelGrace(colors=colors)
  for prop in ['C','LS','LG']:
    for nettype in ['SG','GP']:
      graph=grace.add_graph(Panel)
      graph=format_graph(graph,prop,nettype)

      if prop=='C':
        pointdict=read_summary_file('../../data/randomised_webs/'+nettype+'_Connectance_table.tsv')
      elif prop=='LS':
        pointdict=read_summary_file('../../data/randomised_webs/'+nettype+'_LS_table.tsv') # Salix or galler
      else:
        pointdict=read_summary_file('../../data/randomised_webs/'+nettype+'_LG_table.tsv') # Galler or parasitoid

      graph=populate_graph(graph,pointdict,prop,nettype)

      # for graph in grace.graphs:
      # print graph.get_view()

  grace.multi(rows=3,cols=2,vgap=.08)
  grace.hide_redundant_labels()

  grace.set_row_xaxislabel(colspan=(None,None),row=0,label="Posterior Connectance (C)",just=2,char_size=1,perpendicular_offset=0.05)
  grace.set_row_xaxislabel(colspan=(None,None),row=1,label="Posterior Links per Resource (L/R)",just=2,char_size=1,perpendicular_offset=0.05)
  grace.set_row_xaxislabel(colspan=(None,None),row=2,label="Posterior Links per Consumer (L/C)",just=2,char_size=1,perpendicular_offset=0.05)


  grace.write_file('../../manuscript/figures/Salix_Galler_posterior_properties.eps')

  
if __name__ == '__main__':
  main()
