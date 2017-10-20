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

def format_graph(graph,prop):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  if prop=='C':
    graph.world.xmin=0.5
    graph.world.ymin=0.45
    graph.world.xmax=0.55
    graph.world.ymax=0.55
    graph.yaxis.tick.configure(major=0.01,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    graph.xaxis.tick.configure(major=0.01,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    xtext="Posterior connectance"
    ytext="Detected connectance"
    graph.xaxis.ticklabel.configure(char_size=.75,format='decimal',prec=2)
    graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=2)
  elif prop=='LS':
    graph.world.xmin=33
    graph.world.ymin=30
    graph.world.xmax=39
    graph.world.ymax=36.5
    graph.yaxis.tick.configure(major=2,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    graph.xaxis.tick.configure(major=1,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    xtext="Posterior L/S"
    ytext="Detected L/S"
    graph.xaxis.ticklabel.configure(char_size=.75,format='decimal',prec=0)
    graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=0)
  else:
    graph.world.xmin=0
    graph.world.ymin=0
    graph.world.xmax=1
    graph.world.ymax=1
    graph.yaxis.tick.configure(major=.20,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    graph.xaxis.tick.configure(major=.20,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    xtext="Posterior %Specialists"
    ytext="Detected %Specialists"

    graph.xaxis.ticklabel.configure(char_size=.75,format='decimal',prec=2)
    graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=1)

  graph.xaxis.label.configure(text=xtext,char_size=1,just=2,place='normal')
  graph.yaxis.label.configure(text=ytext,char_size=1,just=2)
  graph.legend.configure(box_linestyle=0,fill=0,fill_pattern=0,char_size=.75,
    loc=(0.55,.25),loctype='world')
  # graph.add_drawing_object(DrawText,text="Threshold",x=150,y=.9,char_size=.75,just=2,loctype='world')
  graph.panel_label.configure(loc='iur',char_size=.75,dx=.02,dy=.02)

  # if graphtype=='occurs':
  #   graph.add_drawing_object(DrawText,text="Salix-Galler",x=35,y=1500,char_size=1.5,just=2,loctype='world')

  return graph

def populate_graph(graph,pointdict,prop):
  col=10
  for proportion in sorted(pointdict,reverse=True):
    dots=graph.add_dataset(pointdict[proportion],type='xydy')
    dots.line.linestyle=0
    # if prop=='C':
    #   dots.legend=str(proportion)+'% of links detected'
    dots.symbol.configure(shape=1,size=.5,fill_color=col)
    col-=1

  # if prop=='C':
  #   print pointdict[0.5]

  return graph

def main():

  grace=MultiPanelGrace(colors=colors)
  for prop in ['C','LS']:
    graph=grace.add_graph(Panel)
    graph=format_graph(graph,prop)

    if prop=='C':
      pointdict=read_summary_file('../../data/randomised_webs/Connectance_table.tsv')
    elif prop=='LS':
      pointdict=read_summary_file('../../data/randomised_webs/LS_table.tsv')
    else:
      pointdict=read_summary_file('../../data/randomised_webs/Specialists_table.tsv')

    graph=populate_graph(graph,pointdict,prop)

    for graph in grace.graphs:
      print graph.get_view()

  grace.multi(rows=2,cols=1,vgap=.07)
    # grace.graphs[0].set_view(0.15,0.15,0.95,0.85)

  grace.write_file('../../manuscript/figures/Salix_Galler_posterior_properties.eps')

  
if __name__ == '__main__':
  main()
