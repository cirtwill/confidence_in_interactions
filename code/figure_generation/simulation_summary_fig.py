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

colors=ColorBrewerScheme('RdYlBu')  # The blue is very beautiful but maybe harder to see.
colors.add_color(130,20,20,'red')
# colors.add_color(120,120,120,'grey')

def read_property_lists(propfile):
  prop_dict={"N_plants":[],"N_gallers":[],"C":[],"L_per_plant":[],"L_per_galler":[],"NODF":[]}
  f=open(propfile,'r')
  for line in f:
    if line.split()[0] not in ['"N_plants"','"L_per_plant"']:
      if len(line.split())==7:
        prop_dict["N_plants"].append(float(line.split()[1]))
        prop_dict["N_gallers"].append(float(line.split()[2]))
        prop_dict["C"].append((float(line.split()[3])))
      else:
        prop_dict["N_plants"].append(15)
        prop_dict["N_gallers"].append(20)
        prop_dict["C"].append(0.5)
      prop_dict["L_per_plant"].append(float(line.split()[-3]))
      prop_dict["L_per_galler"].append(float(line.split()[-2]))
      prop_dict["NODF"].append(float(line.split()[-1]))

  mean_dict={}
  for prop in prop_dict:
    mean_dict[prop]=(np.mean(prop_dict[prop]),np.std(prop_dict[prop]))


  return mean_dict

def configure_yaxes(prop,graph):
  graph.yaxis.tick.configure(onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
  if prop=='N_plants':
    graph.world.ymax=15
    graph.yaxis.tick.configure(major=3)
    graph.yaxis.label.configure(text="N(plants)",char_size=1,just=2,place='normal')
  elif prop=='N_gallers':
    graph.world.ymax=20
    graph.yaxis.tick.configure(major=4)
    graph.yaxis.label.configure(text="N(gallers)",char_size=1,just=2,place='normal')
  elif prop=='C':
    graph.world.ymax=0.5
    graph.yaxis.tick.configure(major=0.1)
    graph.yaxis.ticklabel.configure(prec=1)
    graph.yaxis.label.configure(text="Connectance",char_size=1,just=2,place='normal')
  elif prop=='L_per_plant':
    graph.world.ymax=20
    graph.yaxis.tick.configure(major=4)
    graph.yaxis.label.configure(text="Links per plant",char_size=1,just=2,place='normal')
  elif prop=='L_per_galler':
    graph.world.ymax=15
    graph.yaxis.tick.configure(major=3)
    graph.yaxis.label.configure(text="Links per galler",char_size=1,just=2,place='normal')
  else:
    graph.world.ymax=40
    graph.yaxis.tick.configure(major=10)
    graph.yaxis.label.configure(text='NODF',char_size=1,just=2,place='normal')

  return graph

def format_graph(graph,prop,nettype,site):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=0
  graph.world.xmax=4
  graph.xaxis.tick.set_spec_ticks([1,2,3],[],tick_labels=['Original','Process uncertainty','Detection uncertainty'])
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.5,place='both',major_linewidth=1)
  graph.xaxis.ticklabel.configure(char_size=.75,angle=45)

  # Build y-axes in another function
  graph.world.ymin=0.0
  graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=0)
  graph=configure_yaxes(prop,graph)

  graph.panel_label.configure(placement='iul',char_size=.75,dx=.02,dy=.01)

  return graph

def populate_graph(graph,prop,A_props,B_props,C_props):
  original=graph.add_dataset([(1,A_props[prop][0],A_props[prop][1])],type='xydy')
  original.symbol.configure(shape=1,size=1,fill_color=2,color=2)

  process=graph.add_dataset([(2,B_props[prop][0],B_props[prop][1])],type='xydy')
  process.symbol.configure(shape=3,size=1,fill_color=6,color=1)

  detect=graph.add_dataset([(3,C_props[prop][0],C_props[prop][1])],type='xydy')
  detect.symbol.configure(shape=2,size=1,fill_color=10,color=10)

  return graph

def main():
  grace=MultiPanelGrace(colors=colors)
  grace.add_label_scheme('dummy',['A','B','C','D','E','F'])
  grace.set_label_scheme('dummy')
  A_props=read_property_lists('../../data/simulation_example/A_props.tsv')
  B_props=read_property_lists('../../data/simulation_example/B_props.tsv')
  C_props=read_property_lists('../../data/simulation_example/C_props.tsv')

  for prop in ["N_plants","N_gallers","C","L_per_plant","L_per_galler","NODF"]:
    graph=grace.add_graph(Panel)
    graph=format_graph(graph,prop,nettype,site)

    graph=populate_graph(graph,prop,A_props,B_props,C_props)

    grace.multi(rows=2,cols=3,vgap=.04,hgap=.1)
    grace.hide_redundant_labels()
    grace.set_row_xaxislabel(colspan=(None,None),row=1,label="Type of network",just=2,char_size=1.5,perpendicular_offset=0.05)

    grace.write_file('../../manuscript/figures/simulation_example.eps')

  
if __name__ == '__main__':
  main()
