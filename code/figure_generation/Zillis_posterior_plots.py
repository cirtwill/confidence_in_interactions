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
colors.add_color(130,20,20,'red')
# colors.add_color(120,120,120,'grey')

obsvals={'SG':{'C':0.028,'LS':2.71,'LG':1.47,'NODF':0.5600535},
'GP':{'C':0.078,'LS':9.99,'LG':7.45,'NODF':6.853653}}

def read_summary_file(filename):
  whiskers={0.1:[],0.2:[],0.3:[],0.4:[],0.5:[],0.6:[],0.7:[],0.8:[],0.9:[]}
  obs={0.1:[],0.2:[],0.3:[],0.4:[],0.5:[],0.6:[],0.7:[],0.8:[],0.9:[]}
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
      obs[proportion].append(obs_C)
      whiskers[proportion].append(sample_mean)
  f.close()
  return whiskers, obs

def read_NODF_file(filename):
  whiskers={0.1:[],0.2:[],0.3:[],0.4:[],0.5:[],0.6:[],0.7:[],0.8:[],0.9:[]}
  obs={0.1:[],0.2:[],0.3:[],0.4:[],0.5:[],0.6:[],0.7:[],0.8:[],0.9:[]}
  pointdict={}
  f=open(filename,'r')
  for line in f:
    if line.split()[0]!='"Web"':
      proportion=float(line.split()[2][1:-1]) # Proportion of links observed
      obs_C=float(line.split()[3][1:-1]) # C in the full posterior web
      sample_mean=float(line.split()[4][1:-1]) # Mean C in the detection filter webs
      sample_SD=float(line.split()[5][1:-1])
      sample_SE=sample_SD/10 # SD/sqrt(100 filter webs)
      # CI is +- 1.96SE
      if proportion not in pointdict:
        pointdict[proportion]=[(obs_C,sample_mean,sample_SD)]
      else:
        pointdict[proportion].append((obs_C,sample_mean,sample_SD))
      obs[proportion].append(obs_C)
      whiskers[proportion].append(sample_mean)
  f.close()
  return whiskers, obs

def mean_calculator(listdict):
  dots=[]
  for prop in listdict:
    mean=np.mean(listdict[prop])
    sd=np.std(listdict[prop])
    dots.append((prop,mean,sd))
  return dots

def obs_breaker(obsdict):
  maxdots=[]
  meandots=[]
  mindots=[]
  for prop in sorted(obsdict):
    meandots.append((prop,np.mean(obsdict[prop])))
    maxdots.append((prop,np.mean(obsdict[prop])+np.std(obsdict[prop])))
    if prop==0.5:
      mindots.append((0,np.mean(obsdict[prop])-np.std(obsdict[prop])))
    if prop==0.99:
      mindots.append((1,np.mean(obsdict[prop])-np.std(obsdict[prop])))
    mindots.append((prop,np.mean(obsdict[prop])-np.std(obsdict[prop])))
  return maxdots,meandots,mindots

def Zillis_yaxes(prop,graph):
  graph.world.ymin=0.0
  if prop=='C':
    graph.world.ymax=0.07
    graph.yaxis.tick.configure(major=0.02,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=2)
    ytext="Connectance"
  elif prop=='LS':
    graph.world.ymax=4
    graph.yaxis.tick.configure(major=1,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    ytext="Links per Salix"
  elif prop=='LG':
    graph.world.ymax=7
    graph.yaxis.tick.configure(major=2,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    ytext="Links per galler"
  else:
    graph.world.ymax=2.25
    graph.yaxis.tick.configure(major=.5,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=1)
    ytext='NODF'
  graph.yaxis.label.configure(text=ytext,char_size=1,just=2)

  return graph

def Zillertal_yaxes(prop,graph):
  if prop=='C':
    graph.world.ymin=0
    graph.world.ymax=0.12
    graph.yaxis.tick.configure(major=0.03,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=2)
    ytext="Connectance"
  elif prop=='LS':
    graph.world.ymin=0
    graph.world.ymax=6
    graph.yaxis.tick.configure(major=2,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    ytext="Links per Salix"
  elif prop=='LG':
    graph.world.ymin=0
    graph.world.ymax=11
    graph.yaxis.tick.configure(major=2,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    ytext="Links per galler"
  else:
    graph.world.ymin=0
    graph.world.ymax=4
    graph.yaxis.tick.configure(major=1,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
    ytext='NODF'
  graph.yaxis.label.configure(text=ytext,char_size=1,just=2)

  return graph

def format_graph(graph,prop,nettype,site):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1
  # Xaxis will be strength of filter
  graph.world.xmin=0
  graph.world.xmax=1.0
  graph.world.ymin=0.0
  graph.xaxis.tick.configure(major=0.2,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
  graph.xaxis.ticklabel.configure(char_size=.75,format='decimal',prec=1)
  graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=0)
  if nettype=='SG' and site=='Zillis':
    graph=Zillis_yaxes(prop,graph)
  elif nettype=='SG' and site=='Zillertal':
    graph=Zillertal_yaxes(prop,graph)
  elif nettype=='GP':
    if prop=='C':
      graph.world.ymax=0.2
      graph.yaxis.tick.configure(major=0.05,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
      graph.yaxis.ticklabel.configure(char_size=.75,format='decimal',prec=2)
      ytext="Connectance"
    elif prop=='LS':
      graph.world.ymax=18
      graph.yaxis.tick.configure(major=4,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
      ytext="Links per galler"
    elif prop=='LG':
      graph.world.ymax=23    
      graph.yaxis.tick.configure(major=5,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
      ytext="Links per natural enemy"
    else:
      graph.world.ymax=8     
      graph.yaxis.tick.configure(major=2,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
      ytext="NODF"
    graph.yaxis.label.configure(text=ytext,char_size=1,just=2)
  if nettype=='GP':
    graph.legend.configure(box_linestyle=0,fill=0,fill_pattern=0,char_size=.5,
      loc=(0.527,3),loctype='world')
  else:
    graph.legend.configure(box_linestyle=0,fill=0,fill_pattern=0,char_size=.5,
      loc=(0.527,1.75),loctype='world')
  graph.panel_label.configure(placement='iul',char_size=.75,dx=.02,dy=.01)
  # if prop=='C' and nettype=='SG':
  #   graph.xaxis.label.configure(text='Salix-Galler',place='opposite',char_size=1.25,just=2,loctype='world')
  # elif prop=='C' and nettype=='GP':
    # graph.xaxis.label.configure(text="Galler-Parasitoid",place='opposite',char_size=1.25,just=2,loctype='world')

  return graph

def populate_graph(graph,prop,nettype,obs,whiskers):
  obsval=obsvals[nettype][prop]
  obsline=graph.add_dataset([(-1,obsval),(1,obsval)])
  obsline.symbol.shape=0
  obsline.line.configure(linestyle=2,linewidth=2,color='red')

  maxdots,meandots,mindots=obs_breaker(obs)
  maxline=graph.add_dataset(maxdots,type='xy')
  maxline.symbol.shape=0
  maxline.line.configure(linestyle=0)
  maxline.fill.configure(type=2,pattern=1,color=5)
  minline=graph.add_dataset(mindots,type='xy')
  minline.symbol.shape=0
  minline.line.configure(linestyle=0)
  minline.fill.configure(type=2,pattern=1,color=0)
  meandots=graph.add_dataset(meandots,type='xy')
  meandots.line.linestyle=0
  meandots.symbol.configure(shape=1,size=.5,color=10,fill_color=10)

  whiskerdots=mean_calculator(whiskers)
  dots=graph.add_dataset(whiskerdots,type='xydy')
  dots.line.linestyle=0
  dots.symbol.configure(shape=1,size=.5,fill_color=6)

  obsline2=graph.add_dataset([(-1,obsval),(1,obsval)])
  obsline2.symbol.shape=0
  obsline2.line.configure(linestyle=2,linewidth=2,color='red')

  if prop=="NODF" and nettype=="GP":
    obsline.legend="Observed"
    meandots.legend="Posterior mean"
    dots.legend="Filtered mean"
  if prop=="LS" and nettype=="SG":
    obsline.legend="Observed"
    meandots.legend="Posterior mean"
    dots.legend="Filtered mean"

  return graph

def main():

  for site in ['Zillis']:
  # for site in ['Zillis','Zillertal']:
    for nettype in ['SG','GP']:
      grace=MultiPanelGrace(colors=colors)
      for prop in ['C','LS','NODF','LG']:
        graph=grace.add_graph(Panel)
        graph=format_graph(graph,prop,nettype,site)

        if prop=='C':
          whiskers,obs=read_summary_file('../../data/Zillis_webs/'+nettype+'_Connectance_table_'+site+'.tsv')
        elif prop=='LS':
          whiskers,obs=read_summary_file('../../data/Zillis_webs/'+nettype+'_LS_table_'+site+'.tsv') # Salix or galler
        elif prop=='LG':
          whiskers,obs=read_summary_file('../../data/Zillis_webs/'+nettype+'_LG_table_'+site+'.tsv') # Galler or parasitoid
        else:
          whiskers,obs=read_NODF_file('../../data/Zillis_webs/'+nettype+'_NODF_table_'+site+'.tsv') # Galler or parasitoid

        graph=populate_graph(graph,prop,nettype,obs,whiskers)


      grace.multi(rows=2,cols=2,vgap=.04,hgap=.1)
      grace.hide_redundant_labels()

      for graph in grace.graphs:
        print graph.get_view()

      grace.graphs[0].set_view(0.15, 0.7020588235294118, 0.4970588235294118, 0.95)
      grace.graphs[1].set_view(0.5970588235294119, 0.7020588235294118, 0.95, 0.95)
      grace.graphs[2].set_view(0.15, 0.4141176470588234, 0.4970588235294118, 0.6620588235294117)
      grace.graphs[3].set_view(0.5970588235294119, 0.4141176470588234, 0.95, 0.6620588235294117)
      grace.set_row_xaxislabel(colspan=(None,None),row=1,label="Percent of links detected",just=2,char_size=1.5,perpendicular_offset=0.075)

      grace.write_file('../../manuscript/figures/'+nettype+'_posterior_properties_'+site+'.eps')

  
if __name__ == '__main__':
  main()
