import sys
import os
import networkx as nx 
import numpy as np 
import scipy as sp 

def web_reader(webdir,proportion,post_number,sample_number): 
  gallers=set()
  salixes=set()
  links=set()
  if sample_number==0:
    f=open(webdir+'posterior/P'+str(post_number)+'.web','r')
  else:
    f=open(webdir+'detection_filter/'+str(proportion)+'/P'+str(post_number)+'_'+str(sample_number)+'.web','r')
  for line in f:
    galler=line.split()[0]
    salix=line.split()[1]
    gallers.add(galler)
    salixes.add(salix)
    links.add((galler,salix))
  f.close()

  # G=nx.Graph()
  # for sp in gallers:
  #   G.add_node(sp)
  # for sp in salixes:
  #   G.add_node(sp)
  # for (galler,salix) in links:
  #   G.add_edge(galler,salix)

  return gallers,salixes,links

def calculate_original_properties(gallers,salixes,links):
  # original_properties['L'][post_number]=len(links)
  # original_properties['S'][post_number]=len(G.nodes())
  # original_properties['LS'][post_number]=float(len(G.edges()))/float(len(salixes))
  # original_properties['LG'][post_number]=float(len(G.edges()))/float(len(gallers))
  C=float(len(links))/float(len(salixes)*len(gallers))

  return C

def write_summary_file(summary_dict):
  f=open('../data/randomised_webs/Connectance_table.tsv','w')
  f.write('Web\tProportion\tObs_C\tSample_mean\tSample_SD\n')
  for proportion in summary_dict:
    for post_number in summary_dict[proportion]:
      f.write(str(post_number)+'\t'+str(proportion)+'\t'+'\t'.join(summary_dict[proportion][post_number])+'\n')
  f.close()

def main():

  webdir='../data/randomised_webs/'

  # Starting with C
  summary_dict={0.5:{},0.6:{},0.7:{},0.8:{},0.9:{},0.95:{},0.99:{}}

  for post_number in range(1,101):
    gallers,salixes,links=web_reader(webdir,0,post_number,0)
    obs=calculate_original_properties(gallers,salixes,links)
    for proportion in [0.5,0.6,0.7,0.8,0.9,0.95,0.99]:
    # for proportion in [0.9,0.95,0.99]:
      filtered_C=[]
      for sample_number in range(1,101):
        gallers,salixes,links=web_reader(webdir,proportion,post_number,sample_number)
        C=calculate_original_properties(gallers,salixes,links)
        filtered_C.append(C)
      sample_mean=np.mean(filtered_C)
      sample_SD=np.std(filtered_C)

      summary_dict[proportion][post_number]=(str(obs),str(sample_mean),str(sample_SD))

  print "Ready to write out the connectance file"
  write_summary_file(summary_dict)

  
if __name__ == '__main__':
  main()
