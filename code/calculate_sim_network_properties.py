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

  G=nx.Graph()
  for sp in gallers:
    G.add_node(sp)
  for sp in salixes:
    G.add_node(sp)
  for (galler,salix) in links:
    G.add_edge(galler,salix)

  return gallers,salixes,links,G

def calculate_original_properties(gallers,salixes,links,G):
  # original_properties['L'][post_number]=len(links)
  # original_properties['S'][post_number]=len(G.nodes())
  # original_properties['LS'][post_number]=float(len(G.edges()))/float(len(salixes))
  # original_properties['LG'][post_number]=float(len(G.edges()))/float(len(gallers))
  C=float(len(links))/float(len(salixes)*len(gallers))
  LS=[]
  specialists=[]
  for node in G.nodes():
    links=len(G.neighbors(node))
    LS.append(links)
    if links==1:
      specialists.append(node)
  meanLS=np.mean(LS)
  Spec=float(len(specialists))/float(len(G.nodes()))

  return C,meanLS,Spec

def write_summary_file(summary_dict):
  f=open('../data/randomised_webs/Connectance_table.tsv','w')
  f.write('Web\tProportion\tObs_C\tSample_mean\tSample_SD\n')
  for proportion in summary_dict['C']:
    for post_number in summary_dict['C'][proportion]:
      f.write(str(post_number)+'\t'+str(proportion)+'\t'+'\t'.join(summary_dict['C'][proportion][post_number])+'\n')
  f.close()

  g=open('../data/randomised_webs/LS_table.tsv','w')
  g.write('Web\tProportion\tObs_LS\tSample_mean\tSample_SD\n')
  for proportion in summary_dict['LS']:
    for post_number in summary_dict['LS'][proportion]:
      g.write(str(post_number)+'\t'+str(proportion)+'\t'+'\t'.join(summary_dict['LS'][proportion][post_number])+'\n')
  g.close()

  h=open('../data/randomised_webs/Specialists_table.tsv','w')
  h.write('Web\tProportion\tObs_Percent_Spec\tSample_mean\tSample_SD\n')
  for proportion in summary_dict['Spec']:
    for post_number in summary_dict['Spec'][proportion]:
      h.write(str(post_number)+'\t'+str(proportion)+'\t'+'\t'.join(summary_dict['Spec'][proportion][post_number])+'\n')
  h.close()

def main():

  webdir='../data/randomised_webs/'

  # Starting with C
  summary_dict={'C':{0.5:{},0.6:{},0.7:{},0.8:{},0.9:{},0.95:{},0.99:{}},
  'LS':{0.5:{},0.6:{},0.7:{},0.8:{},0.9:{},0.95:{},0.99:{}},
  'Spec':{0.5:{},0.6:{},0.7:{},0.8:{},0.9:{},0.95:{},0.99:{}}  }

  for post_number in range(1,101):
    gallers,salixes,links,G=web_reader(webdir,0,post_number,0)
    obs=calculate_original_properties(gallers,salixes,links,G)
    for proportion in [0.5]:
    # for proportion in [0.5,0.6,0.7,0.8,0.9,0.95,0.99]:
    # for proportion in [0.9,0.95,0.99]:
      filtered_C=[]
      filtered_LS=[]
      filtered_Spec=[]
      for sample_number in range(1,101):
        gallers,salixes,links,G=web_reader(webdir,proportion,post_number,sample_number)
        C,LS,Spec=calculate_original_properties(gallers,salixes,links,G)
        filtered_C.append(C)
        filtered_LS.append(LS)
        filtered_Spec.append(Spec)
      C_sample_mean=np.mean(filtered_C)
      C_sample_SD=np.std(filtered_C)
      LS_sample_mean=np.mean(filtered_LS)
      LS_sample_SD=np.std(filtered_LS)
      Spec_sample_mean=np.mean(filtered_Spec)
      Spec_sample_SD=np.std(filtered_Spec)

      summary_dict['C'][proportion][post_number]=(str(obs),str(C_sample_mean),str(C_sample_SD))
      summary_dict['LS'][proportion][post_number]=(str(obs),str(LS_sample_mean),str(LS_sample_SD))
      summary_dict['Spec'][proportion][post_number]=(str(obs),str(Spec_sample_mean),str(Spec_sample_SD))

  print "Ready to write out the connectance file"
  write_summary_file(summary_dict)

  
if __name__ == '__main__':
  main()
