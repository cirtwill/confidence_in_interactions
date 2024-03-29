import sys
import os
import networkx as nx 
import numpy as np 
import scipy as sp 

# Think I will have to do NODF in R. vegan function nestednodf(web,order=TRUE,weighted=FALSE) should get it done
# vegan takes community data as a matrix. This may be a problem. 
# Something like the following should convert interaction list to matrix:
# library(igraph);
# edge_list <- read.delim("data_file.txt", header = TRUE, sep = "\t");
# G <- graph.data.frame(edge_list,directed=FALSE);
# A <- as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE,attr="weight");

def web_reader(webdir,proportion,post_number,sample_number,nettype,site): 
  gallers=set()
  salixes=set()
  links=set()
  if sample_number==0:
    if nettype=="SG":
      f=open(webdir+'posterior/P'+str(post_number)+'_'+site+'.web','r')
    else:
      f=open(webdir+'gp_posterior/P'+str(post_number)+'_'+site+'.web','r')
  else:
    if nettype=="SG":
      f=open(webdir+'detection_filter/'+str(proportion)+'/P'+str(post_number)+'_'+site+'_'+str(sample_number)+'.web','r')
    else:
      f=open(webdir+'gp_detection_filter/'+str(proportion)+'/P'+str(post_number)+'_'+site+'_'+str(sample_number)+'.web','r')
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
  LperSalix=[] # Ranges between 18 and 37 in the posterior predicted webs
  LperGaller=[] # Ranges between 40 and 68
  for salix in salixes:
    slist=[]
    for (g,s) in links:
      if s==salix:
        slist.append(g)
    LperSalix.append(len(slist))
  for galler in gallers:
    glist=[]
    for (g,s) in links:
      if g==galler:
        glist.append(s)
    LperGaller.append(len(glist))

  return C,np.mean(LperSalix),np.mean(LperGaller)

def write_summary_file(summary_dict,nettype,site):
  f=open('../../data/Zillis_webs/'+nettype+'_Connectance_table_'+site+'.tsv','w')
  f.write('Web\tProportion\tObs_C\tSample_mean\tSample_SD\n')
  for proportion in summary_dict['C']:
    for post_number in summary_dict['C'][proportion]:
      f.write(str(post_number)+'\t'+str(proportion)+'\t'+'\t'.join(summary_dict['C'][proportion][post_number])+'\n')
  f.close()

  g=open('../../data/Zillis_webs/'+nettype+'_LS_table_'+site+'.tsv','w')
  g.write('Web\tProportion\tObs_LperSalix\tSample_mean\tSample_SD\n')
  for proportion in summary_dict['LS']:
    for post_number in summary_dict['LS'][proportion]:
      g.write(str(post_number)+'\t'+str(proportion)+'\t'+'\t'.join(summary_dict['LS'][proportion][post_number])+'\n')
  g.close()

  h=open('../../data/Zillis_webs/'+nettype+'_LG_table_'+site+'.tsv','w')
  h.write('Web\tProportion\tObs_LperGaller\tSample_mean\tSample_SD\n')
  for proportion in summary_dict['LP']:
    for post_number in summary_dict['LP'][proportion]:
      h.write(str(post_number)+'\t'+str(proportion)+'\t'+'\t'.join(summary_dict['LP'][proportion][post_number])+'\n')
  h.close()

def main():

  webdir='../../data/Zillis_webs/'

  for site in ['Zillis','Zillertal']:
    for nettype in ['SG','GP']:
      print site, nettype
      # Starting with C
      summary_dict={'C':{0.4:{},0.3:{},0.2:{},0.1:{},0.5:{},0.6:{},0.7:{},0.8:{},0.9:{},0.95:{},0.99:{}},
      'LS':{0.4:{},0.3:{},0.2:{},0.1:{},0.5:{},0.6:{},0.7:{},0.8:{},0.9:{},0.95:{},0.99:{}},
      'LP':{0.4:{},0.3:{},0.2:{},0.1:{},0.5:{},0.6:{},0.7:{},0.8:{},0.9:{},0.95:{},0.99:{}}  }

      for post_number in range(1,101):
        gallers,salixes,links=web_reader(webdir,0,post_number,0,nettype,site)
        obsC,obsLS,obsLP=calculate_original_properties(gallers,salixes,links)
        for proportion in [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]:
        # for proportion in [0.5,0.6,0.7,0.8,0.9,0.95,0.99]:
        # for proportion in [0.9,0.95,0.99]:
          filtered_C=[]
          filtered_LS=[]
          filtered_LP=[]
          for sample_number in range(1,101):
            gallers,salixes,links=web_reader(webdir,proportion,post_number,sample_number,nettype,site)
            C,LS,LP=calculate_original_properties(gallers,salixes,links)
            filtered_C.append(C)
            filtered_LS.append(LS)
            filtered_LP.append(LP)
          C_sample_mean=np.mean(filtered_C)
          C_sample_SD=np.std(filtered_C)
          LS_sample_mean=np.mean(filtered_LS)
          LS_sample_SD=np.std(filtered_LS)
          LP_sample_mean=np.mean(filtered_LP)
          LP_sample_SD=np.std(filtered_LP)

          summary_dict['C'][proportion][post_number]=(str(obsC),str(C_sample_mean),str(C_sample_SD))
          summary_dict['LS'][proportion][post_number]=(str(obsLS),str(LS_sample_mean),str(LS_sample_SD))
          summary_dict['LP'][proportion][post_number]=(str(obsLP),str(LP_sample_mean),str(LP_sample_SD))

      print "Ready to write out the connectance file"
      write_summary_file(summary_dict,nettype,site)

  
if __name__ == '__main__':
  main()
