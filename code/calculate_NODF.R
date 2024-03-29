library(vegan)
library(igraph)

for(webtype in c("GP","SG")){
  if(webtype=="SG"){
  path<-"../data/randomised_webs/posterior/"
  propdir<-"../data/randomised_webs/detection_filter/"
  } else {
  path<-"../data/randomised_webs/gp_posterior/"
  propdir<-"../data/randomised_webs/gp_detection_filter/"
  }
  print(webtype)
  props<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  # props<-c(0.8,0.9)
  nesttable<-matrix(nrow=900,ncol=5)
  colnames(nesttable)<-c("Web","Proportion","Obs_NODF","Sample_mean","Sample_SD")
  k<-1      
  for(i in 1:100){
    print(i)
  # Obs NODF only needs calculating once
    webfile<-paste0(path,"P",i,".web",sep='')
    intlist<-read.csv(webfile,sep='\t',header=FALSE)
    G<-graph.data.frame(intlist,directed=TRUE)
    web<-as_adjacency_matrix(G,type="both",names=TRUE,sparse=FALSE)
    NODF<-nestednodf(web,order=TRUE,weighted=FALSE,wbinary=FALSE)$statistic["NODF"]
    # Extract relevant filter webs
    for(j in 1:length(props)){
      filterwebs<-as.character(list.files(path=paste0(propdir,props[j],sep=''),full.names=TRUE,pattern=paste0('P',i,'_',sep='')))
      templist<-matrix(nrow=100,ncol=1)
      # Calculate mean and SD for NODFs
      for(q in 1:length(filterwebs)){
        sampleints<-read.csv(filterwebs[q],sep='\t',header=FALSE)
        newG<-graph.data.frame(sampleints,directed=TRUE)
        newweb<-as_adjacency_matrix(newG,type="both",names=TRUE,sparse=FALSE)
        newNODF<-nestednodf(newweb,order=TRUE,weighted=FALSE,wbinary=FALSE)$statistic["NODF"]
        templist[q,1]<-newNODF
      }
      nesttable[k,1]<-paste0("P",i,".web",sep='')
      nesttable[k,2]<-props[j]
      nesttable[k,3]<-NODF
      nesttable[k,4]<-mean(templist[,1])
      nesttable[k,5]<-sd(templist[,1])
    k<-k+1

  }

  }
  write.table(nesttable,file=paste0('../data/randomised_webs/',webtype,'_NODF_table.tsv',sep=''))

}
