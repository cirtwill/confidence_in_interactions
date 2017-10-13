setwd("C:/Users/kewo0002/confidence_in_interactions/data/Salix_example") 

# Turn our list of cooccurrences and interactions into a matrix
# where positive values are the number of interactions and negative
# values are the number of cooccurrences for species which never interact


# Read in data
adj_list <- read.csv("cooccur_interact_galler_salix.csv")

inter_matrix <- data.frame()

# Run through each salix-galler combination
for(salix in sort(unique(adj_list$Salix))){
  for(galler in sort(unique(adj_list$Rgaller))){
    # if they don't interact, set the matrix value to the negative of the number of times they cooccur
    if(adj_list[which(adj_list$Salix == salix & adj_list$Rgaller == galler),"interact"] == 0){
      inter_matrix[as.character(galler),as.character(salix)] <- -adj_list[which(adj_list$Salix == salix & adj_list$Rgaller == galler),"cooccur"]
    }else{
      # if they do interact, set the matrix value to the number of times they do so
      inter_matrix[as.character(galler),as.character(salix)] <- adj_list[which(adj_list$Salix == salix & adj_list$Rgaller == galler),"interact"]
      
    }
  }
}

# Make a heatmap in R
# (Will probably do this in pygrace for the MS but this is quicker to start with)
library(RColorBrewer) 
library(gplots)
#breaks
bk <- c(-50,-30, -15, -10, -7, -5,-3,-2,-1, -0.5,0.5,1,2,3,5,7,10,15,30,50)
#colors (one less than breaks)
mycols <- colorRampPalette( c("blue","white","red"))(length(bk)-1)
mycols[10]<-"white"
#plot
heatmap(as.matrix(inter_matrix), col=mycols, breaks=bk, scale="none",Rowv = NA, Colv = NA)

write.csv(inter_matrix, file = "galler_salix_heatmap_abs.csv", quote=F)



# Do the same thing with the galler_parasitoid data
adj_list_g_p <- read.csv("cooccur_interact_galler_parasit.csv")

inter_matrix_g_p <- data.frame()

for(para in sort(unique(adj_list_g_p$Rpara))){
  for(galler in sort(unique(adj_list_g_p$Rgaller))){
    if(adj_list_g_p[which(adj_list_g_p$Rpara == para & adj_list_g_p$Rgaller == galler),"interact"] == 0){
      inter_matrix_g_p[as.character(galler),as.character(para)] <- -adj_list_g_p[which(adj_list_g_p$Rpara == para & adj_list_g_p$Rgaller == galler),"cooccur"]
    }else{
      inter_matrix_g_p[as.character(galler),as.character(para)] <- adj_list_g_p[which(adj_list_g_p$Rpara == para & adj_list_g_p$Rgaller == galler),"interact"]
      
    }
  }
}
#breaks
bk <- c(-50,-30, -15, -10, -7, -5,-3,-2,-1, -0.5,0.5,1,2,3,5,7,10,15,30,50)
#colors (one less than breaks
mycols <- colorRampPalette( c("blue","white","red"))(length(bk)-1)
mycols[10]<-"white"
#plot
heatmap(as.matrix(inter_matrix_g_p), col=mycols, breaks=bk, scale="none",Rowv = NA, Colv = NA)
write.csv(inter_matrix_g_p, file = "galler_para_heatmap_abs.csv", quote=F)
