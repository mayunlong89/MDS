#Rscript
#Authors: Yunlong Ma
##The top KEGG pathways significantly enriched by genes associated with coronary artery disease(CAD):

#Part I for preparation

#First set up the directory for R
#A example:
setwd("C:\\Users\\MYL\\Desktop\\Jaccard_distance\\")

#Install R package of proxy and reshape
if(!require("proxy"))install.packages("proxy")
if(!require("reshape"))install.packages("reshape")

#load the librarys of used R packages
library(proxy)
library(reshape)

#read the .txt file that contains the 29 significant enriched pathways
#deal with the file for subsequent analysis
path1=read.delim("pathway-KEGG-CAD.txt",header=FALSE,sep='\t')
path<- path1[,-1]
row.names(path)<-path1[,1]
path2<-as.matrix(path)
number <- seq(1,length(path[,1]),1)

#Part II for calculate Jaccard distance (JD) score

#calculate the JD score among these 29 significant pathways
#Calculate the JD score for two pathways each time
for (i in number) {
  for (j in number){
    
    d1<-intersect(path2[i,],path2[j,])
    d2<-union(path2[i,],path2[j,])
    JD=length(d1)/length(d2)
    JD1=1-JD
    dnb= row.names(path2)[i]
    dna=row.names(path2)[j]
    
    DD <- c(dna,dnb,JD)
    Y =t(DD)
		
    #write out the results of all JD scores for these pathways
    write.table(Y, file="Jaccard_distance_final_CAD.txt", 
                append=TRUE, row.names=FALSE,col.names = FALSE)  
  }
}

#Part III for re-load data

#Re-read the JD score from the Jaccard_distance_final_CAD.txt
Data_11 <- read.table("Jaccard_distance_final_CAD.txt",header = FALSE)

#Read the annotation file of name_KEGG_CAD.txt
anno <- read.table("name_KEGG_CAD.txt",header=TRUE)

#Re-name the collums of Data_11 file
#Annotate the variable name 
#Reshape the formate of data matrix 
colnames(Data_11) <- c("time","variable","value")
Data_11$time <- anno$ID[match(Data_11$time,anno$Name)]
data_cast<-cast(Data_11,time~variable)
data_cast$time <- anno$Name[match(data_cast$time,anno$ID)]

#Part IV for MDS analysis

#Calculate the multidimensional scaling (MDS) value 
#change the data type into dist type
data_cast_2 <-as.dist(data_cast)
voles.mds=cmdscale(data_cast_2,k=9,eig=T)
sum(abs(voles.mds$eig[1:2]))/sum(abs(voles.mds$eig))
sum((voles.mds$eig[1:2])^2)/sum((voles.mds$eig)^2)
MDS_component_1 = voles.mds$points[,1]
MDS_component_2 = voles.mds$points[,2]
plot(MDS_component_1,MDS_component_2,cex=1.8,pch=21,col = "black", bg = "green")
text(MDS_component_1,MDS_component_2, row.names(voles.mds$points),cex=0.6,pos=4

#Part V for visualization

#use yellow and red color to define the signficant of each pathway
dat_1 <- read.table("anno_KEGG_CAD.txt",header=TRUE)
data_temp <- data.frame(dat_1,MDS_component_1,MDS_component_2)
cPal <- colorRampPalette(c('yellow','red'))
datCol<-cPal(18)[as.numeric(cut(dat_1$Zscore,breaks = 18))]

#bubble diagram
symbols(MDS_component_1,MDS_component_2,circle = data_temp$Zscore,inches=0.15,col = "black", bg = datCol)
text(MDS_component_1,MDS_component_2, row.names(voles.mds$points),cex=0.6,pos=4)

#End
