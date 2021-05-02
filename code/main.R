library('R.matlab')
library(xlsx)
library(plotly)
library(tidyr)
library(dplyr)
library(reshape2)
library(tibble)
library(ggplot2)
library(grid)
library(ggfortify)
library(scatterplot3d)
library(profvis)
library(Hmisc)
dir = "/Users/stephanie/UNC-chapel hill/Spring2021/STOR893/893FinalProject/code"
source(paste0(dir,"/merge_function.R"))
source(paste0(dir,"/reorder_cormat.R"))
source(paste0(dir,"/get_upper_tri.R"))
source(paste0(dir,"/traits_category_correlation.R"))

setwd("/Users/stephanie/UNC-chapel hill/Spring2021/STOR893/893FinalProject/data/")

### Load Data ###
filename = system.file('SC/HCP_cortical_DesikanAtlas_SC.mat', package = 'R.matlab')

path <- system.file("mat-files", package = "R.matlab")
pathname <- file.path("SC", "HCP_cortical_DesikanAtlas_SC.mat")
SC <- readMat(pathname)

pathname <- file.path("FC", "HCP_cortical_DesikanAtlas_FC.mat")
FC <- readMat(pathname)

pathname1 <- file.path("TNPCA_Result", "TNPCA_Coeff_HCP_Structural_Connectome.mat")
pathname2 <- file.path("TNPCA_Result", "TNPCA_Coeff_HCP_Functional_Connectome.mat")
TNPCA_Structural <- readMat(pathname1)
TNPCA_Functional <- readMat(pathname2)

traits.description <- read.xlsx("traits/175traits/Details_175_Traits.xls",1) ###the Type of 175traits does not seem to correspond to the description
traits.category <- factor(traits.description$Category)
traits.type <- traits.description$Typle

pathname <- file.path("traits/175traits", "HCP_175Traits.mat")
traits<-readMat(pathname)
rownames(traits$traits.175) <-traits.description$Column_Header

### Data Cleaning & Merge Data ###

# number of NA in traits
traits.na <- data.frame(traits.name = traits.description$Column_Header,
                        percentage = rowSums(is.na(traits$traits.175))/ncol(traits$traits.175))
ggplot(top_n(traits.na,10,percentage), aes(x = reorder(traits.name, -percentage), y = percentage)) + geom_point() +
  scale_fill_gradient2(name = "correlation")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        plot.title = element_text(color="black", size=14, face="bold.italic")) +
  labs(title = paste("Top 10 traits that has most number of missing value") , x = "Traits", y = "percentage") 

## Remove traits w/ more than 40% NA
traits$traits.175 <-traits$traits.175[traits.na$percentage<=0.4,]
nrow(traits$traits.175)

results<- mergedata(SC,FC,TNPCA_Structural,TNPCA_Functional,traits)
data<- results$df

### Exploratory Data Analysis ###

## network embedding

sc.pca <- prcomp(data[,c(results$SC_start:results$SC_end)], center = TRUE,scale. = TRUE)
fc.pca <- prcomp(data[,c(results$FC_start:results$FC_end)], center = TRUE,scale. = TRUE)
tnpca.sc.sdev <- apply(data[,c(results$TNPCA_SC_Score_start:results$TNPCA_SC_Score_end)],2,sd)
tnpca.fc.sdev <- apply(data[,c(results$TNPCA_FC_Score_start:results$TNPCA_FC_Score_end)],2,sd)

var_explained_sc <- data.frame(PC= paste0("PC",1:length(sc.pca$sdev)),
                               var_explained=(sc.pca$sdev)^2/sum((sc.pca$sdev)^2)*100) %>% arrange(desc(var_explained))
var_explained_fc <- data.frame(PC= paste0("PC",1:length(fc.pca$sdev)),
                               var_explained=(fc.pca$sdev)^2/sum((fc.pca$sdev)^2)*100) %>% arrange(desc(var_explained))

var_explained_tnpca_sc <- data.frame(PC= paste0("PC",1:length(tnpca.sc.sdev)),
                               var_explained=(tnpca.sc.sdev)^2/sum((tnpca.sc.sdev)^2)*100) %>% arrange(desc(var_explained))
var_explained_tnpca_fc <- data.frame(PC= paste0("PC",1:length(tnpca.fc.sdev)),
                                     var_explained=(tnpca.fc.sdev)^2/sum((tnpca.fc.sdev)^2)*100) %>% arrange(desc(var_explained))

plot1<-ggplot(var_explained_sc[1:10,] ,aes(x = reorder(PC,-var_explained), y=var_explained, group=1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(color="black", size=14, face="bold.italic"))+
  geom_point(size=4)+
  geom_line()+
  labs(title="Scree plot: PCA on sc", x = "", y = "Percentage of explained variance")

plot2<-ggplot(var_explained_fc[1:10,] ,aes(x = reorder(PC,-var_explained), y=var_explained, group=1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(color="black", size=14, face="bold.italic"))+
  geom_point(size=4)+
  geom_line()+
  labs(title="Scree plot: PCA on fc", x = "", y = "Percentage of explained variance")

plot3<-ggplot(var_explained_tnpca_sc[1:10,] ,aes(x = reorder(PC, -var_explained), y=var_explained, group=1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(color="black", size=14, face="bold.italic"))+
  geom_point(size=4)+
  geom_line()+
  labs(title="Scree plot: TNPCA on sc", x = "", y = "Percentage of explained variance")

plot4<-ggplot(var_explained_tnpca_fc[1:10,] ,aes(x = reorder(PC, -var_explained), y=var_explained, group=1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(color="black", size=14, face="bold.italic"))+
  geom_point(size=4)+
  geom_line()+
  labs(title="Scree plot: TNPCA on fc", x = "", y = "Percentage of explained variance")

grid.newpage()
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pushViewport(viewport(layout = grid.layout(2, 2)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(1, 2))
print(plot3, vp = vplayout(2, 1))
print(plot4, vp = vplayout(2, 2))

### Correlation between PC1 v.s. traits
n <-results$traits_end-results$traits_start+1
idx = c( results$TNPCA_FC_Score_start+25,results$TNPCA_FC_Score_start+47,results$TNPCA_FC_Score_start+57,
         results$TNPCA_FC_Score_start:(results$TNPCA_FC_Score_start+2))
cor.pca.sc.traits<-matrix(0,ncol = n )
cor.pca.fc.traits<-matrix(0,ncol = n )
cor.tnpca.sc.traits<-matrix(0,ncol = n )
cor.tnpca.fc.traits<-matrix(0,ncol = n )
for(i in 1:n){
  tmp <- data.frame(sc.pca$x[,1:3], fc.pca$x[,1:3], data[,idx], traits = data[ ,results$traits_start+i-1]) %>% drop_na() 
  stdev <- apply(tmp,2,sd) 
  if (sum(stdev ==0)>=1){
    print(i)
  }
  cor.pca.sc.traits[i] <- cor(tmp[,c(1,ncol(tmp))])[1,2]
  cor.pca.fc.traits[i] <- cor(tmp[,c(4,ncol(tmp))])[1,2]
  cor.tnpca.sc.traits[i] <- cor(tmp[,c(7,ncol(tmp))])[1,2]
  cor.tnpca.fc.traits[i] <- cor(tmp[,c(10,ncol(tmp))])[1,2]
}
par(mfrow = c(2,2))
plot(1:n,cor.pca.sc.traits, type = "b", main = "cor between sc pc1 and 175 traits", 
     xlab = "traits", ylab = "correlation", ylim = c(-1,1))
plot(1:n,cor.pca.fc.traits, type = 'b', main = "cor between fc pc1 and 175 traits", 
     xlab = "traits", ylab = "correlation", ylim = c(-1,1))
plot(1:n,cor.tnpca.fc.traits, type = "b", main = "cor between sc tnpc1 and 175 traits", 
     xlab = "traits", ylab = "correlation", ylim = c(-1,1))
plot(1:n,cor.tnpca.sc.traits, type = "b", main = "cor between fc tnpc1 and 175 traits", 
     xlab = "traits", ylab = "correlation", ylim = c(-1,1))

which(abs(cor.pca.sc.traits)>0.3)
which(abs(cor.pca.fc.traits)>0.1)
which(abs(cor.tnpca.sc.traits)>0.1)
which(abs(cor.tnpca.fc.traits)>0.1)
## First 3 principal components v.s. traits
i = 51 #SC PC1 can distringuish top 100 values and bottom 100 traits in 51 
tmp <- data.frame(sc.pca$x[,1:3], fc.pca$x[,1:3], data[,idx], data[ ,results$traits_start+i-1]) %>% 
  setNames(c("sc.pc1","sc.pc2","sc.pc3","fc.pc1","fc.pc2","fc.pc3",
             "tn.sc.pc1", "tn.sc.pc2","tn.sc.pc3", "tn.fc.pc1", "tn.fc.pc2", "tn.fc.pc3","traits")) %>% drop_na() 
tmp3d<-tmp[order(tmp[,ncol(tmp)], decreasing = TRUE),] %>%
  filter(row_number() > max(row_number()) - 100 | row_number() <= 100) 
#PCA for SC v.s. traits
plot_ly(tmp3d, x = ~sc.pc1, y = ~sc.pc2, z = ~sc.pc3, color = ~traits, 
        colors = c('#BF382A', '#0C4B8E'))
#PCA for FC v.s. traits
plot_ly(tmp3d, x = ~fc.pc1, y = ~fc.pc2, z = ~fc.pc3, color = ~traits, 
        colors = c('#BF382A', '#0C4B8E'))
#TNPCA for SC v.s. traits  
plot_ly(tmp3d, x = ~tn.sc.pc1, y = ~tn.sc.pc2, z = ~tn.sc.pc3, color = ~traits, 
        colors = c('#BF382A', '#0C4B8E'))
#TNPCA for FC v.s. traits
plot_ly(tmp3d, x = ~tn.fc.pc1, y = ~tn.fc.pc2, z = ~tn.fc.pc3, color = ~traits, 
        colors = c('#BF382A', '#0C4B8E')) 

## T test PC1 v.s. traits##
idx <- seq(1,10,3)
n <- length(idx)
two.sided.t.pvalue<- matrix(0,1,n)

for(j in 1:n){
  two.sided.t.pvalue[j] = t.test(tmp3d[1:100,idx[j]],tmp3d[101:200,idx[j]])$p.value
}


## summary statistics of traits
n <-results$traits_end-results$traits_start+1
plot1<-ggplot(stack(data[,results$traits_start:(results$traits_start+49)]), aes(x = ind, y = values)) +
  geom_boxplot()
plot2<-ggplot(stack(data[,(results$traits_start+50):(results$traits_start+99)]), aes(x = ind, y = values)) +
  geom_boxplot()
plot3<-ggplot(stack(data[,(results$traits_start+100):(results$traits_start+149)]), aes(x = ind, y = values)) +
  geom_boxplot()
plot4<-ggplot(stack(data[,(results$traits_start+150):(results$traits_end)]), aes(x = ind, y = values)) +
  geom_boxplot()


grid.newpage()
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pushViewport(viewport(layout = grid.layout(2, 2)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(1, 2))
print(plot3, vp = vplayout(2, 1))
print(plot4, vp = vplayout(2, 2))

## Correlation matrix btw traits##
Main.category <- unique(traits.category)
traits_category_correlation(data,results$traits_start,traits.category, Main.category[1])
traits_category_correlation(data,results$traits_start,traits.category, Main.category[2])
traits_category_correlation(data,results$traits_start,traits.category, Main.category[3])
traits_category_correlation(data,results$traits_start,traits.category, Main.category[4])
traits_category_correlation(data,results$traits_start,traits.category, Main.category[5])
traits_category_correlation(data,results$traits_start,traits.category, Main.category[6])
traits_category_correlation(data,results$traits_start,traits.category, Main.category[7])
traits_category_correlation(data,results$traits_start,traits.category, Main.category[8])
traits_category_correlation(data,results$traits_start,traits.category, Main.category[9])


### Prediction

