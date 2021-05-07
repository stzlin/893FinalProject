library('R.matlab')
library(xlsx)
library(plotly)
library(car)
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
library(tree)
library(randomForest)
library(class)
library(e1071)
dir = "/Users/stephanie/UNC-chapel hill/Spring2021/STOR893/893FinalProject/code"
source(paste0(dir,"/merge_function.R"))
source(paste0(dir,"/reorder_cormat.R"))
source(paste0(dir,"/get_upper_tri.R"))
source(paste0(dir,"/traits_category_correlation.R"))
source(paste0(dir,"/prediction.R"))

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


### The data description file has been changed so that the index matches the one in traits extraction file.
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
  labs(title = paste("Top 10 traits that has most number of NAs") , x = "Traits", y = "percentage") 


## Remove traits w/ more than 40% NA
traits$traits.175 <-traits$traits.175[traits.na$percentage<=0.4,]
traits.description <- traits.description[traits.na$percentage<=0.4,]
traits.category <- traits.category[traits.na$percentage<=0.4]
traits.type <- traits.type[traits.na$percentage<=0.4]

results<- mergedata(SC,FC,TNPCA_Structural,TNPCA_Functional,traits)
data<- results$df

## Remove traits with constant value
idx <- which(apply(data[,results$traits_start:results$traits_end], 2, function(x){sum(is.na(unique(x))==FALSE)}) ==1)
data <- data[,-(results$traits_start+idx-1)]
results$traits_end<-results$traits_end-1
traits.description <- traits.description[-idx,]
traits.category <- traits.category[-idx]
traits.type <- traits.type[-idx]

ncol(data[,results$traits_start:results$traits_end])
table(traits.type)
tmp_traits<-melt(table(traits.category))

ggplot(data = tmp_traits, aes(x= reorder(traits.category,-value), y = value, fill = value)) + geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient2(low = "white", 
                       high = "purple") +
  theme(axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1), 
        plot.title = element_text(color="black", size=14, face="bold.italic")) +
  labs(title = "Traits Category", x = "", y = "# of traits", fill = "# of traits") 
### Exploratory Data Analysis ###
## network embedding

sc.pca <- prcomp(data[,c(results$SC_start:results$SC_end)], center = TRUE,scale. = TRUE)
fc.pca <- prcomp(data[,c(results$FC_start:results$FC_end)], center = TRUE,scale. = TRUE)

var_explained_sc <- data.frame(PC= paste0("PC",1:length(sc.pca$sdev)),
                               var_explained=(sc.pca$sdev)^2/sum((sc.pca$sdev)^2)*100) %>% arrange(desc(var_explained))
var_explained_fc <- data.frame(PC= paste0("PC",1:length(fc.pca$sdev)),
                               var_explained=(fc.pca$sdev)^2/sum((fc.pca$sdev)^2)*100) %>% arrange(desc(var_explained))


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

grid.newpage()
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pushViewport(viewport(layout = grid.layout(1, 2)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(1, 2))

### Correlation between PC1 v.s. continuous traits
idx_cont_traits <- which(traits.type=="Continuous") + results$traits_start-1
idx_pca = c( results$TNPCA_SC_Score_start:(results$TNPCA_SC_Score_start+2),
         results$TNPCA_FC_Score_start:(results$TNPCA_FC_Score_start+2))
tmp_pca <- data.frame(sc.pca$x[,1:3], fc.pca$x[,1:3], data[,idx_pca])
n <- length(idx_cont_traits)

cor.pca.sc.traits<-matrix(0,ncol = n)
cor.pca.fc.traits<-matrix(0,ncol = n )
cor.tnpca.sc.traits<-matrix(0,ncol = n )
cor.tnpca.fc.traits<-matrix(0,ncol = n )
for(i in 1:n){
  if(traits.type[i] == "Continuous"){
    tmp <- data.frame(tmp_pca, traits = data[,idx_cont_traits[i]]) %>% drop_na() 
    stdev <- apply(tmp,2,sd) 
    if (sum(stdev ==0)>=1){
      print(i)
    }
    cor.pca.sc.traits[i] <- cor(tmp[,c(1,ncol(tmp))])[1,2]
    cor.pca.fc.traits[i] <- cor(tmp[,c(4,ncol(tmp))])[1,2]
    cor.tnpca.sc.traits[i] <- cor(tmp[,c(7,ncol(tmp))])[1,2]
    cor.tnpca.fc.traits[i] <- cor(tmp[,c(10,ncol(tmp))])[1,2]
  }
}
par(mfrow = c(2,2))
plot(1:n,cor.pca.sc.traits, type = "p", main = "PC1 for SC v.s. continuous traits", 
     xlab = "continuous traits", ylab = "correlation", ylim = c(-0.5,0.5))
abline(h = 0.3, col ="red", lty = 2)
abline(h = -0.3, col ="red", lty = 2)
plot(1:n,cor.pca.fc.traits, type = 'p', main = "PC1 for FC v.s. continuous traits", 
     xlab = "continuous traits", ylab = "correlation", ylim = c(-0.5,0.5))
abline(h = 0.3, col ="red", lty = 2)
abline(h = -0.3, col ="red", lty = 2)
plot(1:n,cor.tnpca.sc.traits, type = "p", main = "TNPC1 for SC v.s. continuous traits", 
     xlab = "continuous traits", ylab = "correlation", ylim = c(-0.5,0.5))
abline(h = 0.3, col ="red", lty = 2)
abline(h = -0.3, col ="red", lty = 2)
plot(1:n,cor.tnpca.fc.traits, type = "p", main = "TNPC1 for FC v.s. continuous traits", 
     xlab = "continuous traits", ylab = "correlation", ylim = c(-0.5,0.5))
abline(h = 0.3, col ="red", lty = 2)
abline(h = -0.3, col ="red", lty = 2)

which(abs(cor.pca.sc.traits)>0.3)
which(abs(cor.pca.fc.traits)>0.1)
which(abs(cor.tnpca.sc.traits)>0.3)
which(abs(cor.tnpca.fc.traits)>0.1)
colnames(data)[results$traits_start+which(abs(cor.pca.sc.traits)>0.3)-1]
traits.category[which(abs(cor.pca.sc.traits)>0.3)]
## First 3 principal components v.s. traits
i = 49 #SC PC1 can distringuish top 100 values and bottom 100 traits in 50
tmp <- data.frame(sc.pca$x[,1:3], fc.pca$x[,1:3], data[,idx_pca], data[ ,results$traits_start+i-1]) %>% 
  setNames(c("sc.pc1","sc.pc2","sc.pc3","fc.pc1","fc.pc2","fc.pc3",
             "tn.sc.pc1", "tn.sc.pc2","tn.sc.pc3", "tn.fc.pc1", "tn.fc.pc2", "tn.fc.pc3", colnames(data)[results$traits_start+i-1])) %>% drop_na() 
tmp3d<-tmp[order(tmp[,ncol(tmp)], decreasing = TRUE),] %>%
  filter(row_number() > max(row_number()) - 100 | row_number() <= 100) 
par(mfrow = c(2,2))
#PCA for SC v.s. traits
plot_ly(tmp3d, x = ~sc.pc1, y = ~sc.pc2, z = ~sc.pc3, color = ~GaitSpeed_Comp,  
        colors = c('#BF382A', '#0C4B8E')) %>% add_markers() %>%
  layout(title = 'Motor: GaitSpeed_Comp')

#PCA for FC v.s. traits
plot_ly(tmp3d, x = ~fc.pc1, y = ~fc.pc2, z = ~fc.pc3, color = ~GaitSpeed_Comp, 
        colors = c('#BF382A', '#0C4B8E')) %>% add_markers() %>%
  layout(title = 'Motor: GaitSpeed_Comp')
#TNPCA for SC v.s. traits  
plot_ly(tmp3d, x = ~tn.sc.pc1, y = ~tn.sc.pc2, z = ~tn.sc.pc3, color = ~GaitSpeed_Comp, 
        colors = c('#BF382A', '#0C4B8E')) %>% add_markers() %>%
  layout(title = 'Motor: Endurance_AgeAdj')
#TNPCA for FC v.s. traits
plot_ly(tmp3d, x = ~tn.fc.pc1, y = ~tn.fc.pc2, z = ~tn.fc.pc3, color = ~GaitSpeed_Comp, 
        colors = c('#BF382A', '#0C4B8E')) %>% add_markers() %>%
  layout(title = 'Motor: GaitSpeed_Comp')

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

#Find the number of SC PCA scores we need to keep that can explain at least 60% of total variance.
x_sc=which(cumsum(var_explained_sc$var_explained)>60)[1]
#Find the number of FC PCA scores we need to keep that can explain at least 60% of total variance.
x_fc=which(cumsum(var_explained_fc$var_explained)>60)[1]
#Use SC PCA score, FC PCA score, TNPCA_SC_Score and TNPCA_FC_Score for prediction.
X=data.frame(sc.pca$x[,1:x_sc],fc.pca$x[,1:x_fc],
data[,results$TNPCA_SC_Score_start:results$TNPCA_SC_Score_end],
data[,results$TNPCA_FC_Score_start:results$TNPCA_FC_Score_end])

###SC PC1 can distringuish top 100 values and bottom 100 values in 51th trait. So we first try to predict 51th trait.
i1 = 51
y1=data[,results$traits_start+i1-1]
data1=data.frame(y=y1,X)
#Remove subjects that have NA in 51th trait.
data1=drop_na(data1)

##Regression
#Make a training data and testing data (approximately 2/3 observations as training data)
set.seed(66)
index1=sort(sample(nrow(data1),round(nrow(data1)*2/3)))
train1<-data1[index1,]
test1<-data1[-index1,]

##Linear regression
lm1<-lm(y~.,train1)
testingerror_lm1<-sqrt(sum((predict(lm1,test1)-test1$y)^2)/length(test1$y))
testingerror_lm1

##Single regression tree
fit1.single.full <- tree(y~., train1)
testingerror_singletree1<-sqrt(sum((predict(fit1.single.full,test1)-test1$y)^2)/length(test1$y))
testingerror_singletree1

##Random forest,build the best possible random forest
#ntree
set.seed(123)
fit1.rf <- randomForest(y~., train1, mtry=10, ntree=1000)
plot(fit1.rf, col="red", pch=16, type="p", main="default plot")
#We may need 300 trees.
#Choose "mtry"
rf1.error.p <- 1:50  # set up a vector of length 20
for (p in 1:50)  # repeat the following code inside { } 20 times
{
  fit1.rf <- randomForest(y~., train1, mtry=p, ntree=300)
  rf1.error.p[p] <- fit1.rf$mse[300]  # collecting oob mse based on 250 trees
}
rf1.error.p
plot(1:50, rf1.error.p, pch=16,
     xlab="mtry",
     ylab="OOB mse of mtry")
lines(1:50, rf1.error.p)
#We chose ntree=300 and mtry=k1.
k1=which.min(rf1.error.p)
#Final random forest model.
fit1.rf.final <- randomForest(y~., train1, mtry=k1, ntree=300)
testingerror_rf1<-sqrt(sum((predict(fit1.rf.final,test1)-test1$y)^2)/length(test1))
testingerror_rf1


##Classification
#Reorder data based on 51th trait value
data2=data1[order(data1[,1], decreasing = TRUE),]
#Create High/Low (1/0) label
data2$y=0
data2$y[1:round(dim(data2)[1]/2)]=1
#Encoding the target feature
data2$y=factor(data2$y,levels=c(0,1))
#Make a training data and testing data (approximately 2/3 observations as training data)
set.seed(88)
index2=sort(sample(nrow(data2),round(nrow(data2)*2/3)))
train2<-data2[index2,]
test2<-data2[-index2,]
#Feature scaling
train2[-1]=scale(train2[-1])
test2[-1]=scale(test2[-1])

##knn
sqrt(dim(train2)[1])
#Square root of number of observations in training dataset is around 26.55, therefore we’ll create two models. One with ‘K’ value as 26 and the other model with a ‘K’ value as 27.
knn.26 <- knn(train=train2, test=test2, cl=train2$y, k=26)
knn.27 <- knn(train=train2, test=test2, cl=train2$y, k=27)
#Making confusion matrix
cm_knn26=table(test2$y,knn.26)
#Calculate accuracy
acc_knn26=sum(diag(cm_knn26))/sum(cm_knn26)
cm_knn27=table(test2$y,knn.27)
acc_knn27=sum(diag(cm_knn27))/sum(cm_knn27)
acc_knn26
acc_knn27
#Optimization, find value "K" that model has highest accuracy.
i=1
k.optm=1
for (i in 1:30){
knn.mod <- knn(train=train2, test=test2, cl=train2$y, k=i)
cm_knn=table(test2$y,knn.mod)
k.optm[i] <- sum(diag(cm_knn))/sum(cm_knn)
k=i
}
which.max(k.optm)
#Accuracy plot
plot(k.optm, type="b", xlab="K-Value",ylab="Accuracy level")

##SVM
#Linear kernel
classifier=svm(formula=y~.,data=train2,type='C-classification',kernel='linear')
y_pred=predict(classifier,newdata=test2[-1])
#Making confusion matrix
cm=table(test2[,1],y_pred)
#Calculate accuracy
acc_svm=sum(diag(cm))/sum(cm)
#Radial kernel
classifier=svm(formula=y~.,data=train2,type='C-classification',kernel='radial')
y_pred=predict(classifier,newdata=test2[-1])
#Making confusion matrix
cm=table(test2[,1],y_pred)
#Calculate accuracy
acc_svm_radial=sum(diag(cm))/sum(cm)


##Prediction results for other traits. (testing error for linear regression, single tree, random forest and classification accuracy for knn and svm)
#Also try to only use SC PCA score, FC PCA score, TNPCA_SC_Score and TNPCA_FC_Score for prediction.
X_sc=X=data.frame(sc.pca$x[,1:x_sc],
data[,results$TNPCA_SC_Score_start:results$TNPCA_SC_Score_end])
continuous_variable_index=which((traits.type=='Continuous'))
n=length(continuous_variable_index)
SC_FC_regression=data.frame(matrix(,nrow=n,ncol=3))
colnames(SC_FC_regression)=c('testing_error_lm','testing_error_singletree','testing_error_rf')
SC_regression=data.frame(matrix(,nrow=n,ncol=3))
colnames(SC_regression)=c('testing_error_lm','testing_error_singletree','testing_error_rf')
SC_FC_classification=data.frame(matrix(,nrow=n,ncol=3))
colnames(SC_FC_classification)=c('acc_knn','acc_svm_linear','acc_svm_radial')
SC_classification=data.frame(matrix(,nrow=n,ncol=3))
colnames(SC_classification)=c('acc_knn','acc_svm_linear','acc_svm_radial')
for (i in 1:n){
    pred=prediction(continuous_variable_index[i],X,results)
    pred_SC=prediction(continuous_variable_index[i],X_sc,results)
    print(i)
    SC_FC_regression[i,1]=pred[[1]]
    SC_FC_regression[i,2]=pred[[2]]
    SC_FC_regression[i,3]=pred[[3]]
    SC_FC_classification[i,1]=pred[[4]]
    SC_FC_classification[i,2]=pred[[5]]
    SC_FC_classification[i,3]=pred[[6]]
    SC_regression[i,1]=pred_SC[[1]]
    SC_regression[i,2]=pred_SC[[2]]
    SC_regression[i,3]=pred_SC[[3]]
    SC_classification[i,1]=pred_SC[[4]]
    SC_classification[i,2]=pred_SC[[5]]
    SC_classification[i,3]=pred_SC[[6]]
}
#Find traits that have a high classification accuracy
continuous_variable_index[SC_FC_classification$acc_svm_radial>0.7]
continuous_variable_index[SC_FC_classification$acc_svm_radial>0.65]
