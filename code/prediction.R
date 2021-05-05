prediction = function(i,X,results) {
  y1=data[,results$traits_start+i-1]
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
  
  ##Single regression tree
  library(tree)
  fit1.single.full <- tree(y~., train1)
  testingerror_singletree1<-sqrt(sum((predict(fit1.single.full,test1)-test1$y)^2)/length(test1$y))
  
  ##Random forest,build the best possible random forest
  #ntree
  set.seed(123)
  fit1.rf.final <- randomForest(y~., train1, mtry=40, ntree=300)
  testingerror_rf1<-sqrt(sum((predict(fit1.rf.final,test1)-test1$y)^2)/length(test1))

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
  knn.model <- knn(train=train2, test=test2, cl=train2$y, k=round(sqrt(dim(train2)[1])))
  #Making confusion matrix
  cm_knn=table(test2$y,knn.model)
  #Calculate accuracy
  acc_knn=sum(diag(cm_knn))/sum(cm_knn)
  
  ##SVM
  classifier=svm(formula=y~.,data=train2,type='C-classification',kernel='linear')
  y_pred=predict(classifier,newdata=test2[-1])
  #Making confusion matrix
  cm=table(test2[,1],y_pred)
  #Calculate accuracy
  acc_svm=sum(diag(cm))/sum(cm)
  #See what happpens if we only use SC to do classification
  classifier_sc=svm(formula=y~.,data=train2[,1:162],type='C-classification',kernel='linear')
  y_pred_sc=predict(classifier_sc,newdata=test2[2:162])
  #Making confusion matrix
  cm_sc=table(test2$y,y_pred_sc)
  #Calculate accuracy
  acc_svm_sc=sum(diag(cm_sc))/sum(cm_sc)
  
  return (list(testingerror_lm1,testingerror_singletree1,testingerror_rf1,acc_knn,acc_svm,acc_svm_sc))
}
