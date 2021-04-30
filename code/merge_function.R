mergedata = function(SC,FC,TNPCA_Structural,TNPCA_Functional,traits) {
#Vectorize SC. Note that SC is an upper triangle matrix. When i<=j, SC[i,j]=0. So when we vectorize SC (68*68 matrix), we convert it to a 67+66+...+1=2278 dimension vector.
SC_vectorize=matrix(0,1065,68*67/2)
for (k in 1:1065){
  start=1
  for (i in 1:67){
    SC_vectorize[k,start:(start+67-i)]=SC$hcp.sc.count[i,(i+1):68,k]
    start=start+68-i
  }
}
#Some subjects don't have FC data. Record subjects that don't have FC data.
missing_FC=vector()
#Note that FC is a symmetric matrix and FC[i,i]=1. So when we vectorize FC (68*68 matrix), we also convert it to a 67+66+...+1=2278 dimension vector.
FC_vectorize=matrix(0,1065,68*67/2)
for (k in 239:1065){
  if (dim(FC$hcp.cortical.fc[[k]][[1]])[1]!=0){
    start=1
    for (i in 1:67){
      FC_vectorize[k,start:(start+67-i)]=FC$hcp.cortical.fc[[k]][[1]][i,(i+1):68]
      start=start+68-i
    }
  } else{
    missing_FC=append(missing_FC,k)
  }
}
#Note that FC is a symmetric matrix and FC[i,i]=1. So when we vectorize FC (68*68 matrix), we also convert it to a 67+66+...+1=2278 dimension vector.
identical(FC$subj.list,SC$all.id)
#True
#Merge vectorized SC and vectorized FC.
df=data.frame(id=FC$subj.list,SC_vectorize,FC_vectorize)
#Record start and end location for SC data and FC data in the merged data.
SC_start=2
SC_end=SC_start+dim(SC_vectorize)[2]-1
FC_start=SC_end+1
FC_end=FC_start+dim(FC_vectorize)[2]-1
#Remove subjects who don't have FC data.
df=df[-missing_FC,]
identical(SC$all.id,TNPCA_Structural$sub.id)
#True
TNPCA_SC_Score=TNPCA_Structural$PCA.Coeff[1,1:1065,1:60]
#Remove subjects who don't have FC data.
TNPCA_SC_Score=TNPCA_SC_Score[-missing_FC,]
identical(df$id,as.vector(TNPCA_Functional$network.subject.ids))
#True
TNPCA_FC_Score=TNPCA_Functional$PCA.Coeff[1,1:1058,1:60]
#Merge vectorized SC, vectorized FC, TNPCA_SC_Score and TNPCA_FC_Score.
df=data.frame(df,TNPCA_SC_Score,TNPCA_FC_Score)
#Record start and end location for TNPCA_SC_Score and TNPCA_SC_Score in the merged data.
TNPCA_SC_Score_start=FC_end+1
TNPCA_SC_Score_end=TNPCA_SC_Score_start+dim(TNPCA_SC_Score)[2]-1
TNPCA_FC_Score_start=TNPCA_SC_Score_end+1
TNPCA_FC_Score_end=TNPCA_FC_Score_start+dim(TNPCA_FC_Score)[2]-1
#Keep subjects who have SC and FC data.
traits_data=data.frame(id=traits$hcp.subj.id,t(traits$traits.175))
traits_data=traits_data[traits$hcp.subj.id %in% df$id,]
identical(df$id,as.numeric(traits_data$id))
#True
#Merge all the data. Now the new dataframe contains id, vectorized SC, vectorized FC, TNPCA_SC_Score, TNPCA_FC_Score, and traits.
df=data.frame(df,traits_data[,-1])
#Record start and end location for traits in the merged data.
traits_start=TNPCA_FC_Score_end+1
traits_end=traits_start+175-1
return (list(df=df,SC_start=SC_start,SC_end=SC_end,FC_start=FC_start,FC_end=FC_end,
            TNPCA_SC_Score_start=TNPCA_SC_Score_start,TNPCA_SC_Score_end=TNPCA_SC_Score_end,
            TNPCA_FC_Score_start=TNPCA_FC_Score_start,TNPCA_FC_Score_end=TNPCA_FC_Score_end,
            traits_start=traits_start,traits_end=traits_end))
}