library(glmnet)
#initialize the input matrix to the regression
input.matrix <- matrix(0,nrow=dim(cna)[2],ncol=dim(cna)[1])
input.matrix[t(cna) < -1] <- 1
input.matrix[t(mm.sig) > 0] <- 1
colnames(input.matrix) <- rownames(cna)
rownames(input.matrix) <- colnames(cna)
smok = as.numeric(as.character(smoking[,"patient.number_pack_years_smoked"]));
smok[is.na(smok)] = mean(smok[!is.na(smok)]);
smokhist = as.factor(smoking[,"patient.tobacco_smoking_history"])
smokhist = as.data.frame(model.matrix(~smokhist)[,2:dim(model.matrix(~smokhist))[2]]);
age = as.numeric(as.character(clinical[match(tolower(substr(colnames(cna),1,12)),clinical[,"patient.bcr_patient_barcode"]),"patient.age_at_initial_pathologic_diagnosis"]));
age[is.na(age)] = median(age[!is.na(age)]);

input.matrix = cbind(log(age+1),smokhist,input.matrix);
colnames(input.matrix)[1:5] = c("age","smokhist1","smokhist2","smokhist3","smokhist4")
#k = number of cross validation for signature performance, holdout = number of samples in test set
k = 1000
holdout = 50
dnaintegrity.cv <- matrix(0,nrow=dim(input.matrix)[2],ncol=k)
rownames(dnaintegrity.cv) <- colnames(input.matrix)
cv.lasso.object <- list()
real.train <- matrix(0,nrow=dim(input.matrix)[1]-holdout,ncol=k)
predicted.train <- matrix(0,nrow=dim(input.matrix)[1]-holdout,ncol=k)
real.test <- matrix(0,nrow=holdout,ncol=k)
predicted.test <- matrix(0,nrow=holdout,ncol=k)
for(i in 1:k){
  sample.temp <- sample(1:dim(input.matrix)[1],dim(input.matrix)[1]-holdout)
  cv.lasso.object[[i]] <- cv.glmnet(as.matrix(input.matrix[sample.temp,]),log2(patient.mutburden+1)[sample.temp])
  lasso.object <- glmnet(as.matrix(input.matrix[sample.temp,]),log2(patient.mutburden+1)[sample.temp],family='poisson',lambda = cv.lasso.object[[i]]$lambda.min)
  dnaintegrity.cv[,i] <- lasso.object$beta[,1]
  real.train[,i] <- log2(patient.mutburden+1)[sample.temp]
  predicted.train[,i] <- predict(lasso.object,as.matrix(input.matrix[sample.temp,]))
  real.test[,i] <- log2(patient.mutburden+1)[setdiff(1:dim(input.matrix)[1],sample.temp)]
  predicted.test[,i] <- predict(lasso.object,as.matrix(input.matrix[setdiff(1:dim(input.matrix)[1],sample.temp),]))
  if(i %% 50 == 0){print(i)}
}

#checking to see how well the test/training predictions correlate with the actual mutation burdens
cor.train <- vector(mode='numeric',k)
cor.tester <- vector(mode='numeric',k)
for(i in 1:k){
  cor.tester[i] <- cor(predicted.test[,i],real.test[,i])
  cor.train[i] <- cor(real.train[,i],predicted.train[,i])
}

#Getting Relevant features
feats = list();
num_feats = matrix(0,k,1);
for (i in 1:k){
  feats[[i]] = row.names(coef(cv.lasso.object[[i]]))[as.numeric(coef(cv.lasso.object[[i]]))>0];
  num_feats[i] = length(feats[[i]])-1;
}
AllFeats = unlist(feats);

uniqFeats = unique(AllFeats);
AllFeats.counts = matrix(0,length(uniqFeats),1);
AllFeats.probs = matrix(0,length(uniqFeats),1);
row.names(AllFeats.counts) = uniqFeats;
row.names(AllFeats.probs) = uniqFeats;
for (i in 1:length(uniqFeats)){
  AllFeats.counts[i] = length(grep(uniqFeats[i],AllFeats));
  AllFeats.probs[i] = sum(dbinom(AllFeats.counts[i]:k,k,mean(num_feats/dim(input.matrix)[2])))
}

output = cbind(AllFeats.counts,AllFeats.probs)
colnames(output) = c("counts","probability")
write.csv(output,"~/Desktop/feature_probs.csv")

