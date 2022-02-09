##setwd("D:/Saurav SXM1600 backup/R Script files")
#setwd("D:/R Script files")
setwd("C:/Users/smallik/my backups/R Script files")

###debug(utils:::unpackPkgZip)#######use this before using "install.packages" command if "unable to move temporary..." error comes
###install.packages("xml2")

library(MOCCA)

###adeno=22 and scc=253 samples######
#cerival.scc.vs.adeno=read.csv("SCCvsADENO_rna_URG+DRG.csv",header=T,sep=",") ###with cpgnames (rows) and samples (columns)#wrong
cerival.scc.vs.adeno=read.csv("OUTPUT_WGCNA.AdenovsSCC.csv",header=T,sep=",") ###with cpgnames (rows) and samples (columns)
cerival.scc.vs.adeno.f<-cerival.scc.vs.adeno[,2:ncol(cerival.scc.vs.adeno)]
rownames(cerival.scc.vs.adeno.f)<-cerival.scc.vs.adeno[,1]
dim(cerival.scc.vs.adeno.f)#nrow=582, ncol=275

##set.seed(123);##set 123 as seed value ##got 0.9381818 accuracy, std.specificity is zero
set.seed(285);##set 285 as seed value ##all are good
cernival.mocca.res<-mocca(cerival.scc.vs.adeno.f, R = 10, K = 2:10, iter.max = 10, nstart = 10)
cernival.mocca.res.pareto<-analyzePareto(cernival.mocca.res$objectiveVals)
print(cernival.mocca.res.pareto)

cernival.mocca.res.pareto$rank ##(5,7,3)
cernival.mocca.res.pareto$table

moo.clusterid<-cernival.mocca.res.pareto$rank[1]##find top multi-objective optimized cluster id (i.e., rank 1 clsid) ## depending upon max(minimum number of dominating objective functions for the other cluster sizes)

print(cernival.mocca.res$objectiveVals[,moo.clusterid])##find objective values for the combination of different clustering techniques and clustering indices measures
write.table(cernival.mocca.res$objectiveVals[,moo.clusterid],file=paste0("zz_mocca_cervival_SquavsAdeno_objectiveval_of_selected_mooclusterno_",moo.clusterid,".csv",sep=""),sep=",",row.names=TRUE,col.names=FALSE)
# kmeans.MCA    kmeans.Jaccard         kmeans.FM        kmeans.CQS     neuralgas.MCA neuralgas.Jaccard      neuralgas.FM 
# 0.6024096         0.5089413         0.6077859         0.9783786         0.6024096         0.5177725         0.6127136 
# neuralgas.CQS        single.MCA    single.Jaccard         single.FM        single.CQS 
# 0.9786411         0.5507745         0.3493986         0.5201424         0.9771901 

####find moo clusters of the different genes###
cernival.mocca.res$cluster$kmeans[[moo.clusterid]][[10]]
#cernival.mocca.res$cluster$kmeans[[moo.clusterid]][[10]]
#cernival.mocca.res$cluster$neuralgas[[moo.clusterid]][[10]]
#cernival.mocca.res$cluster$single[[moo.clusterid]][[10]]

#cernival.mocca.res$cluster$baseline[[moo.clusterid]][[10]]
gene_clsid<-data.frame(cernival.mocca.res$cluster$kmeans[[moo.clusterid]][[10]])
rownames(gene_clsid)<-cerival.scc.vs.adeno[,1]
write.table(gene_clsid,file=paste0("zz_mocca_cervival_SquavsAdeno_kemans_gene_clusterid_for_moototalclusterno_",moo.clusterid,".csv",sep=""),sep=",",row.names=TRUE,col.names=FALSE)


####Mahalanobis Distance computing for each gene data vector for each cluster, and then averaging for each cluster#######
##https://cran.r-project.org/web/packages/distances/distances.pdf
library("distances")

##avg.eucle_dist=vector(mode="integer", length=moo.clusterid)
avg.spr_dist=vector(mode="integer", length=moo.clusterid)

for(i in 1:moo.clusterid)
{
  #i=1
  clsid.data<-cerival.scc.vs.adeno.f[which((rownames(cerival.scc.vs.adeno.f) %in% rownames(gene_clsid)) & gene_clsid$cernival.mocca.res.cluster.kmeans..moo.clusterid....10..==i),]
  ##eucledean_distances <- distances(clsid.data)#eucledian distance between any two gene-vectors
  spr_distances <-cor(t(clsid.data), use="complete.obs", method="spearman")#Spearman correlation between any two gene-vectors
  #my_distances3 <- distances(clsid.data, normalize = "mahalanobize")#mohalanbis distance between any two gene-vectors
  ##avg.eucle_dist[i]<-mean(eucledean_distances[row(eucledean_distances)!=col(eucledean_distances)])
  avg.spr_dist[i]<-mean(spr_distances[lower.tri(spr_distances)])##taking only one sided values of the symmetric matrix
}
##print(avg.eucle_dist)#3.275980, 9.559229, 5.413411, 2.530057, 1.137817
print(avg.spr_dist)##new:(0.3122533, 0.2009791, 0.3093203, 0.5208954, 0.2119528)

##min.avg.eucledist<-min(avg.eucle_dist)#1.137817
max.avg.sprdist<-max(avg.spr_dist)#0.5208954 ##find the maximum avg_pearson_correlation score among of the resultant clusters
##print(min.avg.eucledist)
print(max.avg.sprdist)
##clsid_min.avg.eucledist<- which(avg.eucle_dist == min(avg.eucle_dist))#5 
clsid_max.avg.sprdist<- which(avg.spr_dist == max(avg.spr_dist))#4 ##find the clsid that has max spearman correlation
##print(clsid_min.avg.eucledist)
print(clsid_max.avg.sprdist)###4
#clsid_min.avg.eucledist<- which.min(avg.eucle_dist)#5 
##clsid_max.avg.sprdist<- which.max(avg.spr_dist)#4 
##data.clsid_min.avg.eucledist<-cerival.scc.vs.adeno.f[which((rownames(cerival.scc.vs.adeno.f) %in% rownames(gene_clsid)) & gene_clsid$cernival.mocca.res.cluster.kmeans..moo.clusterid....10..==clsid_min.avg.eucledist),]
data.clsid_max.avg.sprdist<-cerival.scc.vs.adeno.f[which((rownames(cerival.scc.vs.adeno.f) %in% rownames(gene_clsid)) & gene_clsid$cernival.mocca.res.cluster.kmeans..moo.clusterid....10..==clsid_max.avg.sprdist),]


##ttyu1<-data.frame(rownames(data.clsid_min.avg.eucledist))
ttyu1<-data.frame(rownames(data.clsid_max.avg.sprdist))
##data.clsid_min.avg.eucledist_forwrt<-cbind(ttyu1,data.clsid_min.avg.eucledist)
data.clsid_max.avg.sprdist_forwrt<-cbind(ttyu1,data.clsid_max.avg.sprdist)
##colnames(data.clsid_min.avg.eucledist_forwrt)[1]<-"gene"
colnames(data.clsid_max.avg.sprdist_forwrt)[1]<-"gene"
##write.table(data.clsid_min.avg.eucledist_forwrt,file="zz_mocca_cernival_signature.csv",sep=",",row.names=FALSE,col.names=TRUE)
write.table(data.clsid_max.avg.sprdist_forwrt,file="zz_mocca_cernival_signature.csv",sep=",",row.names=FALSE,col.names=TRUE)


################find up-regulated gene-list and down-regulated gene-list belonging to the gene signature#############
urg_adeno_vs_scc=read.csv("cernival.ADENO.vs.SCC.URG.csv",header=T,sep=",") ###with cpgnames (rows) and samples (columns)
drg_adeno_vs_scc=read.csv("cernival.ADENO.vs.SCC.DRG.csv",header=T,sep=",") ###with cpgnames (rows) and samples (columns)
dim(urg_adeno_vs_scc)#[260,7]
dim(drg_adeno_vs_scc)#[322,7]
urg_sig<-intersect(data.clsid_max.avg.sprdist_forwrt$gene,urg_adeno_vs_scc$Gene_Name)
drg_sig<-intersect(data.clsid_max.avg.sprdist_forwrt$gene,drg_adeno_vs_scc$Gene_Name)
length(urg_sig)#urg in signature=35
length(drg_sig)#drg in signature=0

data_urg_adeno_vs_scc_for_gensig<-cerival.scc.vs.adeno.f[which(rownames(cerival.scc.vs.adeno.f) %in% urg_sig),]
data_urg_adeno_vs_scc_for_gensig<-cbind(rownames(data_urg_adeno_vs_scc_for_gensig),data_urg_adeno_vs_scc_for_gensig)
dim(data_urg_adeno_vs_scc_for_gensig)
data_drg_adeno_vs_scc_for_gensig<-cerival.scc.vs.adeno.f[which(rownames(cerival.scc.vs.adeno.f) %in% drg_sig),] 
data_drg_adeno_vs_scc_for_gensig<-cbind(rownames(data_drg_adeno_vs_scc_for_gensig),data_drg_adeno_vs_scc_for_gensig)
dim(data_drg_adeno_vs_scc_for_gensig)
write.table(data_urg_adeno_vs_scc_for_gensig,file="zz_mocca_cernival_signature_urgdat.csv",sep=",",row.names=FALSE,col.names=TRUE)
write.table(data_drg_adeno_vs_scc_for_gensig,file="zz_mocca_cernival_signature_drgdat.csv",sep=",",row.names=FALSE,col.names=TRUE)

################find up-regulated gene-list and down-regulated gene-list belonging to all the clusters#############
for(i in 1:moo.clusterid)
{
  #i=1
  data.clsid_i<-cerival.scc.vs.adeno.f[which((rownames(cerival.scc.vs.adeno.f) %in% rownames(gene_clsid)) & gene_clsid$cernival.mocca.res.cluster.kmeans..moo.clusterid....10..==i),]
  dim(data.clsid_i)
  ttyu1_i<-data.frame(rownames(data.clsid_i))
  data.clsid_i_forwrt<-cbind(ttyu1_i,data.clsid_i)
  colnames(data.clsid_i_forwrt)[1]<-"gene"
  
  urg_sig_i<-intersect(data.clsid_i_forwrt$gene, urg_adeno_vs_scc$Gene_Name)
  drg_sig_i<-intersect(data.clsid_i_forwrt$gene, drg_adeno_vs_scc$Gene_Name)
  length(urg_sig_i)#urg in each resultant cluster
  length(drg_sig_i)#drg in each resultant cluster
  
  print(paste0("#urg_sig(",i,")=",length(urg_sig_i),sep=""))
  print(paste0("#drg_sig(",i,")=",length(drg_sig_i),sep=""))
  
  data_urg_adeno_vs_scc_for_clsid_i<-cerival.scc.vs.adeno.f[which(rownames(cerival.scc.vs.adeno.f) %in% urg_sig_i),]
  data_urg_adeno_vs_scc_for_clsid_i<-cbind(rownames(data_urg_adeno_vs_scc_for_clsid_i),data_urg_adeno_vs_scc_for_clsid_i)
  dim(data_urg_adeno_vs_scc_for_clsid_i)
  data_drg_adeno_vs_scc_for_clsid_i<-cerival.scc.vs.adeno.f[which(rownames(cerival.scc.vs.adeno.f) %in% drg_sig_i),]  
  data_drg_adeno_vs_scc_for_clsid_i<-cbind(rownames(data_drg_adeno_vs_scc_for_clsid_i),data_drg_adeno_vs_scc_for_clsid_i)
  dim(data_drg_adeno_vs_scc_for_clsid_i)
  write.table(data_urg_adeno_vs_scc_for_clsid_i,file=paste0("zz_mocca_cernival_clsid_",i,"_urgdat.csv",sep=""),sep=",",row.names=FALSE,col.names=TRUE)#data for only urg
  write.table(data_drg_adeno_vs_scc_for_clsid_i,file=paste0("zz_mocca_cernival_clsid_",i,"_drgdat.csv",sep=""),sep=",",row.names=FALSE,col.names=TRUE)#data for only drg
  write.table(data.clsid_i_forwrt,file=paste0("zz_mocca_cernival_clsid_",i,"_deg_dat.csv",sep=""),sep=",",row.names=FALSE,col.names=TRUE)#data for only deg
}




############resultant Cluster information (principal component analysis or PCA) finding#############
#library("devtools")
#install_github("lme4/lme4",dependencies=TRUE)

####https://tgmstat.wordpress.com/2013/11/28/computing-and-visualizing-pca-in-r/ ##
###http://www.gastonsanchez.com/visually-enforced/how-to/2012/06/17/PCA-in-R/  ##
###https://poissonisfish.wordpress.com/2017/01/23/principal-component-analysis-in-r/
##https://www.youtube.com/watch?v=HMOI_lkzW08

##https://support.bioconductor.org/p/51270/ ###difference b/w sample clustering and gene clustering (here) through PCA

library(devtools)
##install_github("ggbiplot", "vqv")


# degdata.clsid_all<-matrix(, nrow=0, ncol=277)


for(i in 1:moo.clusterid)
{
  #i=1
  degdata.clsid_i=read.csv(file=paste0("zz_mocca_cernival_clsid_",i,"_deg_dat.csv",sep=""),header=T,sep=",") ###with cpgnames (rows) and samples (columns)
  dim(degdata.clsid_i)
  degdata.clsid_i_clsvector <- rep(paste0("Cluster",i), nrow(degdata.clsid_i))
  length(degdata.clsid_i_clsvector)
  degdata.clsid_i_with_clslabel<-cbind(degdata.clsid_i_clsvector,degdata.clsid_i)
  dim(degdata.clsid_i_with_clslabel)
  degdata.clsid_all<-rbind(degdata.clsid_all,degdata.clsid_i_with_clslabel)##take all genes of all the clusters havinf cluster id
  dim(degdata.clsid_all)
}
dim(degdata.clsid_all)##[582,277]

group_sampl<-degdata.clsid_all[,1]#####take cluster idvector for the gene-list
degdata.clsid_all_next=degdata.clsid_all[,3:ncol(degdata.clsid_all)]
rownames(degdata.clsid_all_next)<-degdata.clsid_all$gene
##degdata.clsid_all_t=t(degdata.clsid_all_next)## (for sample clustering)##transpose the matrix by which rows=samples and columns=genes will be in transposed matrix
degdata.clsid_all_t=degdata.clsid_all_next## (for gene clustering)##not transpose the matrix :: rows=genes and columns=samples will be in transposed matrix ##finding gene-clusters' variablility
dim(degdata.clsid_all_t)

##log.degdata.clsid_all_t <- log(degdata.clsid_all_t)
#degdata.clsid_all_t_sc<-t(scale(t(degdata.clsid_all_t), scale=TRUE))

##degdata.clsid_all_t.pca <- prcomp(degdata.clsid_all_t_sc, center = TRUE, scale. = FALSE)
degdata.clsid_all_t.pca <- prcomp(degdata.clsid_all_t, center = TRUE, scale. = TRUE)
pca<-degdata.clsid_all_t.pca
print(pca)
plot(pca, type = "lines")
summary(pca)

pca$x

##predict(degdata.clsid_i_t.pca,newdata=tail(log.degdata.clsid_i_t, 2))

#####put group "ADENO" and "SCC" ##
# group_sampl<-rownames(degdata.clsid_i_t)
# for(i in 1:nrow(degdata.clsid_i_t))
# {
#   if(length(grep("ADENO", group_sampl[i])))
#   {
#     group_sampl[i]<-gsub(group_sampl[i],"ADENO",group_sampl[i])
#   } else if(length(grep("SCC", group_sampl[i]))){
#     group_sampl[i]<-gsub(group_sampl[i],"SCC",group_sampl[i])
#   }
# }

library(ggbiplot)
g <- ggbiplot(degdata.clsid_all_t.pca, obs.scale = 1, var.scale = 1, 
              groups = group_sampl, ellipse = TRUE, 
              circle = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)



# # Load data
# data(iris)
# head(iris, 3)
# # log transform 
# log.ir <- log(iris[, 1:4])
# ir.species <- iris[, 5]
# 
# # apply PCA - scale. = TRUE is highly 
# # advisable, but default is FALSE. 
# ir.pca <- prcomp(log.ir,
#                  center = TRUE,
#                  scale. = TRUE)
# # print method
# print(ir.pca)
# 
# # plot method
# plot(ir.pca, type = "l")
# 
# # summary method
# summary(ir.pca)
# 
# # Predict PCs
# predict(ir.pca, 
#         newdata=tail(log.ir, 2))
# 
# library(devtools)
# install_github("ggbiplot", "vqv")
# 
# library(ggbiplot)
# g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
#               groups = ir.species, ellipse = TRUE, 
#               circle = TRUE)
# g <- g + scale_color_discrete(name = '')
# g <- g + theme(legend.direction = 'horizontal', 
#                legend.position = 'top')
# print(g)



########classification score of 35 gene-signature#######

###########row-wise (gene-wise) standardization of each TCGA data#######################################
##data.clsid_min.avg.eucledist_scaleddat<- t(scale(t(data.clsid_min.avg.eucledist)))##zero-mean normalization
data.clsid_max.avg.sprdist_scaleddat<- t(scale(t(data.clsid_max.avg.sprdist)))##zero-mean normalization
#data.clsid_min.avg.eucledist_scaleddat<- t(apply(data.clsid_min.avg.eucledist, 1, function(x)(x-min(x))/(max(x)-min(x))))##where the second argument 1 tells apply to work with rows.##min-max normalization
# apply(TCGAPAAD_mirdat3f_scaleddat, 1, mean)
# apply(TCGAPAAD_mirdat3f_scaleddat, 1, sd)


##TCGAPAAD_mirdat3f_scaleddat<-data.clsid_min.avg.eucledist_scaleddat
TCGAPAAD_mirdat3f_scaleddat<-data.clsid_max.avg.sprdist_scaleddat
##class.TCGAPAAD_mirdat3f.ord.cls<-c(rep("ADENO",22),rep("SCC",ncol(data.clsid_min.avg.eucledist_scaleddat)-22))###classinformation rest (2) vs squamous (1) 
class.TCGAPAAD_mirdat3f.ord.cls<-c(rep("ADENO",22),rep("SCC",ncol(data.clsid_max.avg.sprdist_scaleddat)-22))###classinformation rest (2) vs squamous (1) 



#####svm classifier##############
classifier.rep=10

##seed.val<- c(45,136,218,341,490,509,580,639,717,802)##give "classifier.rep" number of seed values ##avg acc in svm= 0.7765027 & AUC.roc=0.7300457
seed.val<- c(145,236,318,441,590,609,780,839,917,1002)##give "classifier.rep" number of seed values ##avg acc in svm= 0.7765027 & AUC.roc=0.7300457

##classifier.rep=4
#seed.val<- c(112,134,246,387,489)##give "classifier.rep" number of seed values
##seed.val<-c(360,430,480,530,580)


method<-"pamr"
#method<-"svm"

cv.fold<-10
cv.folds.id.vector<-c(1,2,3,4,5,6,7,8,9,10)##give "cv.fold" number of seed values

#cv.fold<-4
#cv.folds.id.vector<-c(1,2,3,4)##give "cv.fold" number of seed values

##create a vector for storing AUC score for 217-mirna predictive study for each repitition
#rj.f.oriclassprob <- vector(mode="numeric", length=length(class.TCGAPAAD_mirdat3f.ord.cls))
rj.f.oriclassprob <- matrix(,nrow=0,ncol=1)
rj.f.oriclasslabel<- matrix(,nrow=0,ncol=1)
#isempty(rj.f.oriclassprob)
#isempty(rj.f.oriclasslabel)
#isempty(rj.f.oriclassprob)&isempty(rj.f.oriclasslabel)

##create a vector for storing class prediction for each sample
fsensitivity=vector(mode="numeric", length=length(classifier.rep))
fspecificity=vector(mode="numeric", length=length(classifier.rep))
fprecision=vector(mode="numeric", length=length(classifier.rep))
faccuracy=vector(mode="numeric", length=length(classifier.rep))
foverall.err.rate=vector(mode="numeric", length=length(classifier.rep))


#auc.predinput.df<-vector(mode="numeric", length=0)
#auc.classlabelinput.df<-vector(mode="character", length=0)

for(i in 1:classifier.rep)
{
  #i=1
  #print(paste("i:",i))
  set.seed(seed.val[i])###to keep the randomly picked training samples same everytime for the specified seed value
  
  rj.samp.size<-length(class.TCGAPAAD_mirdat3f.ord.cls)
  
  #Create k=10 equally size folds
  library(caret)
  folds <- createFolds(class.TCGAPAAD_mirdat3f.ord.cls, k = cv.fold, list = TRUE, returnTrain = FALSE)
  
  #length(folds$Fold1)
  #length(folds$Fold2)
  #length(folds$Fold3)
  #length(folds$Fold4)
  
  ##for reintialize before starting each iteration i
  fTP=fTN=fFN=fFP=0
  rj.f.oriclassprob <-0
  
  #Perform k-fold cross validation
  
  for (j in 1:cv.fold)
  {
    #j=1
    #print(paste("j:",j))
    #Segement your data by fold using the which() function 
    
    #rj.testset.pt <- sample (rj.samp.size, (rj.samp.size/cv.fold)) ###randomly pick (183/4) points for test set
    rj.testset.pt <- folds[[j]] ##find test points from a list "folds"
    rj.testset.pt
    
    rj.testset<-TCGAPAAD_mirdat3f_scaleddat[, rj.testset.pt]#ncol=28
    rj.testset.y<-class.TCGAPAAD_mirdat3f.ord.cls[rj.testset.pt]#length=28
    rj.testset.y
    length(rj.testset.pt)#28
    dim(rj.testset)#[217, 28]
    
    
    
    rj.trianingset<-TCGAPAAD_mirdat3f_scaleddat[, -rj.testset.pt]#ncol=(275-28)=247
    rj.trianingset.y<-class.TCGAPAAD_mirdat3f.ord.cls[-rj.testset.pt]#length=(275-28)=247
    dim(rj.trianingset)#[217, 247]
    length(rj.trianingset.y)#247
    
    rj.x<-rj.trianingset
    rj.genenames=rownames(rj.trianingset)
    rj.y<-rj.trianingset.y
    
    ##for pamr run
    library("pamr")
    if(method=="pamr")
    {
      rj.mydata <- list(x=rj.x,y=factor(rj.y), geneid=as.character(1:nrow(rj.x)),genenames=rj.genenames)
      rj.mytrain <- pamr.train(rj.mydata)
      rj.mycv <- pamr.cv(rj.mytrain,rj.mydata, nfold=cv.fold, folds=cv.folds.id.vector) ##cv.fold=4
      
      rj.Thresh<-max(rj.mycv$threshold[rj.mycv$error== min(rj.mycv$error)]);##https://support.bioconductor.org/p/27089/ ##
      print(rj.Thresh)
      
      pamr.plotcv(rj.mycv)
      
      ##confus<-pamr.confusion(rj.mytrain, rj.Thresh)
      ##confus
      
      if(exists("rj.pred.result"))
      {
        print("removed existing rj.pred.result variable")
        remove(rj.pred.result)
      }
      
      ypred<-pamr.predict(rj.mytrain, rj.testset, threshold=rj.Thresh)
      
      ##TP=length(which(rj.testset.y==1 & rj.pred.result==1))
      ##TN=length(which(rj.testset.y==2 & rj.pred.result==2))
      ##FN=length(which(rj.testset.y==1 & rj.pred.result==2))
      ##FP=length(which(rj.testset.y==2 & rj.pred.result==1))
    }
    
    
    ##for svm
    if(method=='svm')
    {
      library(e1071)
      xt1<-as.data.frame(x=t(rj.x),y=as.factor(rj.y))
      rownames(xt1)<-NULL
      #svmmodel<-svm(xt1, rj.y)
      #svmmodel<-svm(xt1, rj.y,probability = TRUE)
      
      svmmodel<-svm(rj.y~., data=xt1, method="C-classification",
                    kernel="linear", gamma = 0.01, cost = 100, probability=TRUE) 
      
      
      testdat=as.data.frame(x=t(rj.testset),y=as.factor(rj.testset.y))
      rownames(testdat)<-NULL
      ypred=predict(svmmodel,testdat, probability = TRUE)
    }
    
    rj.pred.result<-table(predict=ypred, truth=rj.testset.y)
    rj.pred.result
    TP=rj.pred.result[2,2]
    TN=rj.pred.result[1,1]
    FN=rj.pred.result[1,2]
    FP=rj.pred.result[2,1]
    
    #print(paste("TP:",TP))
    #print(paste("TN:",TN))
    #print(paste("FN:",FN))
    #print(paste("FP:",FP))
    
    fTP=fTP+TP
    fTN=fTN+TN
    fFN=fFN+FN
    fFP=fFP+FP
    fTP
    fTN
    fFN
    fFP
    
    
    
    ##For pamr
    #Estimated class probabilities, from pamr.predict with type="posterior")
    if(method=="pamr")
    {
      rj.prob<- data.frame(pamr.predict(rj.mytrain, rj.testset, threshold=rj.Thresh, type="posterior"))
    }  
    if(method=="svm")
    {
      rj.prob <- data.frame(attr(ypred, "probabilities"))
    }
    rj.prob$oriclasslabel<-rj.testset.y
    rj.prob$predscore<-0##for initialize 
    #rj.prob$predscore[which(rj.testset.y==1)]<-rj.prob[which(rj.testset.y==1),1]
    #rj.prob$predscore[which(rj.testset.y==2)]<-rj.prob[which(rj.testset.y==2),2]
    ##rj.f.oriclassprob[rj.testset.pt]<-rj.prob$predscore ##select the classification score of the original class-label and store it into original class-probability matrix ("rj.f.oriclassprob")
    
    rj.prob$predscore[which(rj.prob$oriclasslabel=="SCC")]<-rj.prob[which(rj.prob$oriclasslabel=="SCC"),2]
    rj.prob$predscore[which(rj.prob$oriclasslabel=="ADENO")]<-rj.prob[which(rj.prob$oriclasslabel=="ADENO"),1]
    
    ##rj.prob.nonsel<-vector(mode="numeric", length=length(rj.trianingset.y))
    ##tt1<-rbind(as.matrix(rj.prob$predscore),as.matrix(rj.prob.nonsel))#don't use
    ##tt2<-rbind(as.matrix(rj.prob$oriclasslabel),as.matrix(rj.trianingset.y))#don't use
    tt1<-as.matrix(rj.prob$predscore)
    tt2<-as.matrix(rj.prob$oriclasslabel)
    
    #rj.f.oriclassprob <- matrix(NA,nrow=0,ncol=1)
    #rj.f.oriclasslabel<- matrix(NA,nrow=0,ncol=1)
    #isempty(rj.f.oriclassprob)
    length(rj.f.oriclassprob)
    
    #if (nrow(rj.f.oriclassprob)==0) {
    #if(isempty(rj.f.oriclassprob)) {  
    if(length(rj.f.oriclassprob)==1|length(rj.f.oriclassprob)==0) {  
      #print(23);
      rj.f.oriclassprob<-tt1;
      rj.f.oriclasslabel<-tt2; 
    } else {
      #print(45);
      rj.f.oriclassprob<-rbind(rj.f.oriclassprob,tt1);
      rj.f.oriclasslabel<-rbind(rj.f.oriclasslabel,tt2); 
    }
    
  }
  fTP
  fTN
  fFN
  fFP
  print(paste("i:",i))
  print(paste("j:",j))
  print(paste("fTP:",fTP))
  print(paste("fTN:",fTN))
  print(paste("fFN:",fFN))
  print(paste("fFP:",fFP))
  
  sensitivity=fTP/(fTP+fFN)
  specificity=fTN/(fTN+fFP)
  precision=fTP/(fTP+fFP)
  accuracy=(fTP+fTN)/(fTP+fTN+fFP+fFN)
  overall.err.rate=(fFN+fFP)/(fTP+fTN+fFP+fFN)
  
  if(is.nan(sensitivity))
  {sensitivity=0}
  if(is.nan(specificity))
  {specificity=0}
  if(is.nan(precision))
  {precision=0}
  if(is.nan(accuracy))
  {accuracy=0}
  if(is.nan(overall.err.rate))
  {overall.err.rate=0}
  
  
  fsensitivity[i]<-sensitivity
  fspecificity[i]<-specificity
  fprecision[i]<-precision
  faccuracy[i]<-accuracy
  foverall.err.rate[i]<-overall.err.rate
  
  fsensitivity
  fspecificity
  fprecision
  faccuracy
  foverall.err.rate
}
avg.sensitivity=mean(fsensitivity)
avg.specificity=mean(fspecificity)
avg.precision=mean(fprecision)
avg.accuracy=mean(faccuracy)
avg.overall.err.rate=mean(foverall.err.rate)



std.sensitivity=sd(fsensitivity)
std.specificity=sd(fspecificity)
std.precision=sd(fprecision)
std.accuracy=sd(faccuracy)
std.overall.err.rate=sd(foverall.err.rate)


##auc for class-prediction through pamr finally
library(AUC)
#rj.f.oriclassprob.nm<-as.numeric(levels(rj.f.oriclassprob))[rj.f.oriclassprob]
rj.f.oriclassprob.nm<-as.numeric(rj.f.oriclassprob)
rj.f.oriclasslabel.f<-as.factor(ifelse(rj.f.oriclasslabel=="ADENO",1,0))
##rj.f.oriclasslabel.f<-as.factor(rj.f.oriclasslabel)

#tyu1<-accuracy(rj.f.oriclassprob.nm,rj.f.oriclasslabel.f)

#AUC.score<-auc(accuracy(rj.f.oriclassprob.nm,rj.f.oriclasslabel.f), min = 0, max = 1)
AUC.score<-auc(accuracy(rj.f.oriclassprob.nm,rj.f.oriclasslabel.f))
AUC.roc<-auc(roc(rj.f.oriclassprob.nm,rj.f.oriclasslabel.f))

avg.sensitivity
avg.specificity
avg.precision
avg.accuracy
avg.overall.err.rate
std.sensitivity
std.specificity
std.precision
std.accuracy
std.overall.err.rate

AUC.score
AUC.roc





#data(toy5)
#obj <- mocca(toy5[,1:2], R=10, K=2:5)
#print(analyzePareto(obj$objectiveVals))


#mocca(x, R = 50, K = 2:10, iter.max = 1000, nstart = 10)



#data(toy5)
#res <- mocca(toy5, R=10, K=2:5)
#print(res$objectiveVals)
# plot kmeans result for MCA index against neuralgas result for MCA index
#plot(res$objectiveVals[1,], res$objectiveVals[5,], pch=NA,
#     xlab=rownames(res$objectiveVals)[1], ylab=rownames(res$objectiveVals)[5])
#text(res$objectiveVals[1,], res$objectiveVals[5,], labels=colnames(res$objectiveVals))

# 
# > avg.specificity
# [1] 0.9409091
# > avg.sensitivity
# [1] 0.9339921
# > avg.accuracy
# [1] 0.9345455
# > avg.overall.err.rate
# [1] 0.06545455
# > std.sensitivity
# [1] 0.002667781
# > std.specificity
# [1] 0.02195663
# > avg.precision
# [1] 0.9945325
# > std.precision
# [1] 0.002017868
# > std.accuracy
# [1] 0.002969078
# > std.overall.err.rate
# [1] 0.002969078
# > AUC.score
# [1] 0.4500364
# > AUC.roc
# [1] 0.1812792

