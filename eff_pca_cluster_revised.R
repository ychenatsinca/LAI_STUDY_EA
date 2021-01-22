## FUNCTIONS

#################################################################
#**************************
#return the rules of a tree
#**************************
getConds<-function(tree){
  #store all conditions into a list
  conds<-list()
  #start by the terminal nodes and find previous conditions
  id.leafs<-which(tree$status==-1)
  j<-0
  for(i in id.leafs){
    j<-j+1
    prevConds<-prevCond(tree,i)
    conds[[j]]<-prevConds$cond
    while(prevConds$id>1){
      prevConds<-prevCond(tree,prevConds$id)
      conds[[j]]<-paste(conds[[j]]," & ",prevConds$cond)
      if(prevConds$id==1){
        conds[[j]]<-paste(conds[[j]]," => ",tree$prediction[i])
        break()
      }
    }
    
  }
  
  return(conds)
}

#**************************
#find the previous conditions in the tree
#**************************
prevCond<-function(tree,i){
  if(i %in% tree$right_daughter){
    id<-which(tree$right_daughter==i)
    cond<-paste(tree$split_var[id],">",tree$split_point[id])
  }
  if(i %in% tree$left_daughter){
    id<-which(tree$left_daughter==i)
    cond<-paste(tree$split_var[id],"<",tree$split_point[id])
  }
  
  return(list(cond=cond,id=id))
}

#remove spaces in a word
collapse<-function(x){
  x<-sub(" ","_",x)
  
  return(x)
}


#####################################################################################
to.dendrogram <- function(dfrep,rownum=1,height.increment=0.1){
  
  if(dfrep[rownum,'status'] == -1){
    rval <- list()
    
    attr(rval,"members") <- 1
    attr(rval,"height") <- 0.0
    attr(rval,"label") <- dfrep[rownum,'prediction']
    attr(rval,"leaf") <- TRUE
    
  }else{##note the change "to.dendrogram" and not "to.dendogram"
    left <- to.dendrogram(dfrep,dfrep[rownum,'left daughter'],height.increment)
    right <- to.dendrogram(dfrep,dfrep[rownum,'right daughter'],height.increment)
    rval <- list(left,right)
    
    attr(rval,"members") <- attr(left,"members") + attr(right,"members")
    attr(rval,"height") <- max(attr(left,"height"),attr(right,"height")) + height.increment
    attr(rval,"leaf") <- FALSE
    attr(rval,"edgetext") <- paste(dfrep[rownum,'split var'],"\n<",round(dfrep[rownum,'split point'], digits = 2),"=>", sep = " ")
  }
  
  class(rval) <- "dendrogram"
  
  return(rval)
}



#####################################
# Quality Check funtion 
#####################################    
#QC1 mean difference  on LAI
#QC2 measurement uncertainty on LAI  
fun_table.qc <- function(wrk.table, qc1.set=0.5, qc2.set=0.2, area.set=0) {
 
      #print( paste("Total events:", length(wrk.table$eve.for.pix.aff),sep=""))
      #subset the datset with QC and QA score
      #Avariable area
      wrk.table <- subset(wrk.table, ((wrk.table$eve.for.pix.aff)> area.set)  )
      all.eve <- length(wrk.table$eve.for.pix.aff)
      print( paste("Total events:", all.eve,sep=""))
 
      #QC1 mean difference on LAI #QC2 measurement uncertainty on LAI  
      #qc1.set <- 0.5
      #qc2.set <- 0.2
      flag1.table <- subset(wrk.table, (abs(wrk.table$qc1.score)<=qc1.set) )
      tmp.eve <- length(flag1.table$eve.for.pix.aff)
      flag1.eve <- all.eve - length(flag1.table$eve.for.pix.aff)
      print( paste("Flag-1 events (remove Big LAI difference):", flag1.eve ,sep=""))
      
      flag2.table <- subset(wrk.table, ((abs(wrk.table$qc2.score)<qc2.set) & (abs(wrk.table$qc1.score)<qc1.set) ))
      flag2.eve <- tmp.eve - length(flag2.table$eve.for.pix.aff) 
      print( paste("Flag-2 events (Neutral events):", flag2.eve,sep=""))

      wrk.table <-  subset(wrk.table, ((abs(wrk.table$qc2.score)>qc2.set) & (abs(wrk.table$qc1.score)<qc1.set) ))
      return(list(qcqa=wrk.table,flag1=flag1.table,flag2=flag2.table))
      } 
####################################

#set qc/qa parameters
qc1.set <- 0.6
qc2.set <- 0.2
area.set <- 0

#get file-list 
in.files <- list.files(path="./Rda_60/", pattern="EA_*")
   
table.all <- NA
for (id in 1:length(in.files)  ) { 
#for (id in 1:1) { 
   load(paste("./Rda_60/",in.files[id],sep=""))
   table.comb$date.day <- substr(table.comb$date.time,start=7,stop=8)
   
  # combine all table
   table.all <- rbind(table.comb,table.all) 
 
   #scatter.smooth(y=table.comb$eff.size.for, x=(table.comb$date.mon+ as.integer(table.comb$date.day)/30 ),
   #               col="lightgray", span=0.5, degree=1, 
   #               lpars=list(lwd=3,col="gray"), ylim=c(-0.5,0.5), xlab="",ylab="Effect size",  xaxt = "n")
}

# remove the NA cases
table.all <- table.all[complete.cases(table.all), ]
raw.table <- table.all

table.all <- subset(table.all, table.all$eve.for.mws.aff>1.0)

# subset datataset (this is super critical!! we develop this QC flags)
# combine qcqc table and flag2 table (neutral)  
table.tmp <- fun_table.qc(wrk.table=table.all, qc1.set=qc1.set, qc2.set=qc2.set, area.set=area.set)

# assign the positve/negtive value of eff_size as group B/A
table.tmp[['qcqa']]$group <- NA
table.tmp[['qcqa']]$group[table.tmp[['qcqa']]$eff.size.for >  0.2] <- "Postive"
table.tmp[['qcqa']]$group[table.tmp[['qcqa']]$eff.size.for < -0.2] <- "Negative"
table.tmp[['qcqa']]$group[abs(table.tmp[['qcqa']]$eff.size.for) <= 0.2] <- "Neutral"

table.tmp[['flag2']]$group <- "Neutral"

table.all <- rbind(table.tmp[['qcqa']],table.tmp[['flag2']])  
#table.all <- table.tmp[['qcqa']]


table.all$group <- as.factor(table.all$group)


library("dplyr")
library("randomForest")
library("caret")
library("e1071")


# read the Oceanic Nino Index table 
oni.enso.table <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/ENSO_DATA/NOI_1999to2019.txt",sep=",",header=T)


#assign the  index to the TC event 
table.all$oni.enso <- NA
#assign the  index to the TC event 
for (it in 1: length(table.all$group))  {
   # find the correct time step for all events 
   table.all$oni.enso[it] <- oni.enso.table$Ocean_Nino_Index[which(oni.enso.table$Year == table.all$date.yr[it]  
                                                         & oni.enso.table$Mon == table.all$date.mon[it] ) ] 
}


mat1 <-matrix(c(1,1,2,3),2,2)
layout(mat1, widths=c(2,2),height=c(2,2))



hist.breaks <- c(seq(-5,5.0,length=100))
hist.col <- ifelse( hist.breaks < 0, "orange2", "forestgreen") 
 

hist(table.all$eff.size.for, col=hist.col
     ,breaks=hist.breaks, main="All events", xlim=c(-1,1) )
hist(table.all$eff.size.for[ table.all$oni.enso >= 0.5 ], col= hist.col,
     ,breaks=hist.breaks, main="El Nino (positive phase events)",xlim=c(-1,1))
hist(table.all$eff.size.for[ table.all$oni.enso <= -0.5 ], col= hist.col,
     breaks=hist.breaks, main="La Nina (negtive phase events)",xlim=c(-1,1))



eff_cex  <-  sqrt(table.all$eve.for.pix.aff/3.1416)/100
table.all$eff.cex <- eff_cex

table.all$date.mon.day <- (as.numeric(table.all$date.day)/10 + as.numeric(table.all$date.mon))


# prepare table for pca analysis
table.pca <- table.all[, c("eff.size.for",
                         "date.mon.day","eve.for.mws.aff","eve.for.acf.aff","bef.for.acf.aff","bef.for.lai.aff","bef.for.asw.aff","eve.lat.med.aff","eff.cex","oni.enso")]

go_pca <- FALSE

if (go_pca) {
# funtion for pair analysis
upper.panel <- function(x, y){
                points(x,y, cex=table.comb$eff.cex)
              }
# pair analysis
#pairs(   ~eff.size.for+eve.for.mws.aff+eve.for.acf.aff+bef.for.lai.aff+bef.for.asw.aff+eve.lat.med.aff+enso,upper.panel=upper.panel,
#lower.panel=NULL, data=table.all, main="Pair Correlation")


# pca analysis
pca <- prcomp(formula = ~ date.mon.day+eve.for.mws.aff+eve.for.acf.aff+bef.for.acf.aff+bef.for.lai.aff+bef.for.asw.aff+eve.lat.med.aff+eff.cex+oni.enso,
              data = table.pca,
              scale = TRUE)
#show pca
pca
dev.new()
layout(matrix(c(1,2), 1, 2,byrow = TRUE))
#plot pca
plot(pca,
     type="line", #
     main="Eigenvalue Plot of PCA") #

#show variance =1.0
abline(h=1, col="blue") # Kaiser eigenvalue-greater-than-one rule

vars <- (pca$sdev)^2  # 從pca中取出標準差(pca$sdev)後再平方，計算variance(特徵值)
vars

# 計算每個主成分的解釋比例 = 各個主成分的特徵值/總特徵值
props <- vars / sum(vars)
props

cumulative.props <- cumsum(props)  # 累加前n個元素的值
cumulative.props

#show top three PCs
cumulative.props[6]

# 累積解釋比例圖
plot(cumulative.props)

# pca$rotation
pca.data <- pca$x[, 1:6]
pca.data

# 特徵向量(原變數的線性組合)
pca.eigenvector <- pca$rotation[, 1:6]
pca.eigenvector
pca$rotation

first.pca  <- pca.eigenvector[, 1]   #  
second.pca <- pca.eigenvector[, 2]   #  
third.pca  <- pca.eigenvector[, 3]   #  
fourth.pca <- pca.eigenvector[, 4]   # 
fiveth.pca <- pca.eigenvector[, 5]   #  
sixth.pca  <- pca.eigenvector[, 6]   #  


dev.new()
layout(matrix(c(1,2,3,4,5,6), 3, 2,byrow = TRUE))
#par(mfrow=c(2,3))
#
# 第一主成份：由小到大排序原變數的係數
first.pca[order(first.pca, decreasing=FALSE)]
# 使用dotchart，繪製主成份負荷圖
dotchart(first.pca[order(first.pca, decreasing=FALSE)] ,   # 排序後的係數
         main="Loading Plot for PC1",                      # 主標題
         xlab="Variable Loadings",                         # x軸的標題
         col="red")
abline(v=0, col="gray",lty="dashed")

# 第二主成份：由小到大排序原變數的係數
second.pca[order(second.pca, decreasing=FALSE)]
# 使用dotchart，繪製主成份負荷圖
dotchart(second.pca[order(second.pca, decreasing=FALSE)] ,  # 排序後的係數
         main="Loading Plot for PC2",                       # 主標題
         xlab="Variable Loadings",                          # x軸的標題
         col="blue")                                        # 顏色
abline(v=0, col="gray",lty="dashed")

# 第三主成份：由小到大排序原變數的係數
third.pca[order(third.pca, decreasing=FALSE)]
# 使用dotchart，繪製主成份負荷圖
dotchart(third.pca[order(third.pca, decreasing=FALSE)] ,   # 排序後的係數
         main="Loading Plot for PC3",                      # 主標題
         xlab="Variable Loadings",                         # x軸的標題
         col="purple")                                     # 顏色
abline(v=0, col="gray",lty="dashed")

# 第四主成份：由小到大排序原變數的係數
fourth.pca[order(fourth.pca, decreasing=FALSE)]
# 使用dotchart，繪製主成份負荷圖
dotchart(fourth.pca[order(fourth.pca, decreasing=FALSE)] , # 排序後的係數
         main="Loading Plot for PC4",                      # 主標題
         xlab="Variable Loadings",                         # x軸的標題
         col="forestgreen")                                # 顏色
abline(v=0, col="gray",lty="dashed")

# 第五主成份：由小到大排序原變數的係數
fiveth.pca[order(fiveth.pca, decreasing=FALSE)]
# 使用dotchart，繪製主成份負荷圖
dotchart(fiveth.pca[order(fiveth.pca, decreasing=FALSE)] , # 排序後的係數
         main="Loading Plot for PC5",                      # 主標題
         xlab="Variable Loadings",                         # x軸的標題
         col="orange")                                     # 顏色
abline(v=0, col="gray",lty="dashed")

# 第六主成份：由小到大排序原變數的係數
sixth.pca[order(sixth.pca, decreasing=FALSE)]
# 使用dotchart，繪製主成份負荷圖
dotchart(sixth.pca[order(fiveth.pca, decreasing=FALSE)] ,  # 排序後的係數
         main="Loading Plot for PC6",                      # 主標題
         xlab="Variable Loadings",                         # x軸的標題
         col="brown")                                      # 顏色
abline(v=0, col="gray",lty="dashed")

#assign the PCs to the table.comb
table.all$pc1<-pca.data[,1]
table.all$pc2<-pca.data[,2]
table.all$pc3<-pca.data[,3]
table.all$pc4<-pca.data[,4]
table.all$pc5<-pca.data[,5]
table.all$pc6<-pca.data[,6]

#do the pair analysis again
dev.new()
# pair analysis
pairs(~eff.size.for+pc1+pc2+pc3+pc4+pc5+pc6,upper.panel=upper.panel,
,lower.panel=NULL, data=table.all, main="Pair Correlation")

}

# go cluster analysis
go.cluster <- FALSE
if (go.cluster) {
library(cluster)
data <- table.all[,c("pc1","pc2","pc3","pc4","pc5","pc6")] 
# pam = Partitioning Around Medoids
kmeans.cluster <- kmeans(data, centers=3)
#pam.cluster <- pam(data, k=3)
# comparison 
print(table(kmeans.cluster$cluster, table.all$group) )
}


# go for factor analysis
go.factor <- TRUE

if (go.factor) {

library(psych)


# Maximum Likelihood Factor Analysis
# entering raw data and extracting 3 factors,
# with varimax rotation

#prepare table for factor analysis
#table.fta <- table.all[, c("eff.size.for",
#                         "date.mon.day","eve.for.mws.aff","eve.for.acf.aff","bef.for.acf.aff","bef.for.lai.aff","bef.for.asw.aff","eve.lat.med.aff","eff.cex","oni.enso")]

# group the data for Tropical Cyclone Characteristics
table.tcc <- table.all[, c("eve.for.mws.aff","eve.for.acf.aff","eve.lat.med.aff","eff.cex","eve.intensity")]


# group the data for Pre-event Characteristics
table.pc <- table.all[, c("date.mon.day","bef.for.acf.aff","bef.for.lai.aff","bef.for.asw.aff","oni.enso")]

fta.tcc <- principal(table.tcc, nfactors=3, rotate="varimax")
print(fta.tcc)
#fta.tcc$scores[,"RC1"]

fta.pc <- principal(table.pc, nfactors=3, rotate="varimax")
print(fta.pc)

#fta.tcc <- principal(table.tcc, cor=TRUE)
#print(fta.tcc, digits=2, cutoff=.3, sort=TRUE)

#fta.pc <- principal(table.pc, cor=TRUE)
#print(fta.pc, digits=2, cutoff=.3, sort=TRUE)






} 





#go classification of gression tree
go.tree <- TRUE

if (go.tree) {
dev.new()

# Classification Tree with rpart
library(rpart)
library(rpart.plot)
go.pdf = FALSE
if(go.pdf) {pdf(file="Regression_Tree.pdf",width=8,height=6) } else{ print("go x-windows") }
#eve.for.mws.aff","eve.for.acf.aff","bef.for.acf.aff","bef.for.lai.aff","bef.for.asw.aff","date.mon","eff.cex","oni.enso"
# Six most important variables are 
# eve.for.acf.aff, date.mon.day, eve.lat.med.aff, bef.for.lai.aff, bef.for.acf.aff, oni.enso


#parameter seting 
rpart.control(minsplit = 30, minbucket = 20, maxdepth = 5)

#fit.0 <- rpart(group ~ eve.for.acf.aff+
#                       bef.for.lai.aff+bef.for.asw.aff+
#                       eve.lat.med.aff+eve.intensity+
#                       oni.enso, method="class", data=table.all)

fit.0 <- rpart(group ~  eve.for.acf.aff + date.mon.day +  eve.lat.med.aff +
                        bef.for.lai.aff + bef.for.acf.aff + oni.enso,
                        method="class", data=table.all)

#pruning the modeled tree05fit.1 = prune(fit.0, cp = 0.1)
fit.1 = prune(fit.0, cp = 0.03)
rpart.plot(fit.1, uniform=TRUE, box.palette=list("lightgray","forestgreen"))


#rpart.plot(fit.0, uniform=TRUE, box.palette=list("orange2","lightgray","forestgreen"))


#


# grow tree
#fit <- rpart(group ~ pc1 + pc2 + pc3 + pc4,
#   method="class", data=table.all)
# plot tree
#rpart.plot(fit.0, uniform=TRUE,
#         main="Classification Tree for Effect Size")
#text(fit, use.n=TRUE, all=TRUE, cex=.8)

# create attractive postscript plot of tree
#post(fit, file = "rgression_tree.ps",
#   title = "Classification Tree for Effect Size")

if(go.pdf) {dev.off()}
#printcp(fit) # display the results
#dev.new()
#plotcp(fit.0) # visualize cross-validation results
#summary(fit) # detailed summary of splits
}












