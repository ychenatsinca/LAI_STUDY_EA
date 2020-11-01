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
fun_table.qc <- function(wrk.table, qc1.set=0.5, qc2.set=0.2, area.set=1000) {
 
      #subset the datset with QC and QA score
      #Avariable area
      wrk.table <- subset(wrk.table, ((wrk.table$eve.for.pix.aff)> area.set)  )
      #QC1 mean difference  on LAI #QC2 measurement uncertainty on LAI  
      #qc1.set <- 0.5
      #qc2.set <- 0.2
      wrk.table <- subset(wrk.table, (abs(wrk.table$qc1.score)<=qc1.set &  abs(wrk.table$qc2.score)>=qc2.set ))
      return(wrk.table)
      } 
####################################

#set qc/qa parameters
qc1.set <- 0.5
qc2.set <- 0.18
area.set <- 100

#get file-list 
in.files <- list.files(path="./Rda_EA/", pattern="EA_*")
   
table.all <- NA
for (id in 1:length(in.files)  ) { 
#for (id in 1:1) { 
   load(paste("./Rda_EA/",in.files[id],sep=""))
   table.comb$date.day <- substr(table.comb$date.time,start=7,stop=8)
   
   # subset datataset (this is super critical!! we develop this QC flags)
   table.comb <- fun_table.qc(wrk.table=table.comb, qc1.set=qc1.set, qc2.set=qc2.set, area.set=area.set)
   # combine all table
   table.all <- rbind(table.comb,table.all) 
 
   #scatter.smooth(y=table.comb$eff.size.for, x=(table.comb$date.mon+ as.integer(table.comb$date.day)/30 ),
   #               col="lightgray", span=0.5, degree=1, 
   #               lpars=list(lwd=3,col="gray"), ylim=c(-0.5,0.5), xlab="",ylab="Effect size",  xaxt = "n")
}

# remove the NA cases
table.all <- table.all[complete.cases(table.all), ]


library("dplyr")
library("randomForest")
library("caret")
library("e1071")


# read the enso index table 
enso.table <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/ENSO_DATA/NOI_1999to2019.txt",sep=",",header=T)

# assign the positve/negtive value of eff_size as group B/A
table.all$group <- NA
table.all$group[table.all$eff.size.for >0] <- "Growth"
table.all$group[table.all$eff.size.for <0] <- "Damage"
table.all$group <- as.factor(table.all$group)

#assign the  index to the TC event 
table.all$enso <- NA
#assign the  index to the TC event 
for (it in 1: length(table.all$group))  {
   # find the correct time step for all events 
   table.all$enso[it] <- enso.table$Ocean_Nino_Index[which(enso.table$Year == table.all$date.yr[it]  
                                                         & enso.table$Mon == table.all$date.mon[it] ) ] 
}


mat1 <-matrix(c(1,1,2,3),2,2)
layout(mat1, widths=c(2,2),height=c(2,2))



hist.breaks <- c(seq(-5,5.0,length=100))
hist.col <- ifelse( hist.breaks < 0, "brown", "forestgreen") 
 

hist(table.all$eff.size.for, col=hist.col
     ,breaks=hist.breaks, main="All events", xlim=c(-1,1) )
hist(table.all$eff.size.for[ table.all$enso >= 0.5 ], col= hist.col,
     ,breaks=hist.breaks, main="El Nino (positive phase events)",xlim=c(-1,1))
hist(table.all$eff.size.for[ table.all$enso <= -0.5 ], col= hist.col,
     breaks=hist.breaks, main="La Nina (negtive phase events)",xlim=c(-1,1))


table.ana <- select(table.all , group, bef.for.lai.aff, bef.for.asw.aff, eve.lat.med.aff, eve.for.mws.aff, enso)

# Calculate the size of each of the data sets:
data_set_size <- floor(nrow(table.ana))
# Mote Carlo sampling with n/2 of total data_set_size
#data_set_size <- as.integer(data_set_size/2)
#data_set_size = 131 



ld.go <- T

if (ld.go) {
# subset table
imp.table <- NA 

for (i in 1:200) {
# Set random seed to make results reproducible:
set.seed(17*i)
# Generate a random sample of "data_set_size" indexes
indexes <- sample(1:nrow(table.ana), size = data_set_size)


# Assign the data to the correct sets
training <- table.ana[indexes,]


#Apply Classification with Random Forest:
# In the randomForest Package this would mean that you make a classification (type = "classification"), ntree = 2000 means that 2000 prediction trees are grown

rf_class = randomForest(group ~ ., data=training, ntree=200, mtry=2,norm.votes=FALSE,  importance=TRUE,proximity=TRUE, type="classification", na.action=na.omit)
#
# Set your workingdirectory for the output here
#pdf(file='all_trees.pdf', paper = "a4r", width = 10) #Name of the PDF File. In this case it will be a PDF with 2000 pages
#for(k in 100:110){
  # the getTree function is part of the randomForest package; it extracts the structure of a tree from a randomForest object - I hope the biomod2 package has the same functionality
#  tree <- getTree(rf_class, k=k, labelVar = TRUE)
#  str(tree)
#  print(sum(tree[,1] +tree[,2]))
#  d <- to.dendrogram(tree)
#  str(d)
#  plot(d,center=TRUE,edgePar=list(t.cex=.55,p.col=NA,p.lty=0), yaxt = "n", 
#  main = paste("Tree Numbe:",k," Accuracy of classification: ",round(1-rf_class$err.rate[k],digits = 2)*100,"%"))
  #dev.new()
#}
#dev.off()

#dev.new()
#plot importance 
#varImpPlot(rf_class)


#export the importance from the model 

tmp.imp.table <- as.data.frame(importance(rf_class)) 
tmp.imp.table$vars <- row.names(tmp.imp.table)

# combine table 
imp.table <- rbind(tmp.imp.table, imp.table, make.row.names=FALSE)

 

}


# rm NA cases
imp.table <- imp.table[complete.cases(imp.table), ]

# set vars as factor 
imp.table$vars <- as.factor(imp.table$vars)

# adjust the orders 
imp.table$vars <- factor(imp.table$vars , levels=c("eve.for.mws.aff","bef.for.lai.aff", "bef.for.asw.aff","enso","eve.lat.med.aff"))

boxplot( imp.table$MeanDecreaseAccuracy ~ imp.table$vars,
         medcol="red",col="lightblue", horizontal=TRUE)


#import the package
#library("randomForest")
#library("reprtree")

# Perform training:
#rf_class = randomForest(group ~ ., data=training, ntree=100, mtry=2, importance=TRUE)
# Plot the tree 
#reprtree:::plot.getTree(rf_class)


#test <- training
#library(ggplot2)
#pGroup <- predict(rf_class,test,'vote')[,2]
#plotData <- lapply(names(test[,2:5]), function(x){
#  out <- data.frame(
#    var = x,
#    type = c(rep('Actual',nrow(test)),rep('Predicted',nrow(test))),
#    value = c(test[,x],test[,x]),
#    group = c(as.numeric(test$group)-1,pGroup)
#    )
#  out$value <- out$value-min(out$value) #Normalize to [0,1]
#  out$value <- out$value/max(out$value)
#  out
#})
#plotData <- do.call(rbind,plotData)
##qplot(value, group, data=plotData, facets = type ~ var, geom='smooth', span = 0.5)


} #ld_go

dev.new()

# performing regrssion trees analysis based on the random forest analysis (assign the importance variables) 
#  
library(rpart)       # performing regression trees
library(rpart.plot)  # plotting regression trees

library(rsample)     # data splitting 
library(dplyr)       # data wrangling


table.ana$group <- as.character(table.ana$group)
table.ana$group <- as.factor(table.ana$group)

# build the regression tree model
#tc.m1 <- rpart(formula = group ~.,  data = table.ana  )     
#rpart.plot (tc.m1)





