require(kohonen)
require(lsa)
require(cluster)
require(CoImp)
require(VIM)
require(mice)
require(Metrics)
require(caret)
source("buildMapper.R")
source("genData.R")
source("generateNewData.R")
source("cluster_cutoff_at_first_empty_bin.R")

beta = c(0,0.5,0.5,0.5,1,1,1)
mapper = matrix(0,nrow = 500,ncol = 7) #Coefficient estimates for each iteration
mapperInt = matrix(0,nrow = 500,ncol = 7) #Percent confidence interval coverage of each coefficient 
mapperRMSE = matrix(0,nrow = 500,ncol = 7) #RMSE of coefficients 
mapperBias = matrix(0,nrow = 500,ncol = 7) #Bias of coefficients

dt = matrix(0,nrow = 500,ncol = 7)
dtInt = matrix(0,nrow = 500,ncol = 7)
dtRMSE = matrix(0,nrow = 500,ncol = 7)
dtBias = matrix(0,nrow = 500,ncol = 7)

rf = matrix(0,nrow = 500,ncol = 7)
rfInt = matrix(0,nrow = 500,ncol = 7)
rfRMSE = matrix(0,nrow = 500,ncol = 7)
rfBias = matrix(0,nrow = 500,ncol = 7)

ARMSE = matrix(0,nrow = 500,ncol = 3)

simData = genData()

for (i in 1:500) {
print(i)

pointCloud = simData
missingCol = MCAR(as.matrix(pointCloud[,1:6]),perc.miss = 0.25,setseed = i)@db.missing
pointCloud[,1:6] = missingCol

newData = generateNewData(i)

dtMi = mice(pointCloud,method = "cart",printFlag = FALSE,m = 20)
rfMi = mice(pointCloud,method = "rf",printFlag = FALSE,m = 20)

missingFeatures = which(colSums(is.na(pointCloud)) > 0)
categoricalVars = NULL

cfMapper = numeric(7)
cfDT = matrix(0,nrow = 10,ncol = 7)
cfRF = matrix(0,nrow = 20,ncol = 7)

withinCIMAPPER =  numeric(7)
withinCIDT = matrix(0,nrow = 10,ncol = 7)
withinCIRF = matrix(0,nrow = 10,ncol = 7)

ptm = proc.time()

#Begin mapper imputation

  pointCloudComplete = pointCloud[complete.cases(pointCloud),]
  
  bins = buildMapper(pointCloudComplete)
  pointCloudComplete = as.matrix(pointCloudComplete)
  pointCloud = as.matrix(pointCloud)
  
  codes = list()
  maps = list()
  pcComplete = pointCloud
  
  #Train self-organizing maps on mapper bins
  for (k in 1:length(bins)) {
      points = bins[[k]]
      if(length(points) == 1){
        maps[[k]] = NULL
        next
      }
      mygrid = somgrid(floor(sqrt(length(points))),floor(sqrt(length(points))))
      maps[k] = list(supersom(pointCloudComplete[points,],mygrid,rlen = 10000,maxNA.fraction = 1)) 
      codes[[k]] = as.data.frame(maps[[k]]$codes)
    }
  
  #Construct matrix of ranges
  ranges = matrix(0,nrow = length(bins),ncol = dim(pointCloud)[2] * 2 )
  
  for (r in 1:nrow(ranges)){
    k = 1
    for(j in seq(1,ncol(ranges),2)){
      missingIndices = which(is.na(pointCloud[,k]))
      if(length(bins[[r]][-which(!is.na(match(bins[[r]],missingIndices)))]) > 0){
        pts = bins[[r]][-which(!is.na(match(bins[[r]],missingIndices)))]
      }
      ranges[r,j] = range(pointCloudComplete[pts,k])[1]
      ranges[r,j + 1] = range(pointCloudComplete[pts,k])[2]
      k = k + 1
    }
  }

  
  
  missingIndices = which(!complete.cases(pointCloud))
  
  #Begin imputation
  
  for (mis in missingIndices) {
    completeFeatures = which(complete.cases(pointCloud[mis,]))
    mf = which(!complete.cases(pointCloud[mis,]))
    index = which(ranges[,completeFeatures[1] + (completeFeatures[1] - 1)] <= pointCloud[mis,completeFeatures[1]] & ranges[,completeFeatures[1]+completeFeatures[1]] >= pointCloud[mis,completeFeatures[1]])
    completeFeatures = completeFeatures[-1]
    
    #Find bins whose ranges include the complete features of the observation
    for (cf in completeFeatures){
      if(length(index) == 0){
        index = which(ranges[,cf + (cf - 1)] <= pointCloud[mis,cf] & ranges[,cf+cf] >= pointCloud[mis,cf])
      }
      if(cf %in% categoricalVars)
      {
        for(j in 1:nrow(ranges)){
          if(pointCloud[mis,cf] %in% pointCloud[bins[[j]],cf]){
            ind = append(ind,j)
          }
        }
      }else{
        ind = which(ranges[,cf + (cf - 1)] <= pointCloud[mis,cf] & ranges[,cf+cf] >= pointCloud[mis,cf])
      }
      if(length(intersect(index,ind)) == 0){
       next
      }
      else
      {
        index = intersect(index,ind)
      }
    }
    
    len1 = length(index)
    if(len1 > 0){
 
     vals = matrix(0,nrow = length(index),ncol = length(mf))
    
    #Find BMU from each SOM on the set of mapper bins
    
      training = index
      for(N in 1:length(index)){
        points = bins[[training[N]]]
        len2 = length(points)
        if(len2 == 1){
            vals[N,] = pointCloudComplete[points,mf]
          }else{
          mapping = map(maps[[training[N]]],t(pointCloud[mis,]),maxNA.fraction = 1)
          unit = mapping$unit.classif
          vals[N,] = as.numeric(codes[[training[N]]][unit,mf])
        }
      }
        pcComplete[mis,mf] = colMeans(vals)
    }else{
    #If the complete features fall outside the ranges of all bins, impute with the corresponding complete features of the nearest neighbor
     nn = knn.index(rbind(pointCloudComplete[,-mf],pointCloud[mis,complete.cases(pointCloud[mis,])]),k = 1)[dim(pointCloudComplete)[1] + 1]
     pcComplete[mis,mf] = pointCloudComplete[nn,mf]
    }
  }
  print(proc.time() - ptm) 
  l1 = lm(pcComplete[,7] ~ pcComplete[,6] + pcComplete[,5] + pcComplete[,4] + pcComplete[,3] + pcComplete[,2] + pcComplete[,1])
  mapper[i,] = l1$coefficients
  ci1 = confint(l1,level = 0.95)
  withinCIMAPPER[1] = ci1[1,1] <= beta[1] & ci1[1,2] >= beta[1]
  withinCIMAPPER[2] = ci1[7,1] <= beta[2] & ci1[7,2] >= beta[2]
  withinCIMAPPER[3] = ci1[6,1] <= beta[3] & ci1[6,2] >= beta[3]
  withinCIMAPPER[4] = ci1[5,1] <= beta[4] & ci1[5,2] >= beta[4]
  withinCIMAPPER[5] = ci1[4,1] <= beta[5] & ci1[4,2] >= beta[5]
  withinCIMAPPER[6] = ci1[3,1] <= beta[6] & ci1[3,2] >= beta[6]
  withinCIMAPPER[7] = ci1[2,1] <= beta[7] & ci1[2,2] >= beta[7]
  mapperInt[i,] = withinCIMAPPER
  mapperRMSE[i,1] = rmse(beta[1],mapper[i,1])
  mapperRMSE[i,2] = rmse(beta[2],mapper[i,7])
  mapperRMSE[i,3] = rmse(beta[3],mapper[i,6])
  mapperRMSE[i,4] = rmse(beta[4],mapper[i,5])
  mapperRMSE[i,5] = rmse(beta[5],mapper[i,4])
  mapperRMSE[i,6] = rmse(beta[6],mapper[i,3])
  mapperRMSE[i,7] = rmse(beta[7],mapper[i,2])
  mapperBias[i,1] = (beta[1] - mapper[i,1]) * -1
  mapperBias[i,2] = (beta[2] - mapper[i,7]) * -1
  mapperBias[i,3] = (beta[3] - mapper[i,6]) * -1
  mapperBias[i,4] = (beta[4] - mapper[i,5]) * -1
  mapperBias[i,5] = (beta[5] - mapper[i,4]) * -1
  mapperBias[i,6] = (beta[6] - mapper[i,3]) * -1
  mapperBias[i,7] = (beta[7] - mapper[i,2]) * -1
  
  for (l in 1:10) {
    
  dtComplete = complete(dtMi,l)
  rfComplete = complete(rfMi,l)
  
  
  l2 = lm(dtComplete[,7] ~ dtComplete[,6] + dtComplete[,5] + dtComplete[,4] + dtComplete[,3] + dtComplete[,2] + dtComplete[,1])
  l3 = lm(rfComplete[,7] ~ rfComplete[,6] + rfComplete[,5] + rfComplete[,4] + rfComplete[,3] + rfComplete[,2] + rfComplete[,1])
  
  
  ci2 = confint(l2,level = 0.95)
  ci3 = confint(l3,level = 0.95)
  
  cfDT[l,] = l2$coefficients
  cfRF[l,] = l3$coefficients
  
  
  
  withinCIDT[l,1] = ci2[1,1] <= beta[1] & ci2[1,2] >= beta[1]
  withinCIDT[l,2] = ci2[7,1] <= beta[2] & ci2[7,2] >= beta[2]
  withinCIDT[l,3] = ci2[6,1] <= beta[3] & ci2[6,2] >= beta[3]
  withinCIDT[l,4] = ci2[5,1] <= beta[4] & ci2[5,2] >= beta[4]
  withinCIDT[l,5] = ci2[4,1] <= beta[5] & ci2[4,2] >= beta[5]
  withinCIDT[l,6] = ci2[3,1] <= beta[6] & ci2[3,2] >= beta[6]
  withinCIDT[l,7] = ci2[2,1] <= beta[7] & ci2[2,2] >= beta[7]
  
  withinCIRF[l,1] = ci3[1,1] <= beta[1] & ci3[1,2] >= beta[1]
  withinCIRF[l,2] = ci3[7,1] <= beta[2] & ci3[7,2] >= beta[2]
  withinCIRF[l,3] = ci3[6,1] <= beta[3] & ci3[6,2] >= beta[3]
  withinCIRF[l,4] = ci3[5,1] <= beta[4] & ci3[5,2] >= beta[4]
  withinCIRF[l,5] = ci3[4,1] <= beta[5] & ci3[4,2] >= beta[5]
  withinCIRF[l,6] = ci3[3,1] <= beta[6] & ci3[3,2] >= beta[6]
  withinCIRF[l,7] = ci3[2,1] <= beta[7] & ci3[2,2] >= beta[7]
  
}

 
dt[i,] = colSums(cfDT) / 10
dtInt[i,] = colSums(withinCIDT) / 10

dtRMSE[i,1] = rmse(beta[1],dt[i,1])
dtRMSE[i,2] = rmse(beta[2],dt[i,7])
dtRMSE[i,3] = rmse(beta[3],dt[i,6])
dtRMSE[i,4] = rmse(beta[4],dt[i,5])
dtRMSE[i,5] = rmse(beta[5],dt[i,4])
dtRMSE[i,6] = rmse(beta[6],dt[i,3])
dtRMSE[i,7] = rmse(beta[7],dt[i,2])

dtBias[i,1] = (beta[1] - dt[i,1]) * -1
dtBias[i,2] = (beta[2] - dt[i,7]) * -1
dtBias[i,3] = (beta[3] - dt[i,6]) * -1
dtBias[i,4] = (beta[4] - dt[i,5]) * -1
dtBias[i,5] = (beta[5] - dt[i,4]) * -1
dtBias[i,6] = (beta[6] - dt[i,3]) * -1
dtBias[i,7] = (beta[7] - dt[i,2]) * -1

rf[i,] = colSums(cfRF) / 10
rfInt[i,] = colSums(withinCIRF) / 10

rfRMSE[i,1] = rmse(beta[1],rf[i,1])
rfRMSE[i,2] = rmse(beta[2],rf[i,7])
rfRMSE[i,3] = rmse(beta[3],rf[i,6])
rfRMSE[i,4] = rmse(beta[4],rf[i,5])
rfRMSE[i,5] = rmse(beta[5],rf[i,4])
rfRMSE[i,6] = rmse(beta[6],rf[i,3])
rfRMSE[i,7] = rmse(beta[7],rf[i,2])

rfBias[i,1] = (beta[1] - rf[i,1]) * -1
rfBias[i,2] = (beta[2] - rf[i,7]) * -1
rfBias[i,3] = (beta[3] - rf[i,6]) * -1
rfBias[i,4] = (beta[4] - rf[i,5]) * -1
rfBias[i,5] = (beta[5] - rf[i,4]) * -1
rfBias[i,6] = (beta[6] - rf[i,3]) * -1
rfBias[i,7] = (beta[7] - rf[i,2]) * -1

pred1 = newData[,1:6] %*% rev(mapper[i,2:7]) + mapper[i,1]
pred2 = newData[,1:6] %*% rev(dt[i,2:7]) + dt[i,1]
pred3 = newData[,1:6] %*% rev(rf[i,2:7]) + rf[i,1]

ARMSE[i,1] = rmse(newData[,7],pred1)
ARMSE[i,2] = rmse(newData[,7],pred2)
ARMSE[i,3] = rmse(newData[,7],pred3)
}
