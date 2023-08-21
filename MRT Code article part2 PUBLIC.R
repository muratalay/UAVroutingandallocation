#Set Working Directory 
setwd("workingdirectory")
nrep<-1
sink("vPart2_Demonstration_DecSpace_1.txt", append=TRUE)
cat("nrep, nTarget, nAppLoc, MissionTimeLimit:, TotalTimeTour:, wlimit, UAV_Type, Status, OptimalityGap, RunTimeGivenbyGurobi, RunTimeGivenbyR, ObjectiveValue, teta_result, BestAvailableUpperBound, TotalTimeTour, TargetID_Type_LocX_LocY, x_ijkl_variables, yia_variables, sia_results, sia_variables")
cat("\n")
sink()

sink("vPart2_Demonstration_Summary_1.txt", append=TRUE)
cat("nrep, nTarget, nAppLoc, MissionTimeLimit, TotalTimeTour, wlimit, UAV_Type, TotalRadarThreat, result$status, result$objval, result$objbound, result$mipgap, result$runtime, runTime, TargetID_Type_LocX_LocY")
cat("\n")
sink()

#Load Packages
library(Matrix)
library(gurobi)
library(rdist, lib.loc="/user/erdidasd/R_Packages")

#Set Global Parameters
options(max.print = 99999)
SeedAll<-c("define seeds")
#---------------------------------------------------------------------------------

  nTarget<-10
  nAppLoc<-4
nTargetType<-4 
wTargetType<-c(1,1.4,1.6,2) 
nPayloadType<-3 
MissionArea<-1000 
FlightV<- 250 
MissionTimeLimitSet<-c(10,15,20)

#UAV Configuration Settings
UAV_Conf<-c(1,2,3)
UAV_podQuant<-c(10,6,4) 
UAV_Payload<-c(1500,750,350) 

PayloadRangeMaxMinLoadQuant<-matrix(c(14,25,25,3,5,5,125,250,325,6,3,2,6,2,2,4,2,0),nrow=nPayloadType,ncol=6) 
R<-max(PayloadRangeMaxMinLoadQuant[,1]) 

MinProbMat<-matrix(c(0.7,0.8,0,0.8,0.8,0.9,0,0.85,0.9,0,0,0.8),nrow=nPayloadType,ncol=nTargetType)
MaxProbMat<-matrix(c(0.8,0.9,0,0.9,0.9,0.99,0,0.9,0.99,0,0,0.9),ncol=nTargetType,nrow=nPayloadType)

#-------------
FOV_equation<-c(0.3176,-0.01372) 
lambda<-c(0.01,0.03) 

wSet<-c(0.1,0.15,0.2) 
#-------------
expSetResults<-expand.grid(MTL=MissionTimeLimitSet,RT=wSet,UAV_Type=UAV_Conf)
#-------------

  set.seed(SeedAll[nrep])


    nNodes<- nTarget + 1 
    
    i=2
    TargetLoc<-matrix(data=NA,nrow=nTarget+1,ncol=2)
    TargetLoc[1,]=0
    
    DistMtrx<-c(rep(1,nNodes*nNodes))
    VDM=matrix(DistMtrx,nrow=length(DistMtrx),ncol=1)
    while ((sum(VDM[,1]>0 & VDM[,1]<(2*R+5)))>0  ){
    for (i in 2:(nTarget+1))
    {
      j=1
      for (j in 1:2) {
        
        TargetLoc[i,j]=round(runif(1,min=(2*R),max=(MissionArea-2*R))) 
      }}
      BaseLoc<-matrix(0,nrow=1,ncol=2)
      DistMtrx<-pdist(TargetLoc)
      VDM=matrix(DistMtrx,nrow=length(DistMtrx),ncol=1)
    }
    
    TargetLoc<-TargetLoc[order(TargetLoc[,1],decreasing=FALSE),]

    TargetID<-c(0:nTarget)
    TargetType<-c(0,round(runif(nTarget,min=1,max=nTargetType)))
    TargetID_Type_Loc<-cbind(TargetID,TargetType,TargetLoc) 

    wClass<-spMatrix(nTarget,1)
    for(i in 2:length(TargetType)){
    wClass[i-1,1]<-c(wTargetType[TargetType[i]])
    }
    w<-as.vector(wClass[,1]) 
    #-------------
    
    #check locations
    plot(TargetLoc[,1],TargetLoc[,2])

          angle=(360/nAppLoc)
    
    TargetAppLoc<-array(NA,dim = c((nTarget+1),(nAppLoc+1),2))
    i=1
    j=1
    for (i in 1:(nTarget+1)){
      TargetAppLoc[i,j,1]<- TargetLoc[i,1]
      TargetAppLoc[i,j,2]<- TargetLoc[i,2]
    }
      i=1
      j=2
    
    for (i in 2:(nTarget+1))
    {
      
      for (j in 2:(nAppLoc+1))
      {
          
      TargetAppLoc[i,j,1]<- TargetLoc[i,1]+sin(angle*j*pi/180)*R 
        TargetAppLoc[i,j,2]<- TargetLoc[i,2]+cos(angle*j*pi/180)*R
      }
    }
      
      plot(TargetAppLoc[,,1],TargetAppLoc[,,2])
        
    Xdata<-matrix(TargetAppLoc[,,1],nrow=length(TargetAppLoc[,,1]),1)
    Ydata<-matrix(TargetAppLoc[,,2],nrow=length(TargetAppLoc[,,2]),1)
    XYdata<-cbind(Xdata,Ydata)
    
    FromToChart<-pdist(XYdata) 
    FromToChartARRY<-array(NA,dim=c((nTarget+1),(nAppLoc+1),(nTarget+1),(nAppLoc+1)))
    FromToChart[is.nan(FromToChart)]=NA
    i=1
    j=1
    k=1
    l=1
    
    for (i in 1:(nTarget+1))
    {
      for (j in 1:(nAppLoc+1))
      {
        for (k in 1:(nTarget+1))
        {
          for (l in 1:(nAppLoc+1))
          {
            
            FromToChartARRY[i,j,k,l]<- FromToChart[(i+(j-1)*(nTarget+1)),(k+(l-1)*(nTarget+1))]
    
          }
        }
      }
    }
    
    TangDist <- function(Px,Py,Cx,Cy,a){ 
      
        if(!is.nan(sqrt((Px - Cx)^2 + (Py - Cy)^2))){
        b=sqrt((Px - Cx)^2 + (Py - Cy)^2)
        
        if (!is.nan(acos(a / b))) {
        th = acos(a / b) 

        d = atan2(Py - Cy, Px - Cx)  
        d1 = d + th  
        d2 = d - th  
        
        T1x = Cx + a * cos(d1)
        T1y = Cy + a * sin(d1)
        T2x = Cx + a * cos(d2)
        T2y = Cy + a * sin(d2)
        
        if (!is.nan(sqrt((Px - T1x)^2 + (Py - T1y)^2))) {
          k1=sqrt((Px - T1x)^2 + (Py - T1y)^2)
        
          }else k1=NA
        }else k1=NA
      }else k1=NA
    }

    TangDistArray<-array(NA,dim=c((nTarget+1),(nAppLoc+1),(nTarget+1))) 
    i=1
    j=1
    k=1
    for (i in 1:(nTarget+1))
    {
      i
      for (j in 1:(nAppLoc+1))
      {
        j
        for (k in 1:(nTarget+1))
        {
          if (i==k) {
            TangDistArray[i,j,k]=NA 
    
          }else if(k==1){
          
          TangDistArray[i,j,k]<- TangDist(XYdata[(i+(j-1)*(nTarget+1)),1],XYdata[(i+(j-1)*(nTarget+1)),2],XYdata[k,1],XYdata[k,2],0)
          
          }else
          TangDistArray[i,j,k]<- TangDist(XYdata[(i+(j-1)*(nTarget+1)),1],XYdata[(i+(j-1)*(nTarget+1)),2],XYdata[k,1],XYdata[k,2],R)
          
        }
      }
    }
    
    variabledef<-array(NA,dim=c((nTarget+1),(nAppLoc+1),(nTarget+1),(nAppLoc+1))) 
    variabledef<-FromToChartARRY     
    
    X_ijkl_matrix<-c()
    for (i in 1:(nTarget+1)) 
    {
      for (j in 1:(nAppLoc+1))
      {
        for (k in 1:(nTarget+1))
        {
          for (l in 1:(nAppLoc+1))
          {
          if(i==k || (!(i==1) && j==1)){
          variabledef[i,j,k,l]=0
          }else if(!is.na(FromToChartARRY[i,j,k,l]) && !is.na(TangDistArray[i,j,k]) && FromToChartARRY[i,j,k,l]<=TangDistArray[i,j,k]) {
          variabledef[i,j,k,l]=1
          X_ijkl_matrix<-rbind(X_ijkl_matrix,as.matrix(cbind(i,j,k,l,FromToChartARRY[i,j,k,l])))
          }else
          variabledef[i,j,k,l]=0
            
          }}}}
          
          X_ijkl_matrix2<-matrix(cbind(X_ijkl_matrix[,3:4],X_ijkl_matrix[,1:2],X_ijkl_matrix[,5]),ncol=5)
          X_ijkl_matrix3<-rbind(X_ijkl_matrix,X_ijkl_matrix2)
          X_ijkl_matrix4<-unique(X_ijkl_matrix3)
          X_ijkl_matrix<-X_ijkl_matrix4
          X_ijkl_matrix<-X_ijkl_matrix[order(X_ijkl_matrix[,1],decreasing=FALSE),]
          X_ijkl_matrix[,5]<-round(X_ijkl_matrix[,5],4)
          X_ijkl_matrix<-unique(X_ijkl_matrix)
          
    PayloadTargetCompMtrx<-matrix(NA,nrow=nPayloadType*nTarget,ncol=(3+2+2)) 
    l=1
    i=1
      for(i in 1:(nTarget)){
          for(k in 1:nPayloadType){
      PayloadTargetCompMtrx[l,]<-c(i,TargetType[i+1],k,MinProbMat[k,TargetType[i+1]],MaxProbMat[k,TargetType[i+1]],PayloadRangeMaxMinLoadQuant[k,2],PayloadRangeMaxMinLoadQuant[k,1])
      l=l+1
        }
      }
    colnames(PayloadTargetCompMtrx)<-c("TargetID","TargetType","PayloadType","MinProb","MaxProb","MinDist","MaxDist")
    
    PayloadTargetCompMtrx2<-cbind(PayloadTargetCompMtrx, 
                               (PayloadTargetCompMtrx[,4]-PayloadTargetCompMtrx[,5])/ (PayloadTargetCompMtrx[,7]-PayloadTargetCompMtrx[,6]),
                               PayloadTargetCompMtrx[,5]-PayloadTargetCompMtrx[,6]*((PayloadTargetCompMtrx[,4]-PayloadTargetCompMtrx[,5])/ (PayloadTargetCompMtrx[,7]-PayloadTargetCompMtrx[,6])))
    colnames(PayloadTargetCompMtrx2)<-c("TargetID","TargetType","PayloadType","MinProb","MaxProb","MinDist","MaxDist","C3 ratio","C3 rhs")
    
                          
    #---------------------------------------------------------------------------------------------------------------
    #DOE FACTOR DESIGN
    
    for (i in 1:nrow(expSetResults)){
      
      MissionTimeLimit<-expSetResults[i,1]
      wlimit<-expSetResults[i,2]
      UAV_Type<-expSetResults[i,3]

    #DOE FACTOR DESIGN end
      
    #---------------------------------------------------------------------------------------------------------------
                          #START MIP MODEL
                          
    
    #VARIABLE CONSTRUCTION
    
              #Xijkl Binary Variables: Variables for arcs from ith target's jth waypoint to kth target's lth waypoint
              X_ijkl_matrix<-X_ijkl_matrix #already defined previously. Matrix defining the available arcs between waypoints and its euclidean distances.
              colnames(X_ijkl_matrix)<-c("i","j","k","l","EucDist")
              X_ijkl_ord<-as.matrix(cbind(paste0(X_ijkl_matrix[,1],",",X_ijkl_matrix[,2],",",X_ijkl_matrix[,3],",",X_ijkl_matrix[,4] )))
              X_ijkl_names<-apply(X_ijkl_matrix[,1:4],1, function (x) paste0('x',",",x[1],",",x[2],",",x[3],",",x[4]))
              X_ijkl_obj<-rep(0,length(X_ijkl_names)) #objective function coefficients of X binary variables.
              
              #Yia If UCAV send 'i' type of target is hit by 'a' type of Payload
              Y_ia_matrix<-as.matrix(cbind(cbind(rep(TargetID[2:length(TargetID)], each=nPayloadType),1:nPayloadType),as.vector(PayloadRangeMaxMinLoadQuant[,3])))
              colnames(Y_ia_matrix)<-c("TargetID","PayloadType","PayloadLoad")
              Y_ia_ord<-as.matrix(cbind(paste0(Y_ia_matrix[,1],",",Y_ia_matrix[,2])))
              Y_ia_names<-apply(Y_ia_matrix[,1:2],1, function (x) paste0('y',",",x[1],",",x[2]))
              Y_ia_obj<-rep(0,length(Y_ia_names)) #objective function coefficients of Y binary variables.
              
              #Sia Distance where 'i' type of target is hit by 'a' type of Payload
              S_ia_matrix<-as.matrix(cbind(cbind(rep(TargetID[2:length(TargetID)], each=nPayloadType),1:nPayloadType), as.vector(PayloadRangeMaxMinLoadQuant[,2]),as.vector(PayloadRangeMaxMinLoadQuant[,1])))
              colnames(S_ia_matrix)<-c("TargetID","PayloadType","MinDist","MaxDist")
              S_ia_ord<-as.matrix(cbind(paste0(S_ia_matrix[,1],",",S_ia_matrix[,2])))
              S_ia_names<-apply(S_ia_matrix[,1:2],1, function (x) paste0('s',",",x[1],",",x[2]))
              S_ia_obj<-rep(0,length(S_ia_names)) #objective function coefficients of S variables.
              
              #Pia Probability where 'i' type of target is destroyed by 'a' type of Payload
              P_ia_matrix<-as.matrix(cbind(rep(TargetID[2:length(TargetID)], each=nPayloadType),1:nPayloadType,rep(w,each=nPayloadType)))
              colnames(P_ia_matrix)<-c("TargetID","PayloadType","TargetWeights")
              P_ia_ord<-as.matrix(cbind(paste0(P_ia_matrix[,1],",",P_ia_matrix[,2])))
              P_ia_names<-apply(P_ia_matrix[,1:2],1, function (x) paste0('p',",",x[1],",",x[2]))
              P_ia_obj<-P_ia_matrix[,3] #objective function coefficients of P variables.
              
              #Dummy Variable u defined (Only includes target approach locations-no origin)
              u_ij_matrix<-X_ijkl_matrix[,1:2]
              u_ij_matrix<-unique(u_ij_matrix)
              u_ij_matrix<-u_ij_matrix[-1,] #origin haric
              colnames(u_ij_matrix)<-c("i","j")
              u_ij_ord<-as.matrix(cbind(paste0(u_ij_matrix[,1],",",u_ij_matrix[,2])))
              u_ij_names<-apply(u_ij_matrix[,1:2],1, function (x) paste0('u',",",x[1],",",x[2]))
              u_ij_obj<-rep(0,length(u_ij_names)) #objective function coefficients of u dummy variables.
              
              teta_ij_matrix<-c(1)
              teta_ij_names<-c("teta")
              teta_ij_obj<-0 
              
              #All variables - ordered
              all_names<-c(X_ijkl_names,Y_ia_names,S_ia_names,P_ia_names, u_ij_names,teta_ij_names)
              all_obj<-c(X_ijkl_obj,Y_ia_obj,S_ia_obj,P_ia_obj,u_ij_obj,teta_ij_obj)
    
    
    #CONSTRAINT CONSTRUCTION
            #Constraint 3 - Destruction Probability Calculator for Target "i" (of a specific type) hit by Payload type "a"
              
              C3_Siar<-c()
              C3_Piar<-c()
              C3_nRow<-nTarget*nPayloadType
              satir<-1
              for(satir in 1:(C3_nRow)){
                C3_Sia<-spMatrix(1,(length(S_ia_names)), i=1, j=satir, x= (-PayloadTargetCompMtrx2[satir,8]))
                C3_Pia<-spMatrix(1,(length(P_ia_names)), i=1, j=satir, x= 1)
                C3_Siar<-rbind(C3_Siar,C3_Sia)
                C3_Piar<-rbind(C3_Piar,C3_Pia)
              }
              
              C3_Xijklr<-spMatrix(C3_nRow,(length(X_ijkl_names)))
              C3_Yiar<-spMatrix(C3_nRow,(length(Y_ia_names)))
              C3_Uijr<-spMatrix(C3_nRow,(length(u_ij_names)))
              C3_teta<-spMatrix(C3_nRow,(length(teta_ij_names)))
              
              C3<- cbind(C3_Xijklr,C3_Yiar,C3_Siar,C3_Piar,C3_Uijr,C3_teta)
              C3_rhs<-PayloadTargetCompMtrx2[,9]
              C3_dir<-c(rep("<=",nrow(C3)))
           
            #Constraint 4 - Hit Distance Intervals
              
              C4_nRow<-(nPayloadType*nTarget) #number of rows (constraints)
    
              C4.1_Xijkl<-spMatrix(C4_nRow,(length(X_ijkl_names)))
              C4.1_Yia<-diag(-S_ia_matrix[,4],nrow=(length(Y_ia_names)),ncol =(length(Y_ia_names)))
              C4.1_Sia<-diag(1,nrow=(length(S_ia_names)),ncol =(length(S_ia_names))) 
              C4.1_Pia<-spMatrix(C4_nRow,(length(P_ia_names)))
              C4.1_Uij<-spMatrix(C4_nRow,(length(u_ij_names)))    
              C4.1_teta<-spMatrix(C4_nRow,(length(teta_ij_names)))
              
              C4.1<- cbind(C4.1_Xijkl,C4.1_Yia,C4.1_Sia,C4.1_Pia,C4.1_Uij,C4.1_teta)
              C4.1_rhs<-c(rep(0,nrow(C4.1)))
              C4.1_dir<-c(rep("<=",nrow(C4.1)))
          
              C4.2_Xijkl<-spMatrix(C4_nRow,(length(X_ijkl_names))) 
              C4.2_Yia<-diag(-S_ia_matrix[,3],nrow=(length(Y_ia_names)),ncol =(length(Y_ia_names)))
              C4.2_Sia<-diag(1,nrow=(length(S_ia_names)),ncol =(length(S_ia_names))) 
              C4.2_Pia<-spMatrix(C4_nRow,(length(P_ia_names)))
              C4.2_Uij<-spMatrix(C4_nRow,(length(u_ij_names)))
              C4.2_teta<-spMatrix(C4_nRow,(length(teta_ij_names)))
              
              C4.2<- cbind(C4.2_Xijkl,C4.2_Yia,C4.2_Sia,C4.2_Pia,C4.2_Uij,C4.2_teta)
              C4.2_rhs<-c(rep(0,nrow(C4.2)))
              C4.2_dir<-c(rep(">=",nrow(C4.2)))
              
              C4<-rbind(C4.1,C4.2)
              C4_rhs<-c(C4.1_rhs,C4.2_rhs)
              C4_dir<-c(C4.1_dir,C4.2_dir)
              
            #Constraint 5 - Gain only if the target is hit
              
              C5_Xijkl<-Matrix(0,nrow=(nPayloadType*nTarget),ncol=(length(X_ijkl_names)))
              C5_Yia<-diag(-1,nrow=(nPayloadType*nTarget),ncol=length(Y_ia_names))
              C5_Sia<-Matrix(0,nrow=(length(S_ia_names)),ncol =(length(S_ia_names)))
              C5_Pia<-diag(1,nrow=(nPayloadType*nTarget),ncol=length(P_ia_names))
              C5_Uij<-Matrix(0,nrow=(nPayloadType*nTarget), ncol=(length(u_ij_names)))
              C5_teta<-Matrix(0,nrow=(nPayloadType*nTarget), ncol=(length(teta_ij_names)))
              
              C5<- cbind(C5_Xijkl,C5_Yia,C5_Sia,C5_Pia,C5_Uij,C5_teta)
              C5_rhs<-c(rep(0,length(P_ia_names)))
              C5_dir<-c(rep("<=",nrow(C5)))
              
            #Constraint 6 - First Visit a Target, Then Hit it. Visit and Hit only once at most.
              
              C6_Xijklr<-c()
              C6_Yiar<-c()
    
              for(nti in 1:nTarget){
                C6_Xijkl<-spMatrix(1,(length(X_ijkl_names)), i=rep(1, (length(which(X_ijkl_matrix[,3]==(nti+1))))), j=which(X_ijkl_matrix[,3]==(nti+1)), x= c(rep(-1,length(which(X_ijkl_matrix[,3]==(nti+1))))))
                C6_Yia<-spMatrix(1,(length(Y_ia_names)), i=rep(1, (length(which(Y_ia_matrix[,1]==(nti))))), j=which(Y_ia_matrix[,1]==(nti)), x= c(rep(1, length(which(Y_ia_matrix[,1]==(nti))))))
                
                C6_Xijklr<-rbind(C6_Xijklr,C6_Xijkl)
                C6_Yiar<-rbind(C6_Yiar,C6_Yia)
              }
              
              C6_Siar<-spMatrix(nTarget,(length(S_ia_names)))
              C6_Piar<-spMatrix(nTarget,(length(P_ia_names)))
              C6_Uijr<-spMatrix(nTarget,(length(u_ij_names)))
              C6_teta<-spMatrix(nTarget,(length(teta_ij_names)))
              
              C6<- cbind(C6_Xijklr,C6_Yiar,C6_Siar,C6_Piar,C6_Uijr,C6_teta)
              C6_rhs<-c(rep(0,nTarget))
              C6_dir<-c(rep("=",nrow(C6)))
              
            #Constraint 7 - At most 1 Payload per target
              
              C7_Yiar<-c()
    
              for(nti in 1:nTarget){
                
                C7_Yia<-spMatrix(1,(length(Y_ia_names)), i=rep(1, (length(which(Y_ia_matrix[,1]==(nti))))), j=which(Y_ia_matrix[,1]==(nti)), x= c(rep(1, length(which(Y_ia_matrix[,1]==(nti))))))
                C7_Yiar<-rbind(C7_Yiar,C7_Yia)
    
              }
              
              C7_Xijklr<-spMatrix(nTarget,(length(X_ijkl_names)))
              C7_Siar<-spMatrix(nTarget,(length(S_ia_names)))
              C7_Piar<-spMatrix(nTarget,(length(P_ia_names)))
              C7_Uijr<-spMatrix(nTarget,(length(u_ij_names)))
              C7_teta<-spMatrix(nTarget,(length(teta_ij_names)))
              
              C7<- cbind(C7_Xijklr,C7_Yiar,C7_Siar,C7_Piar,C7_Uijr,C7_teta)
              C7_rhs<-c(rep(1,nTarget))
              C7_dir<-c(rep("<=",nrow(C7)))
              
            #Constraint 8 - Total Payload loaded should be less than total number of pods
              C8_Xijkl<-spMatrix(1,(length(X_ijkl_names)))
              C8_Yia<-spMatrix(1,(length(Y_ia_names)), i=rep(1, (length(Y_ia_names))), j=(1:(length(Y_ia_names))), x= c(rep(1, (length(Y_ia_names)))))
              C8_Sia<-spMatrix(1,(length(S_ia_names)))
              C8_Pia<-spMatrix(1,(length(P_ia_names)))
              C8_Uij<-spMatrix(1,(length(u_ij_names)))
              C8_teta<-spMatrix(1,(length(teta_ij_names)))
              
              C8<- cbind(C8_Xijkl,C8_Yia,C8_Sia,C8_Pia,C8_Uij,C8_teta)
              C8_rhs<-c(UAV_podQuant[UAV_Type])
              C8_dir<-c("<=")
              
            #Constraint 9 - Total payload limit
              C9_Xijkl<-spMatrix(1,(length(X_ijkl_names)))
              C9_Yia<-spMatrix(1,(length(Y_ia_names)), i=rep(1, (length(Y_ia_names))), j=(1:(length(Y_ia_names))), x= c(Y_ia_matrix[,3]))
              C9_Sia<-spMatrix(1,(length(S_ia_names)))
              C9_Pia<-spMatrix(1,(length(P_ia_names)))
              C9_Uij<-spMatrix(1,(length(u_ij_names)))
              C9_teta<-spMatrix(1,(length(teta_ij_names)))
              
              C9<- cbind(C9_Xijkl,C9_Yia,C9_Sia,C9_Pia,C9_Uij,C9_teta)
              C9_rhs<-c(UAV_Payload[UAV_Type])
              C9_dir<-c("<=")
              
              
            #Constraint 10 - At most Q quantities can be used for each target
              
              C10_Yiar<-c()
    
              for(nai in 1:nPayloadType){
                
                C10_Yia<-spMatrix(1,(length(Y_ia_names)), i=rep(1, (length(which(Y_ia_matrix[,2]==(nai))))), j=which(Y_ia_matrix[,2]==(nai)), x= c(rep(1, length(which(Y_ia_matrix[,2]==(nai))))))
                C10_Yiar<-rbind(C10_Yiar,C10_Yia)
    
              }
              
              C10_Xijklr<-spMatrix(nPayloadType,(length(X_ijkl_names)))
              C10_Siar<-spMatrix(nPayloadType,(length(S_ia_names)))
              C10_Piar<-spMatrix(nPayloadType,(length(P_ia_names)))
              C10_Uijr<-spMatrix(nPayloadType,(length(u_ij_names)))
              C10_teta<-spMatrix(nPayloadType,(length(teta_ij_names)))
              
              C10<- cbind(C10_Xijklr,C10_Yiar,C10_Siar,C10_Piar,C10_Uijr,C10_teta)
              C10_rhs<-c(PayloadRangeMaxMinLoadQuant[,(4+UAV_Type-1)])
              C10_dir<-c(rep("<=",nrow(C10)))
              
              
            #Balance Equations
              #Constraint 11 - Incoming Active arcs to a target must be equal to outgoing active arcs
              
              wPoints<-unique(X_ijkl_matrix[,1:2]) #All waypoints
              ntw<-1
              
              C11_Xijkl<-c()
              for(ntw in 1:nrow(wPoints)){
                
                wPoint<-wPoints[ntw,]
                inTowPoint<-which(X_ijkl_matrix[,3]==wPoint[1] &X_ijkl_matrix[,4]==wPoint[2] ) #indices of the incoming arcs to wPoint
                outFromwPoint<-which(X_ijkl_matrix[,1]==wPoint[1] &X_ijkl_matrix[,2]==wPoint[2] ) #indices of the outgoing arcs from wPoint
                
                C11_Xijkl_1<-spMatrix(1,(length(X_ijkl_names)), i=rep(1, length(inTowPoint)),j=inTowPoint, x= rep(1, length(inTowPoint)))
                C11_Xijkl_2<-spMatrix(1,(length(X_ijkl_names)), i=rep(1, length(outFromwPoint)), j=outFromwPoint, x= c(rep(-1, length(outFromwPoint))))
                
                C11_Xijkl<-rbind(C11_Xijkl,(C11_Xijkl_1+C11_Xijkl_2))
                
              }
              
              C11_Yia<-spMatrix(nrow(wPoints),(length(Y_ia_names)))
              C11_Sia<-spMatrix(nrow(wPoints),(length(S_ia_names)))
              C11_Pia<-spMatrix(nrow(wPoints),(length(P_ia_names)))
              C11_Uij<-spMatrix(nrow(wPoints),(length(u_ij_names)))
              C11_teta<-spMatrix(nrow(wPoints),(length(teta_ij_names)))
              
              
              C11<- cbind(C11_Xijkl,C11_Yia,C11_Sia,C11_Pia,C11_Uij,C11_teta)
              C11_rhs<-c(rep(0,nrow(C11)))
              C11_dir<-c(rep("=",nrow(C11)))
              
              #Constraint 12 - Subtour Elimination_1 - Dummy
              
              wPoints<-unique(X_ijkl_matrix[,1:2]) #All waypoints
              wPoints<-wPoints[-1,] #All waypoints except base(origin)
              ntw<-1
              
              C12_Xijkl<-c()
              C12_Uij<-c()
              counter<-0
              for(ntw in 1:(nrow(wPoints))){
                
                wPointk<-wPoints[ntw,]
                wPointl<-matrix(wPoints[wPoints[,1]!=wPoints[ntw,1]],ncol=2)
                ntkl<-1
                
                for(ntkl in 1:length(wPointl[,1])){
                
                  if(length(which(X_ijkl_matrix[,1]==wPointk[1] &X_ijkl_matrix[,2]==wPointk[2] & X_ijkl_matrix[,3]==wPointl[ntkl,1] &X_ijkl_matrix[,4]==wPointl[ntkl,2]))!=0){
                
                Xkl<-which(X_ijkl_matrix[,1]==wPointk[1] &X_ijkl_matrix[,2]==wPointk[2] & X_ijkl_matrix[,3]==wPointl[ntkl,1] &X_ijkl_matrix[,4]==wPointl[ntkl,2] ) #indices of the outgoing arcs from wPoint
                Xlk<-which(X_ijkl_matrix[,3]==wPointk[1] &X_ijkl_matrix[,4]==wPointk[2] & X_ijkl_matrix[,1]==wPointl[ntkl,1] &X_ijkl_matrix[,2]==wPointl[ntkl,2] ) #indices of the incoming arcs to wPoint
                
                uk<-which(u_ij_matrix[,1]==wPointk[1] & u_ij_matrix[,2]==wPointk[2]) #indices of the outgoing arcs from wPoint
                ul<-which(u_ij_matrix[,1]==wPointl[ntkl,1] & u_ij_matrix[,2]==wPointl[ntkl,2]) #indices of the incoming arcs to wPoint
                
                C12_Xijkl_1<-spMatrix(1,(length(X_ijkl_names)), i=1,j=Xkl, x= nNodes-1)
                C12_Xijkl_2<-spMatrix(1,(length(X_ijkl_names)), i=1,j=Xlk, x= nNodes-3)
                
                C12_Uij_1<-spMatrix(1,(length(u_ij_names)), i=1,j=uk, x= 1)
                C12_Uij_2<-spMatrix(1,(length(u_ij_names)), i=1,j=ul, x= -1)
                
                C12_Xijkl<-rbind(C12_Xijkl,(C12_Xijkl_1+C12_Xijkl_2))
                C12_Uij<-rbind(C12_Uij,(C12_Uij_1+C12_Uij_2))
                
                counter<-counter+1}
                }
              }
              counter
              C12_Yia<-spMatrix(counter,(length(Y_ia_names)))
              C12_Sia<-spMatrix(counter,(length(S_ia_names)))
              C12_Pia<-spMatrix(counter,(length(P_ia_names)))
              C12_teta<-spMatrix(counter,(length(teta_ij_names)))
              
              C12<- cbind(C12_Xijkl,C12_Yia,C12_Sia,C12_Pia,C12_Uij,C12_teta)
              C12_rhs<-c(rep((nNodes-2),counter))
              C12_dir<-c(rep("<=",counter))
              
            #Constraint 13 - Source Limit (Baseden Cikislar en fazla 1)
              
              C13_Xijkl<-spMatrix(1,(length(X_ijkl_names)), i=rep(1, length(which(X_ijkl_matrix[,1]==1))), j=(1:(length(which(X_ijkl_matrix[,1]==1)))), x=rep(1, length(which(X_ijkl_matrix[,1]==1))))
              C13_Yia<-spMatrix(1,(length(Y_ia_names)))
              C13_Sia<-spMatrix(1,(length(S_ia_names))) 
              C13_Pia<-spMatrix(1,(length(P_ia_names)))
              C13_Uij<-spMatrix(1,(length(u_ij_names))) 
              C13_teta<-spMatrix(1,(length(teta_ij_names)))  
              C13<-cbind(C13_Xijkl,C13_Yia,C13_Sia,C13_Pia,C13_Uij,C13_teta)
              C13_rhs<-c(1)
              C13_dir<-c("<=")
              
              

            #Constraint 2 - Mission Time Limit
              C2_Xijkl<-spMatrix(1,(length(X_ijkl_names)), i=rep(1, (length(X_ijkl_names))), j=(1:(length(X_ijkl_names))), x= c(X_ijkl_matrix[,5]/FlightV))
              C2_Yia<-spMatrix(1,(length(Y_ia_names)), i=rep(1, (length(Y_ia_names))), j=(1:(length(Y_ia_names))), x= c(rep((2*max(PayloadRangeMaxMinLoadQuant[,1])/FlightV), (length(Y_ia_names)))))
              C2_Sia<-spMatrix(1,(length(S_ia_names)), i=rep(1, (length(S_ia_names))), j=(1:(length(S_ia_names))), x= c(rep((-2/FlightV), (length(S_ia_names)))))
              C2_Pia<-spMatrix(1,(length(P_ia_names)), i=rep(1, (length(P_ia_names))), j=(1:(length(P_ia_names))), x= c(rep(0, (length(P_ia_names)))))
              C2_Uij<-spMatrix(1,(length(u_ij_names)), i=rep(1, (length(u_ij_names))), j=(1:(length(u_ij_names))), x= c(rep(0, (length(u_ij_names)))))
              C2_teta<-spMatrix(1,(length(teta_ij_names)), i=rep(1, (length(teta_ij_names))), j=(1:(length(teta_ij_names))), x= c(rep(0, (length(teta_ij_names)))))
              C2<- cbind(C2_Xijkl,C2_Yia,C2_Sia,C2_Pia,C2_Uij,C2_teta)
              C2_rhs<-c(MissionTimeLimit)
              C2_dir<-c("<=")
              


                #Constraint 14 - Radar Detection Threat
                C14_Xijkl<-spMatrix(1,(length(X_ijkl_names)), i=rep(1, (length(X_ijkl_names))), j=(1:(length(X_ijkl_names))), x= c(lambda[1]*X_ijkl_matrix[,5]/FlightV))
                C14_Yia<-spMatrix(1,(length(Y_ia_names)), i=rep(1, (length(Y_ia_names))), j=(1:(length(Y_ia_names))), x= c(rep(FOV_equation[1]*lambda[2], (length(Y_ia_names)))))
                C14_Sia<-spMatrix(1,(length(S_ia_names)), i=rep(1, (length(S_ia_names))), j=(1:(length(S_ia_names))), x= c(rep(FOV_equation[2]*lambda[2], (length(S_ia_names)))))
                C14_Pia<-spMatrix(1,(length(P_ia_names)), i=rep(1, (length(P_ia_names))), j=(1:(length(P_ia_names))), x= c(rep(0, (length(P_ia_names)))))
                C14_Uij<-spMatrix(1,(length(u_ij_names)), i=rep(1, (length(u_ij_names))), j=(1:(length(u_ij_names))), x= c(rep(0, (length(u_ij_names)))))
                C14_teta<-spMatrix(1,(length(teta_ij_names)), i=rep(1, (length(teta_ij_names))), j=(1:(length(teta_ij_names))), x= c(rep(-1, (length(teta_ij_names)))))
                C14<-cbind(C14_Xijkl,C14_Yia,C14_Sia,C14_Pia,C14_Uij,C14_teta)
                C14_rhs<-0
                C14_dir<-c("=")
                
                #Constraint 15 - Radar Detection Threat Limit
                C15_Xijkl<-spMatrix(1,(length(X_ijkl_names)), i=rep(1, (length(X_ijkl_names))), j=(1:(length(X_ijkl_names))), x= c(rep(0, (length(X_ijkl_names)))))
                C15_Yia<-spMatrix(1,(length(Y_ia_names)), i=rep(1, (length(Y_ia_names))), j=(1:(length(Y_ia_names))), x= c(rep(0, (length(Y_ia_names)))))
                C15_Sia<-spMatrix(1,(length(S_ia_names)), i=rep(1, (length(S_ia_names))), j=(1:(length(S_ia_names))), x= c(rep(0, (length(S_ia_names)))))
                C15_Pia<-spMatrix(1,(length(P_ia_names)), i=rep(1, (length(P_ia_names))), j=(1:(length(P_ia_names))), x= c(rep(0, (length(P_ia_names)))))
                C15_Uij<-spMatrix(1,(length(u_ij_names)), i=rep(1, (length(u_ij_names))), j=(1:(length(u_ij_names))), x= c(rep(0, (length(u_ij_names)))))
                C15_teta<-spMatrix(1,(length(teta_ij_names)), i=rep(1, (length(teta_ij_names))), j=(1:(length(teta_ij_names))), x= c(rep(1, (length(teta_ij_names)))))
                C15<-cbind(C15_Xijkl,C15_Yia,C15_Sia,C15_Pia,C15_Uij,C15_teta)
                C15_rhs<-c(-log(1-wlimit))
                C15_dir<-c("<=")
                
      # BUILD MODEL
              
              model <- list()
              model$modelname <- 'UCAV Expected Damage'
              model$modelsense <- 'max'
              
              # initialize data for variables
              model$lb       <- c(rep(0, length(all_names)-length(u_ij_names)-length(teta_ij_names)), rep(1,length(u_ij_names) ),0)
              model$ub       <- c(rep(1, length(X_ijkl_names)),rep(1, length(Y_ia_names)),rep(max(PayloadTargetCompMtrx2[,7]), length(S_ia_names)),rep(1, length(P_ia_names)),rep((nNodes-1), length(u_ij_names)),1e20)
              model$vtype    <- c(rep("B", length(X_ijkl_names)),rep("B", length(Y_ia_names)), rep("C", length(S_ia_names)),rep("C", length(P_ia_names)),rep("C", length(u_ij_names)),rep("C", length(teta_ij_names))) #rep("B", length(wu_names)),"C","C","C","C","C")
              model$obj      <- all_obj
              model$varnames <- all_names
              
              # build constraint matrix
              model$A        <- rbind(C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15)
              model$rhs      <- c(C2_rhs,C3_rhs,C4_rhs,C5_rhs,C6_rhs,C7_rhs,C8_rhs,C9_rhs,C10_rhs,C11_rhs,C12_rhs,C13_rhs,C14_rhs,C15_rhs)
              model$sense    <- c(C2_dir,C3_dir,C4_dir,C5_dir,C6_dir,C7_dir,C8_dir,C9_dir,C10_dir,C11_dir,C12_dir,C13_dir,C14_dir,C15_dir)
              
        # gurobi_write(model, 'C:/Users/..../UCAV.mps')
              
              # set parameters
              params <- list()
              #params$TimeLimit<-10
              #params$MIPGap<-0
              #params$MIPFocus<-2
              
              print("rep:")
              print(rep)
              print("nTarget:")
              print(nTarget)
              print("nAppLoc:")
              print(nAppLoc)
              print("T:")
              print(MissionTimeLimit)
              
              ptm<-proc.time()
              result <- gurobi(model,params)
              
              runTime<-(proc.time()-ptm)
              runTime<-sum(runTime[1:2])

              decvr<-result$x
              all_names[which(decvr>0.0001)]
              decvr[which(decvr>0.0001)]

              # x_ijkl results
              x_ijkl_results<-decvr[1:length(X_ijkl_names)]
              x_ijkl_variables<-all_names[1:length(X_ijkl_names)]
              
              x_ijkl_variables<-x_ijkl_variables[which(x_ijkl_results>0.0001)]
              x_ijkl_results<-x_ijkl_results[which(x_ijkl_results>0.0001)]
              
              #yia results
              yia_results<-decvr[(length(X_ijkl_names)+1):(length(X_ijkl_names)+length(Y_ia_names))]
              yia_variables<-all_names[(length(X_ijkl_names)+1):(length(X_ijkl_names)+length(Y_ia_names))]
              # 
              yia_variables<-yia_variables[which(yia_results>0.0001)]
              yia_results<-yia_results[which(yia_results>0.0001)]
              
              #Sia results
              sia_results<-decvr[(length(X_ijkl_names)+length(Y_ia_names)+1):(length(X_ijkl_names)+length(Y_ia_names)+length(S_ia_names))]
              sia_variables<-all_names[(length(X_ijkl_names)+length(Y_ia_names)+1):(length(X_ijkl_names)+length(Y_ia_names)+length(S_ia_names))]
              # 
              sia_variables<-sia_variables[which(sia_results>0.0001)]
              sia_results<-sia_results[which(sia_results>0.0001)]
              
              #pia results
              pia_results<-decvr[(length(X_ijkl_names)+length(Y_ia_names)+length(S_ia_names)+1):(length(X_ijkl_names)+length(Y_ia_names)+length(S_ia_names)+length(P_ia_names))]
              pia_variables<-all_names[(length(X_ijkl_names)+length(Y_ia_names)+length(S_ia_names)+1):(length(X_ijkl_names)+length(Y_ia_names)+length(S_ia_names)+length(P_ia_names))]
              # 
              pia_variables<-pia_variables[which(pia_results>0.0001)]
              pia_results<-pia_results[which(pia_results>0.0001)]
              
              #teta_result
              teta_result<-tail(decvr,1)
              
              plot(TargetAppLoc[,,1],TargetAppLoc[,,2])
              for(i in 1:(nTarget+1)){
              text(TargetID_Type_Loc[i,3],TargetID_Type_Loc[i,4]-3*R,labels=(TargetID_Type_Loc[i,1]), col="red")
              }

              #Total Time of the Tour
              TotalTimeTour<-sum((X_ijkl_matrix[,5]/FlightV)*decvr[1:length(X_ijkl_names)])+ sum((R-sia_results)*2/FlightV)
              TotalRadarThreat<-teta_result
              
              #ResultReports
              
              sink("vPart2_Demonstration_Summary_1.txt", append=TRUE)
              cat(c(nrep, nTarget, nAppLoc, MissionTimeLimit, TotalTimeTour, wlimit, UAV_Type, TotalRadarThreat, result$status, result$objval, result$objbound, result$mipgap, result$runtime, runTime, c(TargetID_Type_Loc)))
              cat("\n")
              sink()
              
              sink("vPart2_Demonstration_DecSpace_1.txt", append=TRUE)
              cat(c(nrep, nTarget, nAppLoc, MissionTimeLimit,TotalTimeTour, wlimit, UAV_Type, result$status,(result$mipgap)*100,result$runtime,runTime,result$objval,teta_result,result$objbound,TotalTimeTour, c(TargetID_Type_Loc), x_ijkl_variables,yia_variables,sia_results,sia_variables))
              cat("\n")
              sink()
              