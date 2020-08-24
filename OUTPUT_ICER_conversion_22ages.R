
apportion<-function(small.set)
{
  
  if(is.matrix(small.set))
  {output<-cbind(
    small.set[,1], #0-<1 yrs LR
    small.set[,2], #1-2 yrs LR
    small.set[,2], #2-4 yrs LR
    small.set[,3], #5-11 yrs LR
    small.set[,3], #12-14 yrs LR
    small.set[,4], #15-16 yrs LR
    small.set[,4], #17-24 yrs LR
    small.set[,4], #25-44 yrs LR
    small.set[,5], #45-64 yrs LR
    small.set[,6], #65-74 LR
    small.set[,6], #75+ LR
    small.set[,7], #0-<1 yrs HR
    small.set[,8], #1-2 yrs HR
    small.set[,8], #2-4 yrs HR
    small.set[,9], #5-11 yrs HR
    small.set[,9], #12-14 yrs HR
    small.set[,10], #15-16 yrs HR
    small.set[,10], #17-24 yrs HR
    small.set[,10], #25-44 yrs HR
    small.set[,11], #45-64 yrs HR
    small.set[,12], #65-74 HR
    small.set[,12] #75+ HR
  )
  ;
  colnames(output)<-c('0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR')
  }
  
  if(is.vector(small.set)) {
    output<-c(
      small.set[1], #0-<1 yrs LR
      small.set[2], #1-2 yrs LR
      small.set[2], #2-4 yrs LR
      small.set[3], #5-11 yrs LR
      small.set[3], #12-14 yrs LR
      small.set[4], #15-16 yrs LR
      small.set[4], #17-24 yrs LR
      small.set[4], #25-44 yrs LR
      small.set[5], #45-64 yrs LR
      small.set[6], #65-74 LR
      small.set[6], #75+ LR
      small.set[7], #0-<1 yrs HR
      small.set[8], #1-2 yrs HR
      small.set[8], #2-4 yrs HR
      small.set[9], #5-11 yrs HR
      small.set[9], #12-14 yrs HR
      small.set[10], #15-16 yrs HR
      small.set[10], #17-24 yrs HR
      small.set[10], #25-44 yrs HR
      small.set[11], #45-64 yrs HR
      small.set[12], #65-74 HR
      small.set[12])#75+ HR
    
    names(output)<-c('0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR')
  }
  #Computes the number of health outcome from incidence table by sampling among the 
  #1000 samples of the risk table
  
  return(output)
}



CONVERT_INCIDENCE_TO_HEALTH_OUTCOME<<-function(incidence1,incidence2,outcome,strain)
{
  #checking list length
  inc.length<-length(incidence1)
  missing<-length(incidence1[1][sapply(incidence1,is.null)]) #number of missing
  tab.diff<<-inc.length-missing
  tab.temp1<-incidence1[1:tab.diff]
  tab.temp2<-incidence2[1:tab.diff]
  
  #strain<-c('H1N1','H3N2','B')
  setwd('/Users/Natasha/Dropbox/UKfluworkGIT/DavidsCode/RSave/')
  if(strain==1) {load(death.risk.tables[2]);load(hosp.risk.tables[2]);load(GP.risk.tables[2])}
  if(strain==2) {load(death.risk.tables[3]);load(hosp.risk.tables[3]);load(GP.risk.tables[3])}
  if(strain==3) {load(death.risk.tables[1]);load(hosp.risk.tables[1]);load(GP.risk.tables[1])}
  
  #Computes the number of symptomatic/febrile cases from incidence table
  if(outcome=="cases")
  {
    risk.sample<-rtriangle(n.samples,a=0.309,b=0.513,c=0.396) 
    risk.sample.ages<-matrix(rep(risk.sample, 22), ncol = 22, byrow = F)
    #vector not age stratified.
    #generate the percentile of febrile cases from triangular distribution
  } else { 
    #other outcomes are calculated using a risk table saved in RSave
    risk.tab <- loadRData(paste0(paste('tab_risk', outcome, strain.name[strain],sep="_"),'.R'))
    if(dim(incidence1[[1]])[1]==2500)
    {risk.sample<-rbind(risk.tab[sample.int(n.samples/2.5, n.samples/2.5),], risk.tab[sample.int(n.samples/2.5, n.samples/2.5),], risk.tab[sample.int(n.samples/5, n.samples/5),])}
    if(dim(incidence1[[1]])[1]==3000)
    {risk.sample<-rbind(risk.tab[sample.int(n.samples/3, n.samples/3),], risk.tab[sample.int(n.samples/3, n.samples/3),], risk.tab[sample.int(n.samples/3, n.samples/3),])
      }else{
    risk.sample<-rbind(risk.tab[sample.int(n.samples/2, n.samples/2),], risk.tab[sample.int(n.samples/2, n.samples/2),])}
    risk.sample.ages<-apportion(risk.sample)
  }
  
  risk1<-lapply(1:tab.diff, FUN = function(i) tab.temp1[[i]][1:n.samples,]*risk.sample.ages)
  risk2<-lapply(1:tab.diff, FUN = function(i) tab.temp2[[i]][1:n.samples,]*risk.sample.ages)
  return(list(risk1, risk2))
}

CONVERT_INCIDENCE_TO_QALY<<-function(incidence1,incidence2,strain,discount)
{
  #QALY loss from febrile cases
  QALY.loss.AJ<-read.csv(paste('/Users/Natasha/Dropbox/UKfluworkGIT/DavidsCode/',"DATA/AJ_QALY_list.csv",sep="/"))
  QALY.loss.cases<-sample(QALY.loss.AJ$v,n.samples,replace=TRUE)
  
  #QALY loss from hospitalisations, using distributions from Marc's paper
  QALY.loss.hosp<-rnorm(n.samples,mean=0.018,sd=0.0018); rownames(QALY.loss.hosp)<-NULL
  
  #Utilize conversion function. Calculate the different health outcomes
  table.cases<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="cases",strain); rownames(table.cases)<-NULL
  table.hosp<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="hosp",strain); rownames(table.hosp)<-NULL
  table.death<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="death",strain); rownames(table.death)<-NULL
  
  if(discount==3.5) death.QALY.loss<-apportion(c(22.61,23.00,22.32,18.84,12.71,6.52,24.07,24.05,23.18,20.33,14.62,7.28))
  if(discount==1.5) death.QALY.loss<-apportion(c(37.06,37.61,35.19,27.16,16.09,7.41,40.77,40.62,37.84,30.40,19.10,8.40))
  if(discount==0) death.QALY.loss<-apportion(c(60.93,61.62,54.96,38.23,19.75,8.24,69.96,69.45,61.64,44.47,24.14,9.46))
  
  QALY.loss1<-NULL
  #For Status QUo vaccine program
  QALY.cases1<-lapply(1:tab.diff, FUN = function(i) sweep(table.cases[[1]][[i]],MARGIN=1,QALY.loss.cases,'*'))
  QALY.hosp1<-lapply(1:tab.diff, FUN = function(i) sweep(table.hosp[[1]][[i]],MARGIN=1,QALY.loss.hosp,'*'))
  QALY.death1<-lapply(1:tab.diff, FUN = function(i) sweep(table.death[[1]][[i]],MARGIN=2,death.QALY.loss,'*'))
  
  
  QALY.loss2<-NULL
  #For Status QUo vaccine program
  QALY.cases2<-lapply(1:tab.diff, FUN = function(i) sweep(table.cases[[2]][[i]],MARGIN=1,QALY.loss.cases,'*'))
  QALY.hosp2<-lapply(1:tab.diff, FUN = function(i) sweep(table.hosp[[2]][[i]],MARGIN=1,QALY.loss.hosp,'*'))
  QALY.death2<-lapply(1:tab.diff, FUN = function(i) sweep(table.death[[2]][[i]],MARGIN=2,death.QALY.loss,'*'))
  
  QALY.loss1<-Map('+',Map('+',QALY.cases1,QALY.hosp1),QALY.death1)
  QALY.loss2<-Map('+',Map('+',QALY.cases2,QALY.hosp2),QALY.death2)
  
  return(list(QALY.loss1,QALY.loss2,QALY.cases1,QALY.cases2,QALY.hosp1,QALY.hosp2,QALY.death1,QALY.death2))
}

QALYdifferences <<- function(incidence1,incidence2,strain,discount)
{
  temp.sample.QALY <- CONVERT_INCIDENCE_TO_QALY(incidence1,incidence2,strain,discount);
  return(Map('-',temp.sample.QALY[[1]],temp.sample.QALY[[2]])) #should return number of QALY's gained in new program
}

GPHospAvert<<- function(incidence1,incidence2,strain){
  
  table.GP<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="GP",strain=strain)
  table.hosp<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="hosp",strain=strain)
  table.death<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="death",strain=strain)
  table.cases<-CONVERT_INCIDENCE_TO_HEALTH_OUTCOME(incidence1,incidence2,outcome="cases",strain=strain)
  
  #GP cases averted by new program
  Diff.GP <-mapply('-',table.GP[[1]],table.GP[[2]],SIMPLIFY=FALSE);
  
  #hospital cases averted by new program
  Diff.Hosp <- mapply('-',table.hosp[[1]],table.hosp[[2]],SIMPLIFY=FALSE);
  
  #formatted in a list with 2 pieces, first is status quo and second is intervention
  avg.GP1<-table.GP[[1]];
  avg.GP2<-table.GP[[2]];
  avg.hosp1<-table.hosp[[1]];
  avg.hosp2<-table.hosp[[2]];
  avg.death1<-table.death[[1]];
  avg.death2<-table.death[[2]];
  avg.cases1<-table.cases[[1]];
  avg.cases2<-table.cases[[2]];
  return(list(Diff.GP,Diff.Hosp,avg.GP1,avg.GP2,avg.hosp1,avg.hosp2, avg.death1,avg.death2, avg.cases1, avg.cases2))
}