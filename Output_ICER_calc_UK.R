library(fluEvidenceSynthesis)
library(plyr)
library(pander)
library(data.table)
library(triangle)
library(car)
library(proto)
library(ggplot2)
library(MASS)
library(colorRamps)
library(RColorBrewer)
library(tableone)
library(dplyr)
library(tidyr)
library(broom)
library(triangle)
library(matrixStats)
library(ggforce)
library(ks)
library(pryr)
library(ggthemes)
rm(list = ls())

#install 
##################################################
####### FILE LOADING FUNCTIONS FOR GRAPHS ########
#################################################

#load workspace
load("~/Dropbox/UKfluworkGIT/Outputfunctions/UK Files/UK_general_parameters.RData")


#loaders for all files
#source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/INTL_SENS files/SUPPORT_ICER_loaders.R')
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/UK Files/OUTPUT_ICER_sampler_UK.R')
vstrategy<-dget('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/UK Files/FUNC_cov_strategy_uk.R')
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens/Output_vaccine_doses_INTL.R')
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/UK Files/OUTPUT_ICER_conversion_22ages.R')

#load all conversion functions
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/UK Files/OUTPUT_icer_conversion_uk.R')
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/UK Files/FUNC_odin_UK.R')

#source('~/Dropbox/UKfluworkGIT/Outputfunctions/UK Files/Output_CEA_outcomes_UK.R')
#load workspace

#remainder items
y.name<<-c(1995:2014)

setwd('/Users/Natasha/Dropbox/UKfluworkGIT/DavidsCode/RSave/')
death.risk.tables<-list.files(pattern=glob2rx('tab_risk_death_*')); #loads in the order B,H1,H3
hosp.risk.tables<-list.files(pattern=glob2rx('tab_risk_hosp*'));
GP.risk.tables<-list.files(pattern=glob2rx('tab_risk_GP*'))

setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputfunctions', 'UK Files'));

####INTL FUNCTION Start=================================================
strainpull<-1
country<-4
i.country<-country
#setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))

#seventy<-c(coverageB[[18]]$V7, 0.7)
#seventy2<-c(seventy[4:length(seventy)], rep(0.7,3))
#seventy2[seventy2>0.7]<-0.7
#seventy.short<-seventy[4:24][seq(1,dim(coverageB[[18]])[1],4)]
#seventy.cov<-list(seventy.short, seventy.short,seventy.short,seventy.short,seventy.short,seventy.short,seventy.short,seventy.short,seventy.short,seventy.short,seventy.short,seventy.short,seventy.short,seventy.short, seventy2, seventy2, seventy2, seventy2, seventy2)

#thirty<-c(coverageB[[18]]$V9, 0.3)
#thirty2<-c(thirty[4:length(thirty)], rep(0.3,3))
#thirty2[thirty2 >0.3]<-0.3
#thirty.short<-thirty2[1:24][seq(1,dim(coverageB[[18]])[1],4)]
#thirty.cov<-list(thirty.short, thirty.short,thirty.short,thirty.short,thirty.short,thirty.short,thirty.short,thirty.short,thirty.short,thirty.short,thirty.short,thirty.short,thirty.short,thirty.short, thirty2, thirty2, thirty2, thirty2, thirty2)

#fiftyfive<-c(coverageB[[18]]$V12[4:23],0.5500)
#fiftyfive2<-c(fiftyfive, rep(0.55,3))
#fiftyfive2[fiftyfive2>0.55]<-0.55
#fiftyfive.short<-fiftyfive2[1:24][seq(1,dim(coverageB[[18]])[1],4)]
#fiftyfive.cov<-list(fiftyfive.short, fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short,fiftyfive.short, fiftyfive2, fiftyfive2, fiftyfive2, fiftyfive2, fiftyfive2)


setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/UK Files")
library(drc)
load('cov.function55.RData')
load('cov.function70.RData')
load('cov.function30.RData')

#SENS.samp.pars.uk(strain=2, num.samp = n.samples, version='u1')

#program 1 is status quo
#program 4 is secondary school vaccination only
#program 6 is prisec and preschool vaccination (all 3)
n.samples<-1000
#strainpull<-2

for(strainpull in 1:length(strain.name))
{
  for(strategy in 1:8)
  {
  #end vaccine program loop
  singleSENS.samp.uk(strain=strainpull, i.country = country, program=strategy,num.samp = n.samples, new.cov=thirty.cov,version='u1');
  singleSENS.samp.uk(strain=strainpull, i.country = country, program=strategy,num.samp = n.samples, new.cov=fiftyfive.cov, version='u1');
  singleSENS.samp.uk(strain=strainpull, i.country = country, program=strategy,num.samp = n.samples, new.cov=seventy.cov, version='u1');
  }
}


#############################################################################################
####### 6 QALY Conversion Equations ########
############################################################################################
strain<-strainpull
i.country<-4
setwd('/Users/Natasha/Dropbox/UKfluworkGIT/DavidsCode/RSave/')
if(strain==1) {load(death.risk.tables[2]);load(hosp.risk.tables[2]);load(GP.risk.tables[2])}
if(strain==2) {load(death.risk.tables[3]);load(hosp.risk.tables[3]);load(GP.risk.tables[3])}
if(strain==3) {load(death.risk.tables[1]);load(hosp.risk.tables[1]);load(GP.risk.tables[1])}
setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))



cov.sens.loader<-function(strain, end.cov, version)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB', 'coverage', end.cov))
  
  #Rearrange tables to correct order of interventions------------------
  sq.tables<-list.files(pattern=glob2rx(paste0('GB','SensStratStatusQuo', end.cov, strain.name[strain], version)))
  sens.tables<-list.files(pattern=glob2rx(paste0('GB','SensStrat*','*chool', end.cov, strain.name[strain], version)))
  sens.tables<-c(sens.tables, list.files(pattern=glob2rx(paste0('GB','SensStrat*','*ary', end.cov, strain.name[strain], version))))
  
  
  load(sq.tables[1]); avg.incidence1<-total.keep; #status quo
  load(sens.tables[2]); avg.incidence2<-total.keep; #preschool
  load(sens.tables[3]); avg.incidence3<-total.keep; #primary
  load(sens.tables[4]); avg.incidence4<-total.keep; #secondary
  load(sens.tables[1]); avg.incidence5<-total.keep; #secondary
  load(sens.tables[5]); avg.incidence6<-total.keep; #secondary
  load(sens.tables[6]); avg.incidence7<-total.keep; #secondary
  load(sens.tables[7]); avg.incidence8<-total.keep; #secondary
 

  var.name<-paste0('uk.table',strain.name[strain])
  #('WTP','Preschool','Primary','Secondary','Preschool+Primary','PrePrimeSeconday','Preschool+Secondary','Primary+Secondary')
  assign(var.name, list(avg.incidence1,avg.incidence2,avg.incidence3, avg.incidence4, avg.incidence5,
                        avg.incidence6, avg.incidence7,avg.incidence8
                        ), envir = .GlobalEnv)

}

##################################
####### 7 CALCULATE ICERs ########
#################################

####### 7.1 Single strain ICER calculation ########
Dataset1<-NULL
Dataset2<-NULL
Dataset3<-NULL
Dataset4<-NULL

#cov.level<-list(thirty.cov, fiftyfive.cov, seventy.cov)
cov.level<-list(mL3, mL55, mL7)

dcount.l<-c(0,1.5,3.5)
cov.f<-c(0.3, 0.55, 0.7)
for(cov.var in 1:length(cov.level))
{
GP.costs.l<-Hosp.costs.l<-sum.costs.l<-aged.cases.l<-QALY.tot.l<-QALY.death.l<-cumi.l<-add.dose<-avt.dose<-list()

for(dd in 1:length(dcount.l))
{
for(new.strat in 2:8)
  {
  cov.sens.loader(strain=2, end.cov=cov.f[cov.var], version='v1')
  
  uk.AgeStuc.CEA.outcomes(Dataset1 = uk.tableH3N2[[1]], Dataset2 = uk.tableH3N2[[new.strat]], strain = 2, dcount = dd, program1=1, program2=new.strat, i.cov=2, new.cov=cov.level[[cov.var]], version='v1')
}
}


for(dd in 1:length(dcount.l))
{
  for(new.strat in 2:4)
  {
    cov.sens.loader(strain=1, end.cov=cov.f[cov.var], version='v1')
    
    uk.AgeStuc.CEA.outcomes(Dataset1 = uk.tableH1N1[[1]], Dataset2 = uk.tableH1N1[[new.strat]], strain = 1, dcount = dd, program1=1, program2=new.strat,i.cov=1, new.cov=cov.level[[cov.var]], version='v1')
  }
}

for(dd in 1:length(dcount.l))
{
  for(new.strat in 2:8)
  {
    cov.sens.loader(strain=3, end.cov=cov.f[cov.var], version='v1')
    
    uk.AgeStuc.CEA.outcomes(Dataset1 = uk.tableB[[1]], Dataset2 = uk.tableB[[new.strat]], strain = 3, dcount = dd, program1=1, program2=new.strat, i.cov=3, new.cov=cov.level[[cov.var]], version='v1')
  }
}
}

##################################################################
########For 55% coverage
##################################################################

inv.names<-c("I1: 2-4 yrs old"  ,                 "I2: 5-11 yrs old"       ,              "I3: 12-16 yrs old"      ,             "I4: 2-11 yrs old" ,"I7: 2-16 yrs old", "I5: 2-4 & 12-16 yrs old"   ,      "I6: 5-16 yrs old")
tab.diff<-19
library(scales)
cols<-c("#FF61CC", "#F8766D","#CD9600","#7CAE00", "#00BE67","#C77CFF","#00BFC4","#00A9FF" )


#####################Cost-Effectiveness Plane
wtp<-20000
  
cep.graph <-function(strain, cov.var, dcount){  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB', 'coverage',cov.f[cov.var]))
  cea.tables<-list.files(pattern=glob2rx(paste0('GBICER.strat', strain.name[strain], cov.f[cov.var], 'p', '*', 'd',dcount.l[dcount])))
  
  cos.f<-qds.f<-qds.m<-cos.m<-cep<-cep1<-contour_95<-mean.ff<-mean.po<-not.ce<-NULL
  
  for(ff in 1:length(inv.names))
  {#ff=interventions. Strains are already chosen
    contour_95<-NULL
    load(cea.tables[[ff]]) #H1N1
    
    cer<-ICER.set2 #grab ICERs
    
    qa<-Reduce('+', cer$`QALY diff`)/tab.diff
    co<-Reduce('+', cer$`cost diff`)/tab.diff
    
    m.qa<-rowSums(qa)
    m.co<-rowSums(co)
    
    
    qds<-data.frame(cbind(rep(fname[i.country], dim(qa)[1]), rep(inv.names[[ff]], dim(qa)[1]), rep(dcount.l[dcount], dim(qa)[1]), qa, m.qa))
    cos<-data.frame(cbind(rep(fname[i.country], dim(qa)[1]) ,rep(inv.names[[ff]], dim(co)[1]), rep(dcount.l[dcount], dim(co)[1]), co, m.co))
    
    colnames(qds)<-c('country','program','discount','0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR','total')
    
    colnames(qds)<-colnames(cos)
    
    ifelse(ff==1, qds.f<-qds, qds.f<-rbind(qds.f, qds))
    ifelse(ff==1, cos.f<-cos, cos.f<-rbind(cos.f, cos))
    
    xy.points<-data.frame(m.qa, m.co)
    names(xy.points)<-c('x', 'y')
    kd <- ks::kde(xy.points, compute.cont=TRUE)
    contour_95 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                        z=estimate, levels=cont["10%"])[[1]])
    contour_95 <- data.frame(contour_95, rep(inv.names[[ff]], length(contour_95$x)))
    colnames(contour_95)<-c('level','x','y','program')
    
    p1<-geom_path(data=contour_95, aes(x, y, color=program))
    contour.name<-paste0('cont',strain, cov.var, dcount.l[dcount], ff)
    assign(contour.name,  p1)
  }

colnames(qds.f)<-c('country','Intervention','discount','0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR','total')
qds.i<-melt(qds.f, id.vars = c('country', 'Intervention','discount'))
colnames(cos.f)<-colnames(qds.f)
cos.i<-melt(cos.f, id.vars = c('country', 'Intervention','discount'))
#setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')

llg<-cbind(qds.i[,1:4],as.numeric(as.character(qds.i$value)), as.numeric(as.character(cos.i$value))) 
colnames(llg)<-c('country','Intervention','discount', 'Age','qaly','cost')
llg<-na.omit(llg)

totals<-llg[llg$Age=='total',]
#library('plyr')
#mean.p<-ddply(totals, .(Intervention), summarize, Rate1=median(qaly), Rate2=median(cost))
#mean.po<-mean.p[order(mean.p$Rate2, decreasing = T),]
##   continent n_uniq_countries

mean.p<-ddply(as.data.frame(llg), .(Intervention), summarize, Rate1=median(as.numeric(as.character(qaly))), Rate2=median(as.numeric(as.character(cost))), costLB=quantile(as.numeric(as.character(cost)), 0.025), costUB=quantile(as.numeric(as.character(cost)), 0.975), qalyLB=quantile(as.numeric(as.character(qaly)), 0.025), qalyUB=quantile(as.numeric(as.character(qaly)), 0.975)   )

mean.sd<-ddply(as.data.frame(llg), .(Intervention), summarize, cost.sd=sd(as.numeric(as.character(cost))), qal.sd=sd(as.numeric(as.character(qaly))))

mean.po<-mean.p[order(mean.p$Rate2, decreasing = T),]

not.ce<-not.ce1<-not.ce2<-not.ce3<-not.ce4<-counter<-NULL

rudi.icer<-diff(rbind(c(0,0),cbind(rev(mean.po[,2]), rev(mean.po[,3]))))
rudi.icer.ub<-diff(rbind(c(0,0),cbind(rev(mean.po$costUB), rev(mean.po$qalyUB))))
rudi.icer.lb<-diff(rbind(c(0,0),cbind(rev(mean.po$costLB), rev(mean.po$qalyLB))))
rudi.icer.sd<-diff(rbind(c(0,0),cbind(rev(mean.sd$cost.sd), rev(mean.sd$qal.sd))))

div<-rudi.icer[,2]/rudi.icer[,1]
div.ub<-rudi.icer.ub[,1]/rudi.icer.ub[,2]
div.lb<-rudi.icer.lb[,1]/rudi.icer.lb[,2]
div.sd<-rudi.icer.sd[,1]/rudi.icer.sd[,2]


mean.f<-data.frame(cbind(apply(mean.po[which(div >0&div <20000),], 2, rev), round(div[which(div >0&div <20000)]), round(div.ub[which(div >0&div <20000)]), round(div.lb[which(div >0&div <20000)])))
levels(mean.f$Intervention)= c(levels(mean.f$Intervention),'SQ: Risk Groups, 65+ yrs old')
sq<-matrix(c('SQ: Risk Groups, 65+ yrs old',rep(0,9)), nrow=1)
mean.ff<-data.frame(rbind(sq, as.matrix(mean.f[order(mean.f$Rate2, decreasing=F),])))

not.ce<-mean.po[which(div <0 | div >20000),]
colnames(mean.ff)<-colnames(mean.f)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB', 'icer', 'icerUB', 'icerLB')

if(length(not.ce)!=0)
{
  counter<-2
  ##second icer eval
  rudi.icer2<-diff(cbind(rev(as.numeric(as.character(mean.ff[,2]))), rev(as.numeric(as.character(mean.ff[,3])))))
  rudi.icer.ub2<-diff(cbind(rev(as.numeric(as.character(mean.ff$costUB))), rev(as.numeric(as.character(mean.ff$qalyUB)))))
  rudi.icer.lb2<-diff(cbind(rev(as.numeric(as.character(mean.ff$costLB))), rev(as.numeric(as.character(mean.ff$qalyLB)))))
  
  
  div2<-rudi.icer2[,2]/rudi.icer2[,1]
  div.ub2<-rudi.icer.ub2[,1]/rudi.icer.ub2[,2]
  div.lb2<-rudi.icer.lb2[,1]/rudi.icer.lb2[,2]
  
  mean.f2<-data.frame(cbind(apply(mean.f[which(div2 >0&div2 <20000),],2,rev), round(div2[which(div2 >0&div2 <20000)]), round(div.ub2[which(div2 >0&div2 <20000)]), round(div.lb2[which(div2 >0&div2 <20000)])))
  
  levels(mean.f2$Intervention)= c("I1: 2-4 yrs old"  ,                 "I2: 5-11 yrs old"       ,              "I3: 12-16 yrs old"      ,             "I4: 2-11 yrs old" ,"I7: 2-16 yrs old", "I5: 2-4 & 12-16 yrs old"   ,      "I6: 5-16 yrs old" , 'SQ: Risk Groups, 65+ yrs old') 
  mean.ff2<-data.frame(rbind(sq, as.matrix(mean.f2[order(mean.f2$Rate2, decreasing=F),])))
  
  not.ce2<-rbind(not.ce, mean.f[which(rudi.icer2[,2]/rudi.icer2[,1] <0 | rudi.icer2[,2]/rudi.icer2[,1] >20000),])
  colnames(mean.ff2)<-c('Intervention','Rate1','Rate2','icer')
  }

if(any(dim(not.ce2)!=dim(not.ce)))
{
  counter<-counter+1
  ##3rd CE check
  #mean.po<-mean.p[order(mean.p$Rate2, decreasing = T),]
  rudi.icer3<-diff(cbind(rev(as.numeric(as.character(mean.ff2[,2]))), rev(as.numeric(as.character(mean.ff2[,3])))))
  div3<-rudi.icer3[,2]/rudi.icer3[,1]
  mean.f3<-cbind(mean.f2[which(div3 >0&div3 <20000),][,1:3], div3[which(div3 >0&div3 <20000)])
  levels(mean.f3$Intervention)= c("I1: 2-4 yrs old"  ,                 "I2: 5-11 yrs old"       ,     "I3: 12-16 yrs old", "I4: 2-11 yrs old" ,"I7: 2-16 yrs old", "I5: 2-4 & 12-16 yrs old"   ,      "I6: 5-16 yrs old" , 'SQ: Risk Groups, 65+ yrs old') 
  mean.ff3<-data.frame(rbind(sq, as.matrix(mean.f3[order(mean.f3$Rate2, decreasing=F),])))
  
  not.ce3<-rbind(not.ce2, mean.f[which(rudi.icer3[,2]/rudi.icer3[,1] <0 | rudi.icer3[,2]/rudi.icer3[,1] >20000),])
  colnames(mean.ff3)<-c('Intervention','Rate1','Rate2','icer')
}

if(any(dim(not.ce3)!=dim(not.ce2)))
{
  counter<-counter+1
  ##4th CE check
  #mean.po<-mean.p[order(mean.p$Rate2, decreasing = T),]
  rudi.icer4<-diff(cbind(rev(as.numeric(as.character(mean.ff3[,2]))), rev(as.numeric(as.character(mean.ff3[,3])))))
  div3<-rudi.icer4[,2]/rudi.icer4[,1]
  mean.f4<-cbind(mean.f3[which(div4>0&div4 <20000),][,1:3], div4[which(div4 >0&div4 <20000)])
  levels(mean.f4$Intervention)= c("I1: 2-4 yrs old"  ,                 "I2: 5-11 yrs old"       ,              "I3: 12-16 yrs old"      ,             "I4: 2-11 yrs old" ,"I7: 2-16 yrs old", "I5: 2-4 & 12-16 yrs old"   ,      "I6: 5-16 yrs old" , 'SQ: Risk Groups, 65+ yrs old') 
  mean.ff4<-data.frame(rbind(sq, as.matrix(mean.f4[order(mean.f4$Rate2, decreasing=F),])))
  
  not.ce4<-rbind(not.ce3, mean.f[which(rudi.icer4[,2]/rudi.icer4[,1] <0 | rudi.icer4[,2]/rudi.icer4[,1] >20000),])
  colnames(mean.ff4)<-c('Intervention','Rate1','Rate2','icer')
}

######CE evaluation

if(counter>1)
{
mean.ff<-get(paste0('mean.ff', counter))
not.ce<-get(paste0('not.ce', counter))
}


####
cep<-cep1<-NULL
#llg[llg$cost>1e6|llg$cost < (-1e6),]<-NA
cep<-ggplot(totals, aes(qaly, cost, color=Intervention))+
  scale_y_continuous(breaks=scales::pretty_breaks(n = 7))+scale_x_continuous(breaks=scales::pretty_breaks(n=6))+geom_hline(aes(yintercept=0))+geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+geom_abline(slope = 15000, intercept=0, linetype='dotdash')+geom_abline(slope = 20000, intercept=0, linetype='longdash')+labs(subtitle=paste(paste0('discount=',dcount.l[dcount],'%'), paste0('strain=',strain.name[strain]), paste0('coverage=',cov.f[cov.var])),x='QALYs gained', y='Net Costs (£)', fill='Intervention')+theme(plot.subtitle=element_text(size=10, hjust=0.5, color="black"))+guides(fill=guide_legend(title="Intervention"))


#+ylab('Incremental Cost in million ₤/year')+xlab('QALY Difference')+labs(subtitle=paste(paste0('discount=',dcount.l[dcount],'%'), paste0('strain=',strain.name[strain]), paste0('coverage=',cov.f[cov.var])), fill='Intervention')

for(fff in 1:length(inv.names))
{add.on<-get(paste0('cont',strain, cov.var, dcount.l[dcount], fff))
  cep<-cep+add.on}
 

       cep1<-cep+geom_point(data=mean.ff, aes(x=as.numeric(mean.ff$Rate1), y=as.numeric(mean.ff$Rate2)),size=3, shape=19, stroke=1, inherit.aes = T, show.legend=F)+geom_line(data=mean.ff, aes(x=as.numeric(mean.ff$Rate1), y=as.numeric(mean.ff$Rate2)), inherit.aes = F)+geom_point(data=not.ce, aes(x=not.ce$Rate1, not.ce$Rate2), size=3, shape=10, inherit.aes = T, show.legend=F)+labs(fill='Intervention')+theme_bw()
       #)


#setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB', 'coverage', cov.f[cov.var]))
#ggsave(paste0('CEP',strain.name[strain],cov.f[cov.var],dcount,".png"),plot = cep, width=11, height=8.5, units='in',device='png')
return(cep1)

}


cep.grapher.short<-function(strain,cov.var, dcount, wtp)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB', 'coverage',cov.f[cov.var]))
  cea.tables<-list.files(pattern=glob2rx(paste0('GBICER.strat', strain.name[strain],cov.f[cov.var],'p', '*', 'd',dcount)))
  i.country<-4
  cos.f<-qds.f<-qds.m<-cos.m<-NULL
  tab.diff<-14
  for(ff in 1:length(inv.names))
  {#ff=interventions. Strains are already chosen
    load(cea.tables[[ff]]) #H1N1
    
    cer<-ICER.set2 #grab ICERs
    
    #m.qa<-median(qa)
    #m.co<-median(co)
    
    qa1<-Reduce('+', cer$`net QALY`[1:14])/tab.diff
    co1<-Reduce('+', cer$`net cost`[1:14])/tab.diff
    qds<-data.frame(cbind(rep(fname[i.country], length(qa1)), qa1, rep(inv.names[[ff]], length(qa1))))
    cos<-data.frame(cbind(rep(fname[i.country], length(co1)), co1, rep(inv.names[[ff]], length(co1))))
    
    #qa<-unlist(cer$`net QALY`[1:14])
    #co<-unlist(cer$`net cost`[1:14])
    
    #qds<-data.frame(cbind(rep(fname[i.country], length(qa)), qa, rep(inv.names[[ff]], length(qa))))
    #cos<-data.frame(cbind(rep(fname[i.country], length(co)), co, rep(inv.names[[ff]], length(co))))
    
    colnames(qds)<-colnames(cos)
    
    ifelse(ff==1, qds.f<-qds, qds.f<-rbind(qds.f, qds))
    ifelse(ff==1, cos.f<-cos, cos.f<-rbind(cos.f, cos))
    
    meds.f<-cbind(as.numeric(as.character(qds[,2])), as.numeric(as.character(cos[,2])))
    
    xy.points<-data.frame(meds.f)
    names(xy.points)<-c('x', 'y')
    kd <- ks::kde(xy.points, compute.cont=TRUE)
    contour_95 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                        z=estimate, levels=cont["10%"])[[1]])
    contour_95$y[length(contour_95$y)]<-contour_95$y[1]
    contour_95 <- data.frame(contour_95, rep(inv.names[[ff]], length(contour_95$x)))
    colnames(contour_95)<-c('level','x','y','program')
    
    p1<-geom_path(data=contour_95, aes(x, y, color=program))
    contour.name<-paste0('cont',strain.name[strain],cov.var,dcount.l[dcount], ff)
    assign(contour.name,  p1)
  }
  
  
  colnames(qds.f)<-  c('country', 'net QALYS', 'Intervention')
  colnames(cos.f)<- c('country', 'net costs', 'Intervention')
  
  counter<-c()
  cost.utility<-as.matrix(cbind(qds.f[,1:2], cos.f[,2], qds.f[,3]))
  cost.utility<-data.frame(na.omit(unname(cost.utility)))
  colnames(cost.utility)<-c('country', 'QALYS', 'Costs', 'Intervention')
  #row.names(cost.utility) <- c()
  #attributes(cost.utility)$row.names<-NULL
  #mean.p<-ddply(cost.utility, .(country, Intervention), summarize, Rate1=mean(as.numeric(as.character(QALYS))), Rate2=mean(as.numeric(as.character(Costs))), sd.q=sd(as.numeric(as.character(QALYS))),sd.c=sd(as.numeric(as.character(Costs))))
  
  #mean.po<-mean.p[order(as.numeric(as.character(mean.p$Rate2)), decreasing = T),]
  ##   continent n_uniq_countries
  
  
  mean.p<-ddply(as.data.frame(cost.utility), .(Intervention), summarize, Rate1=mean(as.numeric(as.character(QALYS))), Rate2=mean(as.numeric(as.character(Costs))), costLB=quantile(as.numeric(as.character(Costs)), 0.025), costUB=quantile(as.numeric(as.character(Costs)), 0.975), qalyLB=quantile(as.numeric(as.character(QALYS)), 0.025), qalyUB=quantile(as.numeric(as.character(QALYS)), 0.975), cost.sd=sd(as.numeric(as.character(Costs))), qal.sd=sd(as.numeric(as.character(QALYS))))
  
  mean.po<-mean.p[order(mean.p$Rate2, decreasing = F),]
  
  not.ce<-not.ce1<-not.ce2<-not.ce3<-not.ce4<-counter<-NULL
  
  rudi.icer<-diff(rbind(c(0,0),cbind(mean.po$Rate1, mean.po$Rate2)))*-1
  rudi.icer.ub<-diff(rbind(c(0,0),cbind(mean.po$costUB, mean.po$qalyUB)))*-1
  rudi.icer.lb<-diff(rbind(c(0,0),cbind(mean.po$costLB, mean.po$qalyLB)))*-1
  rudi.icer.sd<-diff(rbind(c(0,0),cbind(mean.po$cost.sd, mean.po$qal.sd)))*-1
  
  div<-rudi.icer[,2]/rudi.icer[,1]
  div.ub<-rudi.icer.ub[,1]/rudi.icer.ub[,2]
  div.lb<-rudi.icer.lb[,1]/rudi.icer.lb[,2]
  div.sd<-rudi.icer.sd[,1]/rudi.icer.sd[,2]
  
  flip<-as.matrix(mean.po[which(div >0&div <wtp),])
  
  if(dim(flip)[1]<2)
  {
    mean.f<-data.frame(t(as.matrix(c(flip[1:9], 
                                     round(div[which(div >0&div <wtp)]), 
                                     round(div.ub[which(div >0&div <wtp)]), #quantile
                                     round(div.lb[which(div >0&div <wtp)]), #quantile
                                     round(div[which(div >0&div <wtp)]-1.96*div.sd[which(div >0&div <wtp)]), round(div[which(div >0&div <wtp)]+1.96*div.sd[which(div >0&div <wtp)]))))) #sd
  }else{
    mean.f<-data.frame(cbind(as.matrix(flip[,1:9]),
                             round(div[which(div >0&div <wtp)]), 
                             round(div.ub[which(div >0&div <wtp)]), #quantile
                             round(div.lb[which(div >0&div <wtp)]), #quantile
                             round(div[which(div >0&div <wtp)]-1.96*div.sd[which(div >0&div <wtp)]), round(div[which(div >0&div <wtp)]+1.96*div.sd[which(div >0&div <wtp)]))) #sd
  }  
  
  sq<-matrix(c('SQ: Risk Groups, 65+ yrs old',rep(0,13)), nrow=1)
  
  if(dim(mean.f)[2]>0)
  { colnames(mean.f)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
  levels(mean.f$Intervention)= c(levels(mean.f$Intervention),'SQ: Risk Groups, 65+ yrs old')
  mean.ff<-data.frame(as.matrix(mean.f[order(mean.f$Rate2, decreasing=F),]))
  colnames(mean.ff)<-colnames(mean.f)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
  }else{
    mean.ff<-NULL
  }
  
  not.ce<-mean.po[which(div <0 | div >wtp),]
  counter<-1
  
  if(dim(not.ce)[1]!=0)
  {
    counter<-counter+1
    
    ##second icer eval
    rudi.icer2<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff[,2])), as.numeric(as.character(mean.ff[,3])))))*-1
    rudi.icer.ub2<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff$costUB)), as.numeric(as.character(mean.ff$qalyUB)))))*-1
    rudi.icer.lb2<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff$costLB)), as.numeric(as.character(mean.ff$qalyLB)))))*-1
    rudi.icer.sd2<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff$cost.sd)), as.numeric(as.character(mean.ff$qal.sd)))))*-1
    #rudi.icer.sd<-diff(rev(as.numeric(as.character(mean.ff$sd))))
    
    
    div2<-rudi.icer2[,2]/rudi.icer2[,1]
    div.ub2<-rudi.icer.ub2[,1]/rudi.icer.ub2[,2]
    div.lb2<-rudi.icer.lb2[,1]/rudi.icer.lb2[,2]
    div.sd2<-rudi.icer.sd2[,1]/rudi.icer.sd2[,2]
    
    flip2<-as.matrix(mean.f[which(div2 >0&div2 <wtp),])
    
    if(dim(flip2)[1]<2)
    {
      
      mean.f2<-data.frame(t(as.matrix(c(flip2[1:9], 
                                        round(div2[which(div2 >0&div2 <wtp)]), 
                                        round(div.ub2[which(div2 >0&div2 <wtp)]), #quantile
                                        round(div.lb2[which(div2 >0&div2 <wtp)]), #quantile
                                        round(div2[which(div2 >0&div2 <wtp)]-1.96*div.sd2[which(div2 >0&div2 <wtp)]), round(div2[which(div2 >0&div2 <wtp)]+1.96*div.sd2[which(div2 >0&div2 <wtp)]))))) #sd
    }else{
      mean.f2<-data.frame(cbind(as.matrix(flip2[,1:9]), 
                                round(div2[which(div2 >0&div2 <wtp)]), 
                                round(div.ub2[which(div2 >0&div2 <wtp)]), #quantile
                                round(div.lb2[which(div2 >0&div2 <wtp)]), #quantile
                                round(div2[which(div2 >0&div2 <wtp)]-1.96*div.sd2[which(div2 >0&div2 <wtp)]), round(div2[which(div2 >0&div2 <wtp)]+1.96*div.sd2[which(div2 >0&div2 <wtp)]))) #sd
    }
    
    
    
    if(dim(mean.f2)[2]==14)
    {
      colnames(mean.f2)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
      levels(mean.f2$Intervention)= c(levels(mean.f2$Intervention),'SQ: Risk Groups, 65+ yrs old')
      mean.ff2<-data.frame(as.matrix(mean.f2[order(mean.f2$Rate2, decreasing=F),]))
      colnames(mean.ff2)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
    }else{
      mean.ff2<-NULL
    }
    
    not.ce2<-rbind(not.ce, mean.f[which(div2 <0 | div2 >wtp),1:9])
    
    
    
    
    if(dim(not.ce2)[1]!=3)
    {
      if(any(dim(not.ce2)!=dim(not.ce)))
      {
        counter<-counter+1
        ##3rd CE check
        
        rudi.icer3<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff2[,2])), as.numeric(as.character(mean.ff2[,3])))))*-1
        rudi.icer.ub3<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff2$costUB)), as.numeric(as.character(mean.ff2$qalyUB)))))*-1
        rudi.icer.lb3<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff2$costLB)), as.numeric(as.character(mean.ff2$qalyLB)))))*-1
        rudi.icer.sd3<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff2$cost.sd)), as.numeric(as.character(mean.ff2$qal.sd)))))*-1
        #rudi.icer.sd<-diff(rev(as.numeric(as.character(mean.ff$sd))))
        
        
        div3<-rudi.icer3[,2]/rudi.icer3[,1]
        div.ub3<-rudi.icer.ub3[,1]/rudi.icer.ub3[,2]
        div.lb3<-rudi.icer.lb3[,1]/rudi.icer.lb3[,2]
        div.sd3<-rudi.icer.sd3[,1]/rudi.icer.sd3[,2]
        
        flip3<-as.matrix(mean.f2[which(div3 >0&div3 <wtp),])
        
        mean.f3<-data.frame(cbind(as.matrix(flip3[,1:9]), 
                                  round(div3[which(div3 >0&div3 <wtp)]), 
                                  round(div.ub3[which(div3 >0&div3 <wtp)]), #quantile
                                  round(div.lb3[which(div3 >0&div3 <wtp)]), #quantile
                                  round(div3[which(div3 >0&div3 <wtp)]-1.96*div.sd3[which(div3 >0&div3 <wtp)]), round(div3[which(div3 >0&div3 <wtp)]+1.96*div.sd3[which(div3 >0&div3 <wtp)]))) #sd
        
        
        if(dim(mean.f3)[2]==14)
        { colnames(mean.f3)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
        levels(mean.f3$Intervention)= c(levels(mean.f3$Intervention),'SQ: Risk Groups, 65+ yrs old')
        mean.ff3<-data.frame(rbind(sq, as.matrix(mean.f3[order(mean.f2$Rate2, decreasing=F),])))
        colnames(mean.ff3)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
        }else
        {mean.ff3<-NULL}
        not.ce3<-rbind(not.ce2, mean.f2[which(div3 <0 | div3 >wtp),1:9])
        
        
        
        if(dim(not.ce3)[1]!=4)
        {
          if(any(dim(not.ce3)!=dim(not.ce2)))
          {
            counter<-counter+1
            ##3rd CE check
            
            rudi.icer4<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff3[,2])), as.numeric(as.character(mean.ff3[,3])))))*-1
            rudi.icer.ub4<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff3$costUB)), as.numeric(as.character(mean.ff3$qalyUB)))))*-1
            rudi.icer.lb4<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff3$costLB)), as.numeric(as.character(mean.ff3$qalyLB)))))*-1
            rudi.icer.sd4<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff3$cost.sd)), as.numeric(as.character(mean.ff3$qal.sd)))))*-1
            #rudi.icer.sd<-diff(rev(as.numeric(as.character(mean.ff$sd))))
            
            
            div4<-rudi.icer4[,2]/rudi.icer4[,1]
            div.ub4<-rudi.icer.ub4[,1]/rudi.icer.ub4[,2]
            div.lb4<-rudi.icer.lb4[,1]/rudi.icer.lb4[,2]
            div.sd4<-rudi.icer.sd4[,1]/rudi.icer.sd4[,2]
            
            flip4<-as.matrix(mean.f3[which(div4 >0&div4 <wtp),])
            
            #if(dimflip3[1]!=0)
            #{
            mean.f4<-data.frame(cbind(as.matrix(c(flip4[1:9], 
                                                  round(div4[which(div4 >0&div4 <wtp)]), 
                                                  round(div.ub4[which(div4 >0&div4 <wtp)]), #quantile
                                                  round(div.lb4[which(div4 >0&div4 <wtp)]), #quantile
                                                  round(div4[which(div4 >0&div4 <wtp)]-1.96*div.sd4[which(div4 >0&div4 <wtp)]), round(div4[which(div4 >0&div4 <wtp)]+1.96*div.sd4[which(div4 >0&div4 <wtp)]))))) #sd
            
            
            #}
            
            if(dim(mean.f4)[2]==14)
            { colnames(mean.f4)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
            levels(mean.f4$Intervention)= c(levels(mean.f4$Intervention),'SQ: Risk Groups, 65+ yrs old')
            mean.ff4<-data.frame(rbind(sq, as.matrix(mean.f4[order(mean.f3$Rate2, decreasing=F),])))
            colnames(mean.ff4)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
            }else
            {mean.ff4<-NULL}
            
            not.ce4<-rbind(not.ce3, mean.f3[which(div4 <0 | div4 >wtp),1:9])
            
            
            if(dim(not.ce4)[1]!=5)
            {
              if(any(dim(not.ce5)!=dim(not.ce4)))
              {
                counter<-counter+1
                ##3rd CE check
                
                rudi.icer5<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff4[,2])), as.numeric(as.character(mean.ff4[,3])))))*-1
                rudi.icer.ub5<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff4$costUB)), as.numeric(as.character(mean.ff4$qalyUB)))))*-1
                rudi.icer.lb5<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff4$costLB)), as.numeric(as.character(mean.ff4$qalyLB)))))*-1
                rudi.icer.sd5<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff4$cost.sd)), as.numeric(as.character(mean.ff4$qal.sd)))))*-1
                #rudi.icer.sd<-diff(rev(as.numeric(as.character(mean.ff$sd))))
                
                
                div5<-rudi.icer5[,2]/rudi.icer5[,1]
                div.ub5<-rudi.icer.ub5[,1]/rudi.icer.ub5[,2]
                div.lb5<-rudi.icer.lb5[,1]/rudi.icer.lb5[,2]
                div.sd5<-rudi.icer.sd5[,1]/rudi.icer.sd5[,2]
                
                flip5<-as.matrix(mean.f4[which(div5 >0&div5 <wtp),])
                
                #if(dimflip3[1]!=0)
                #{
                mean.f5<-data.frame(cbind(as.matrix(c(flip5[1:9], 
                                                      round(div5[which(div5 >0&div5 <wtp)]), 
                                                      round(div.ub5[which(div5 >0&div5 <wtp)]), #quantile
                                                      round(div.lb5[which(div5 >0&div5 <wtp)]), #quantile
                                                      round(div5[which(div5 >0&div5 <wtp)]-1.96*div.sd5[which(div5 >0&div5 <wtp)]), round(div5[which(div5 >0&div5 <wtp)]+1.96*div.sd5[which(div5 >0&div5 <wtp)]))))) #sd
                
                
                #}
                
                if(dim(mean.f5)[2]==14)
                { colnames(mean.f5)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
                levels(mean.f5$Intervention)= c(levels(mean.f5$Intervention),'SQ: Risk Groups, 65+ yrs old')
                mean.ff5<-data.frame(rbind(sq, as.matrix(mean.f5[order(mean.f4$Rate2, decreasing=F),])))
                colnames(mean.ff5)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
                }else
                {mean.ff5<-NULL}
                not.ce5<-rbind(not.ce4, mean.f4[which(div5 <0 | div5 >wtp),1:9])
                
                if(dim(not.ce5)[1]!=6)
                {
                  if(any(dim(not.ce5)!=dim(not.ce4)))
                  {
                    counter<-counter+1
                    ##3rd CE check
                    
                    rudi.icer6<-rbind(c(0,0),cbind(as.numeric(as.character(mean.ff5[,2])), as.numeric(as.character(mean.ff5[,3]))))*-1
                    rudi.icer.ub6<-rbind(c(0,0),cbind(as.numeric(as.character(mean.ff5$costUB)), as.numeric(as.character(mean.ff5$qalyUB))))*-1
                    rudi.icer.lb6<-rbind(c(0,0),cbind(as.numeric(as.character(mean.ff5$costLB)), as.numeric(as.character(mean.ff5$qalyLB))))*-1
                    rudi.icer.sd6<-rbind(c(0,0),cbind(as.numeric(as.character(mean.ff5$cost.sd)), as.numeric(as.character(mean.ff5$qal.sd))))*-1
                    #rudi.icer.sd<-diff(rev(as.numeric(as.character(mean.ff$sd))))
                    
                    
                    div6<-rudi.icer6[,2]/rudi.icer6[,1]
                    div.ub6<-rudi.icer.ub6[,1]/rudi.icer.ub6[,2]
                    div.lb6<-rudi.icer.lb6[,1]/rudi.icer.lb6[,2]
                    div.sd6<-rudi.icer.sd6[,1]/rudi.icer.sd6[,2]
                    
                    flip6<-as.matrix(mean.f5[which(div6 >0&div6 <wtp),])
                    
                    #if(dimflip3[1]!=0)
                    #{
                    mean.f6<-data.frame(t(as.matrix(c(flip6[1:9], 
                                                      round(div6[which(div6 >0&div6 <wtp)]), 
                                                      round(div.ub6[which(div6 >0&div6 <wtp)]), #quantile
                                                      round(div.lb6[which(div6 >0&div6 <wtp)]), #quantile
                                                      round(div6[which(div6 >0&div6 <wtp)]-1.96*div.sd6[which(div6 >0&div6 <wtp)]), round(div6[which(div6 >0&div6 <wtp)]+1.96*div.sd6[which(div6 >0&div6 <wtp)]))))) #sd
                    
                    
                    #}
                    
                    if(dim(mean.f6)[2]==14)
                    { colnames(mean.f6)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
                    levels(mean.f6$Intervention)= c(levels(mean.f6$Intervention),'SQ: Risk Groups, 65+ yrs old')
                    mean.ff6<-data.frame(rbind(sq, as.matrix(mean.f6[order(mean.f5$Rate2, decreasing=F),])))
                    colnames(mean.ff6)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
                    }else
                    {mean.ff6<-NULL}
                    not.ce6<-rbind(not.ce5, mean.f5[which(div6 <0 | div6 >wtp),1:9])}} 
              }
            }
          }
        }
      }
    }
  }
  
  ######CE evaluation
  
  if(counter>1)
  {
    mean.ff<-get(paste0('mean.ff', counter))
    not.ce<-get(paste0('not.ce', counter))
  }
  
  
  ####axis.text.x = element_text(face="bold", color="#993333", 
  #size=14
  cep<-cep1<-NULL
  cols<-c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC")
  
  cep<-ggplot(cost.utility, aes(as.numeric(as.character(QALYS)), as.numeric(as.character(Costs)), color=Intervention))+theme(axis.text.x = element_text(angle = 45, hjust = 1, size=11), plot.subtitle =element_text(angle = 45, hjust = 1, size=10) )+scale_y_continuous(breaks=scales::pretty_breaks(n = 7), limits = c(0,1.5e8))+scale_x_continuous(breaks=scales::pretty_breaks(n=4), limits=c(-2000,61000))+geom_hline(aes(yintercept=0))+geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+geom_abline(slope = 15000, intercept=0, linetype='dotdash')+geom_abline(slope = wtp, intercept=0, linetype='longdash')+labs(subtitle=paste(paste0('discount=',dcount.l[dcount],'%'), paste0(strain.name[strain]), paste0('coverage=',cov.f[cov.var]) ), fill='Intervention')+theme_bw()+scale_colour_manual(values=cols)
  
  
  for(fff in 1:length(inv.names))
  {add.on<-get(paste0('cont',strain.name[strain],cov.var,dcount.l[dcount], fff))
  cep<-cep+add.on}
  #cep.name<-paste0('cep',strain, cov.var, dcount.l[dcount])
  
  points<-cbind(c(0,as.numeric(as.character(mean.ff[,2]))), c(0,as.numeric(as.character(mean.ff[,3]))))
  names(points)<-c('Rate1', 'Rate2')
  
  if(is.null(mean.ff))
  {
    levels(not.ce$Intervention)=c(levels(not.ce$Intervention),'SQ: Risk Groups, 65+ yrs old')
    not.ce<-data.frame(rbind(sq[1:9], as.matrix(not.ce[order(not.ce$Rate2, decreasing=F),])))
    
    cep1<-cep+geom_point(data=not.ce, aes(x=as.numeric(as.character(Rate1)), y=as.numeric(as.character(Rate2))),size=3, shape=10, inherit.aes = T, show.legend=F)+theme_bw(base_size = 14)+geom_line(data=as.data.frame(points), aes(x=V1, y=V2), inherit.aes = F)
  }else{ 
    
    if(dim(mean.ff)[1]!=7) 
    {
      colnames(mean.ff)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
      cep1<-cep+geom_point(data=mean.ff, aes(x=as.numeric(as.character(Rate1)), y=as.numeric(as.character(Rate2))),size=3, shape=19, stroke=1, inherit.aes = T, show.legend=F)+geom_path(data=mean.ff, aes(x=as.numeric(Rate1), y=as.numeric(Rate2)), inherit.aes = F)+geom_point(data=not.ce, aes(x=as.numeric(as.character(not.ce$Rate1)), y=as.numeric(as.character(not.ce$Rate2))), size=3, shape=10, inherit.aes = T, show.legend=F)+labs(fill='Intervention')+theme_bw(base_size = 11)+geom_line(data=as.data.frame(points), aes(x=V1, y=V2), inherit.aes = F)
    }else{
      colnames(mean.ff)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
      cep1<-cep+geom_point(data=mean.ff, aes(x=as.numeric(as.character(Rate1)), y=as.numeric(as.character(Rate2))),size=3, shape=19, stroke=1, inherit.aes = T, show.legend=F)+theme_bw(base_size = 11)+geom_line(data=as.data.frame(points), aes(x=V1, y=V2), inherit.aes = F)
    }
  }
  
  
  #cep1<-cep+geom_point(data=mean.ff, aes(x=as.numeric(mean.ff$Rate1), y=as.numeric(mean.ff$Rate2)),size=3, shape=19, stroke=1, inherit.aes = T, show.legend=F)+geom_line(data=mean.ff, aes(x=as.numeric(mean.ff$Rate1), y=as.numeric(mean.ff$Rate2)), inherit.aes = F)+geom_point(data=not.ce, aes(x=not.ce$Rate1, not.ce$Rate2), size=3, shape=10, inherit.aes = T, show.legend=F)+labs(fill='Intervention')+theme_bw()
  
  #)
  
  #assign(cep.name, cep1)
  #cep<-cep1<-NULL
  #setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB', 'coverage', cov.f[cov.var]))
  #ggsave(paste0('CEP',strain.name[strain],cov.f[cov.var],dcount,".png"),plot = cep, width=11, height=8.5, units='in',device='png')
  
  not.ce<<-not.ce
  mean.ff<<-mean.ff
  return(cep1)
} 

cep.grapher.short(1,2,3, 20000)

cep.grapher<-function(strain,cov.var, dcount, wtp)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB', 'coverage',cov.f[cov.var]))
  cea.tables<-list.files(pattern=glob2rx(paste0('GBICER.strat', strain.name[strain],cov.f[cov.var],'p', '*', 'd',dcount)))
  i.country<-4
  cos.f<-qds.f<-qds.m<-cos.m<-NULL
  
  for(ff in 1:length(inv.names))
  {#ff=interventions. Strains are already chosen
    load(cea.tables[[ff]]) #H1N1
    
    cer<-ICER.set2 #grab ICERs
    
    #m.qa<-median(qa)
    #m.co<-median(co)
    
    qa1<-Reduce('+', cer$`net QALY`[1:19])/tab.diff
    co1<-Reduce('+', cer$`net cost`[1:19])/tab.diff
    qds<-data.frame(cbind(rep(fname[i.country], length(qa1)), qa1, rep(inv.names[[ff]], length(qa1))))
    cos<-data.frame(cbind(rep(fname[i.country], length(co1)), co1, rep(inv.names[[ff]], length(co1))))
    
    #qa<-unlist(cer$`net QALY`[1:14])
    #co<-unlist(cer$`net cost`[1:14])
    
    #qds<-data.frame(cbind(rep(fname[i.country], length(qa)), qa, rep(inv.names[[ff]], length(qa))))
    #cos<-data.frame(cbind(rep(fname[i.country], length(co)), co, rep(inv.names[[ff]], length(co))))
    
    colnames(qds)<-colnames(cos)
    
    ifelse(ff==1, qds.f<-qds, qds.f<-rbind(qds.f, qds))
    ifelse(ff==1, cos.f<-cos, cos.f<-rbind(cos.f, cos))
    
    meds.f<-cbind(as.numeric(as.character(qds[,2])), as.numeric(as.character(cos[,2])))
    
    xy.points<-data.frame(meds.f)
    names(xy.points)<-c('x', 'y')
    kd <- ks::kde(xy.points, compute.cont=TRUE)
    contour_95 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                        z=estimate, levels=cont["10%"])[[1]])
    contour_95$y[length(contour_95$y)]<-contour_95$y[1]
    contour_95 <- data.frame(contour_95, rep(inv.names[[ff]], length(contour_95$x)))
    colnames(contour_95)<-c('level','x','y','program')
    
    p1<-geom_path(data=contour_95, aes(x, y, color=program))
    contour.name<-paste0('cont',strain.name[strain],cov.var,dcount.l[dcount], ff)
    assign(contour.name,  p1)
  }
  
  
  colnames(qds.f)<-  c('country', 'net QALYS', 'Intervention')
  colnames(cos.f)<- c('country', 'net costs', 'Intervention')
  
  counter<-c()
  cost.utility<-as.matrix(cbind(qds.f[,1:2], cos.f[,2], qds.f[,3]))
  cost.utility<-data.frame(na.omit(unname(cost.utility)))
  colnames(cost.utility)<-c('country', 'QALYS', 'Costs', 'Intervention')
  #row.names(cost.utility) <- c()
  #attributes(cost.utility)$row.names<-NULL
  #mean.p<-ddply(cost.utility, .(country, Intervention), summarize, Rate1=mean(as.numeric(as.character(QALYS))), Rate2=mean(as.numeric(as.character(Costs))), sd.q=sd(as.numeric(as.character(QALYS))),sd.c=sd(as.numeric(as.character(Costs))))
  
  #mean.po<-mean.p[order(as.numeric(as.character(mean.p$Rate2)), decreasing = T),]
  ##   continent n_uniq_countries
  
  
mean.p<-ddply(as.data.frame(cost.utility), .(Intervention), summarize, Rate1=mean(as.numeric(as.character(QALYS))), Rate2=mean(as.numeric(as.character(Costs))), costLB=quantile(as.numeric(as.character(Costs)), 0.025), costUB=quantile(as.numeric(as.character(Costs)), 0.975), qalyLB=quantile(as.numeric(as.character(QALYS)), 0.025), qalyUB=quantile(as.numeric(as.character(QALYS)), 0.975), cost.sd=sd(as.numeric(as.character(Costs))), qal.sd=sd(as.numeric(as.character(QALYS))))
  
  mean.po<-mean.p[order(mean.p$Rate2, decreasing = F),]
  
  not.ce<-not.ce1<-not.ce2<-not.ce3<-not.ce4<-counter<-NULL
  
  rudi.icer<-diff(rbind(c(0,0),cbind(mean.po$Rate1, mean.po$Rate2)))*-1
  rudi.icer.ub<-diff(rbind(c(0,0),cbind(mean.po$costUB, mean.po$qalyUB)))*-1
  rudi.icer.lb<-diff(rbind(c(0,0),cbind(mean.po$costLB, mean.po$qalyLB)))*-1
  rudi.icer.sd<-diff(rbind(c(0,0),cbind(mean.po$cost.sd, mean.po$qal.sd)))*-1
  
  div<-rudi.icer[,2]/rudi.icer[,1]
  div.ub<-rudi.icer.ub[,1]/rudi.icer.ub[,2]
  div.lb<-rudi.icer.lb[,1]/rudi.icer.lb[,2]
  div.sd<-rudi.icer.sd[,1]/rudi.icer.sd[,2]
  
  flip<-as.matrix(mean.po[which(div >0&div <wtp),])
  
  if(dim(flip)[1]<2)
  {
    mean.f<-data.frame(t(as.matrix(c(flip[1:9], 
                                     round(div[which(div >0&div <wtp)]), 
                                     round(div.ub[which(div >0&div <wtp)]), #quantile
                                     round(div.lb[which(div >0&div <wtp)]), #quantile
                                     round(div[which(div >0&div <wtp)]-1.96*div.sd[which(div >0&div <wtp)]), round(div[which(div >0&div <wtp)]+1.96*div.sd[which(div >0&div <wtp)]))))) #sd
  }else{
    mean.f<-data.frame(cbind(as.matrix(flip[,1:9]),
                             round(div[which(div >0&div <wtp)]), 
                             round(div.ub[which(div >0&div <wtp)]), #quantile
                             round(div.lb[which(div >0&div <wtp)]), #quantile
                             round(div[which(div >0&div <wtp)]-1.96*div.sd[which(div >0&div <wtp)]), round(div[which(div >0&div <wtp)]+1.96*div.sd[which(div >0&div <wtp)]))) #sd
  }  
  
  sq<-matrix(c('SQ: Risk Groups, 65+ yrs old',rep(0,13)), nrow=1)
  
  if(dim(mean.f)[2]>0)
  { colnames(mean.f)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
  levels(mean.f$Intervention)= c(levels(mean.f$Intervention),'SQ: Risk Groups, 65+ yrs old')
  mean.ff<-data.frame(as.matrix(mean.f[order(mean.f$Rate2, decreasing=F),]))
  colnames(mean.ff)<-colnames(mean.f)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
  }else{
    mean.ff<-NULL
  }
  
  not.ce<-mean.po[which(div <0 | div >wtp),]
  counter<-1
  
  if(dim(not.ce)[1]!=0)
  {
    counter<-counter+1
    
    ##second icer eval
    rudi.icer2<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff[,2])), as.numeric(as.character(mean.ff[,3])))))*-1
    rudi.icer.ub2<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff$costUB)), as.numeric(as.character(mean.ff$qalyUB)))))*-1
    rudi.icer.lb2<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff$costLB)), as.numeric(as.character(mean.ff$qalyLB)))))*-1
    rudi.icer.sd2<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff$cost.sd)), as.numeric(as.character(mean.ff$qal.sd)))))*-1
    #rudi.icer.sd<-diff(rev(as.numeric(as.character(mean.ff$sd))))
    
    
    div2<-rudi.icer2[,2]/rudi.icer2[,1]
    div.ub2<-rudi.icer.ub2[,1]/rudi.icer.ub2[,2]
    div.lb2<-rudi.icer.lb2[,1]/rudi.icer.lb2[,2]
    div.sd2<-rudi.icer.sd2[,1]/rudi.icer.sd2[,2]
    
    flip2<-as.matrix(mean.f[which(div2 >0&div2 <wtp),])
    
    if(dim(flip2)[1]<2)
    {
      
      mean.f2<-data.frame(t(as.matrix(c(flip2[1:9], 
                                        round(div2[which(div2 >0&div2 <wtp)]), 
                                        round(div.ub2[which(div2 >0&div2 <wtp)]), #quantile
                                        round(div.lb2[which(div2 >0&div2 <wtp)]), #quantile
                                        round(div2[which(div2 >0&div2 <wtp)]-1.96*div.sd2[which(div2 >0&div2 <wtp)]), round(div2[which(div2 >0&div2 <wtp)]+1.96*div.sd2[which(div2 >0&div2 <wtp)]))))) #sd
    }else{
      mean.f2<-data.frame(cbind(as.matrix(flip2[,1:9]), 
                                round(div2[which(div2 >0&div2 <wtp)]), 
                                round(div.ub2[which(div2 >0&div2 <wtp)]), #quantile
                                round(div.lb2[which(div2 >0&div2 <wtp)]), #quantile
                                round(div2[which(div2 >0&div2 <wtp)]-1.96*div.sd2[which(div2 >0&div2 <wtp)]), round(div2[which(div2 >0&div2 <wtp)]+1.96*div.sd2[which(div2 >0&div2 <wtp)]))) #sd
    }
    
    
    
    if(dim(mean.f2)[2]==14)
    {
      colnames(mean.f2)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
      levels(mean.f2$Intervention)= c(levels(mean.f2$Intervention),'SQ: Risk Groups, 65+ yrs old')
      mean.ff2<-data.frame(as.matrix(mean.f2[order(mean.f2$Rate2, decreasing=F),]))
      colnames(mean.ff2)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
    }else{
      mean.ff2<-NULL
    }
    
    not.ce2<-rbind(not.ce, mean.f[which(div2 <0 | div2 >wtp),1:9])
    
    
    
    
    if(dim(not.ce2)[1]!=3)
    {
      if(any(dim(not.ce2)!=dim(not.ce)))
      {
        counter<-counter+1
        ##3rd CE check
        
        rudi.icer3<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff2[,2])), as.numeric(as.character(mean.ff2[,3])))))*-1
        rudi.icer.ub3<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff2$costUB)), as.numeric(as.character(mean.ff2$qalyUB)))))*-1
        rudi.icer.lb3<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff2$costLB)), as.numeric(as.character(mean.ff2$qalyLB)))))*-1
        rudi.icer.sd3<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff2$cost.sd)), as.numeric(as.character(mean.ff2$qal.sd)))))*-1
        #rudi.icer.sd<-diff(rev(as.numeric(as.character(mean.ff$sd))))
        
        
        div3<-rudi.icer3[,2]/rudi.icer3[,1]
        div.ub3<-rudi.icer.ub3[,1]/rudi.icer.ub3[,2]
        div.lb3<-rudi.icer.lb3[,1]/rudi.icer.lb3[,2]
        div.sd3<-rudi.icer.sd3[,1]/rudi.icer.sd3[,2]
        
        flip3<-as.matrix(mean.f2[which(div3 >0&div3 <wtp),])
        
        mean.f3<-data.frame(cbind(as.matrix(flip3[,1:9]), 
                                  round(div3[which(div3 >0&div3 <wtp)]), 
                                  round(div.ub3[which(div3 >0&div3 <wtp)]), #quantile
                                  round(div.lb3[which(div3 >0&div3 <wtp)]), #quantile
                                  round(div3[which(div3 >0&div3 <wtp)]-1.96*div.sd3[which(div3 >0&div3 <wtp)]), round(div3[which(div3 >0&div3 <wtp)]+1.96*div.sd3[which(div3 >0&div3 <wtp)]))) #sd
        
        
        if(dim(mean.f3)[2]==14)
        { colnames(mean.f3)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
        levels(mean.f3$Intervention)= c(levels(mean.f3$Intervention),'SQ: Risk Groups, 65+ yrs old')
        mean.ff3<-data.frame(rbind(sq, as.matrix(mean.f3[order(mean.f2$Rate2, decreasing=F),])))
        colnames(mean.ff3)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
        }else
        {mean.ff3<-NULL}
        not.ce3<-rbind(not.ce2, mean.f2[which(div3 <0 | div3 >wtp),1:9])
    

    
    if(dim(not.ce3)[1]!=4)
    {
      if(any(dim(not.ce3)!=dim(not.ce2)))
      {
        counter<-counter+1
        ##3rd CE check
        
        rudi.icer4<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff3[,2])), as.numeric(as.character(mean.ff3[,3])))))*-1
        rudi.icer.ub4<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff3$costUB)), as.numeric(as.character(mean.ff3$qalyUB)))))*-1
        rudi.icer.lb4<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff3$costLB)), as.numeric(as.character(mean.ff3$qalyLB)))))*-1
        rudi.icer.sd4<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff3$cost.sd)), as.numeric(as.character(mean.ff3$qal.sd)))))*-1
        #rudi.icer.sd<-diff(rev(as.numeric(as.character(mean.ff$sd))))
        
        
        div4<-rudi.icer4[,2]/rudi.icer4[,1]
        div.ub4<-rudi.icer.ub4[,1]/rudi.icer.ub4[,2]
        div.lb4<-rudi.icer.lb4[,1]/rudi.icer.lb4[,2]
        div.sd4<-rudi.icer.sd4[,1]/rudi.icer.sd4[,2]
        
        flip4<-as.matrix(mean.f3[which(div4 >0&div4 <wtp),])
        
        #if(dimflip3[1]!=0)
        #{
        mean.f4<-data.frame(cbind(as.matrix(c(flip4[1:9], 
                                          round(div4[which(div4 >0&div4 <wtp)]), 
                                          round(div.ub4[which(div4 >0&div4 <wtp)]), #quantile
                                          round(div.lb4[which(div4 >0&div4 <wtp)]), #quantile
                                          round(div4[which(div4 >0&div4 <wtp)]-1.96*div.sd4[which(div4 >0&div4 <wtp)]), round(div4[which(div4 >0&div4 <wtp)]+1.96*div.sd4[which(div4 >0&div4 <wtp)]))))) #sd
        
        
        #}
        
        if(dim(mean.f4)[2]==14)
        { colnames(mean.f4)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
        levels(mean.f4$Intervention)= c(levels(mean.f4$Intervention),'SQ: Risk Groups, 65+ yrs old')
        mean.ff4<-data.frame(rbind(sq, as.matrix(mean.f4[order(mean.f3$Rate2, decreasing=F),])))
        colnames(mean.ff4)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
        }else
        {mean.ff4<-NULL}
        
        not.ce4<-rbind(not.ce3, mean.f3[which(div4 <0 | div4 >wtp),1:9])
        
    
    if(dim(not.ce4)[1]!=5)
    {
      if(any(dim(not.ce5)!=dim(not.ce4)))
      {
        counter<-counter+1
        ##3rd CE check
        
        rudi.icer5<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff4[,2])), as.numeric(as.character(mean.ff4[,3])))))*-1
        rudi.icer.ub5<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff4$costUB)), as.numeric(as.character(mean.ff4$qalyUB)))))*-1
        rudi.icer.lb5<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff4$costLB)), as.numeric(as.character(mean.ff4$qalyLB)))))*-1
        rudi.icer.sd5<-diff(rbind(c(0,0),cbind(as.numeric(as.character(mean.ff4$cost.sd)), as.numeric(as.character(mean.ff4$qal.sd)))))*-1
        #rudi.icer.sd<-diff(rev(as.numeric(as.character(mean.ff$sd))))
        
        
        div5<-rudi.icer5[,2]/rudi.icer5[,1]
        div.ub5<-rudi.icer.ub5[,1]/rudi.icer.ub5[,2]
        div.lb5<-rudi.icer.lb5[,1]/rudi.icer.lb5[,2]
        div.sd5<-rudi.icer.sd5[,1]/rudi.icer.sd5[,2]
        
        flip5<-as.matrix(mean.f4[which(div5 >0&div5 <wtp),])
        
        #if(dimflip3[1]!=0)
        #{
        mean.f5<-data.frame(cbind(as.matrix(c(flip5[1:9], 
                                          round(div5[which(div5 >0&div5 <wtp)]), 
                                          round(div.ub5[which(div5 >0&div5 <wtp)]), #quantile
                                          round(div.lb5[which(div5 >0&div5 <wtp)]), #quantile
                                          round(div5[which(div5 >0&div5 <wtp)]-1.96*div.sd5[which(div5 >0&div5 <wtp)]), round(div5[which(div5 >0&div5 <wtp)]+1.96*div.sd5[which(div5 >0&div5 <wtp)]))))) #sd
        
        
        #}
        
        if(dim(mean.f5)[2]==14)
        { colnames(mean.f5)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
        levels(mean.f5$Intervention)= c(levels(mean.f5$Intervention),'SQ: Risk Groups, 65+ yrs old')
        mean.ff5<-data.frame(rbind(sq, as.matrix(mean.f5[order(mean.f4$Rate2, decreasing=F),])))
        colnames(mean.ff5)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
        }else
        {mean.ff5<-NULL}
        not.ce5<-rbind(not.ce4, mean.f4[which(div5 <0 | div5 >wtp),1:9])
    
    if(dim(not.ce5)[1]!=6)
    {
      if(any(dim(not.ce5)!=dim(not.ce4)))
      {
        counter<-counter+1
        ##3rd CE check
        
        rudi.icer6<-rbind(c(0,0),cbind(as.numeric(as.character(mean.ff5[,2])), as.numeric(as.character(mean.ff5[,3]))))*-1
        rudi.icer.ub6<-rbind(c(0,0),cbind(as.numeric(as.character(mean.ff5$costUB)), as.numeric(as.character(mean.ff5$qalyUB))))*-1
        rudi.icer.lb6<-rbind(c(0,0),cbind(as.numeric(as.character(mean.ff5$costLB)), as.numeric(as.character(mean.ff5$qalyLB))))*-1
        rudi.icer.sd6<-rbind(c(0,0),cbind(as.numeric(as.character(mean.ff5$cost.sd)), as.numeric(as.character(mean.ff5$qal.sd))))*-1
        #rudi.icer.sd<-diff(rev(as.numeric(as.character(mean.ff$sd))))
        
        
        div6<-rudi.icer6[,2]/rudi.icer6[,1]
        div.ub6<-rudi.icer.ub6[,1]/rudi.icer.ub6[,2]
        div.lb6<-rudi.icer.lb6[,1]/rudi.icer.lb6[,2]
        div.sd6<-rudi.icer.sd6[,1]/rudi.icer.sd6[,2]
        
        flip6<-as.matrix(mean.f5[which(div6 >0&div6 <wtp),])
        
        #if(dimflip3[1]!=0)
        #{
        mean.f6<-data.frame(t(as.matrix(c(flip6[1:9], 
                                          round(div6[which(div6 >0&div6 <wtp)]), 
                                          round(div.ub6[which(div6 >0&div6 <wtp)]), #quantile
                                          round(div.lb6[which(div6 >0&div6 <wtp)]), #quantile
                                          round(div6[which(div6 >0&div6 <wtp)]-1.96*div.sd6[which(div6 >0&div6 <wtp)]), round(div6[which(div6 >0&div6 <wtp)]+1.96*div.sd6[which(div6 >0&div6 <wtp)]))))) #sd
        
        
        #}
        
        if(dim(mean.f6)[2]==14)
        { colnames(mean.f6)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
        levels(mean.f6$Intervention)= c(levels(mean.f6$Intervention),'SQ: Risk Groups, 65+ yrs old')
        mean.ff6<-data.frame(rbind(sq, as.matrix(mean.f6[order(mean.f5$Rate2, decreasing=F),])))
        colnames(mean.ff6)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
        }else
        {mean.ff6<-NULL}
        not.ce6<-rbind(not.ce5, mean.f5[which(div6 <0 | div6 >wtp),1:9])}} 
            }
          }
        }
      }
    }
  }
}
  
  ######CE evaluation
  
  if(counter>1)
  {
    mean.ff<-get(paste0('mean.ff', counter))
    not.ce<-get(paste0('not.ce', counter))
  }
  
  
  ####axis.text.x = element_text(face="bold", color="#993333", 
  #size=14
  cep<-cep1<-NULL
  cols<-c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC")
  
  cep<-ggplot(cost.utility, aes(as.numeric(as.character(QALYS)), as.numeric(as.character(Costs)), color=Intervention))+theme(axis.text.x = element_text(angle = 45, hjust = 1, size=11), plot.subtitle =element_text(angle = 45, hjust = 1, size=10) )+scale_y_continuous(breaks=scales::pretty_breaks(n = 7), limits = c(0,1.5e8))+scale_x_continuous(breaks=scales::pretty_breaks(n=4), limits=c(-2000,61000))+geom_hline(aes(yintercept=0))+geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+geom_abline(slope = 15000, intercept=0, linetype='dotdash')+geom_abline(slope = wtp, intercept=0, linetype='longdash')+labs(subtitle=paste(paste0('discount=',dcount.l[dcount],'%'), paste0(strain.name[strain]), paste0('coverage=',cov.f[cov.var]) ), fill='Intervention')+theme_bw()+scale_colour_manual(values=cols)
  
  
  for(fff in 1:length(inv.names))
  {add.on<-get(paste0('cont',strain.name[strain],cov.var,dcount.l[dcount], fff))
  cep<-cep+add.on}
  #cep.name<-paste0('cep',strain, cov.var, dcount.l[dcount])
  
  points<-cbind(c(0,as.numeric(as.character(mean.ff[,2]))), c(0,as.numeric(as.character(mean.ff[,3]))))
  names(points)<-c('Rate1', 'Rate2')
  
  if(is.null(mean.ff))
  {
    levels(not.ce$Intervention)=c(levels(not.ce$Intervention),'SQ: Risk Groups, 65+ yrs old')
    not.ce<-data.frame(rbind(sq[1:9], as.matrix(not.ce[order(not.ce$Rate2, decreasing=F),])))
    
    cep1<-cep+geom_point(data=not.ce, aes(x=as.numeric(as.character(Rate1)), y=as.numeric(as.character(Rate2))),size=3, shape=10, inherit.aes = T, show.legend=F)+theme_bw(base_size = 14)+geom_line(data=as.data.frame(points), aes(x=V1, y=V2), inherit.aes = F)
  }else{ 
    
    if(dim(mean.ff)[1]!=7) 
    {
      colnames(mean.ff)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
      cep1<-cep+geom_point(data=mean.ff, aes(x=as.numeric(as.character(Rate1)), y=as.numeric(as.character(Rate2))),size=3, shape=19, stroke=1, inherit.aes = T, show.legend=F)+geom_path(data=mean.ff, aes(x=as.numeric(Rate1), y=as.numeric(Rate2)), inherit.aes = F)+geom_point(data=not.ce, aes(x=as.numeric(as.character(not.ce$Rate1)), y=as.numeric(as.character(not.ce$Rate2))), size=3, shape=10, inherit.aes = T, show.legend=F)+labs(fill='Intervention')+theme_bw(base_size = 11)+geom_line(data=as.data.frame(points), aes(x=V1, y=V2), inherit.aes = F)
    }else{
      colnames(mean.ff)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
      cep1<-cep+geom_point(data=mean.ff, aes(x=as.numeric(as.character(Rate1)), y=as.numeric(as.character(Rate2))),size=3, shape=19, stroke=1, inherit.aes = T, show.legend=F)+theme_bw(base_size = 11)+geom_line(data=as.data.frame(points), aes(x=V1, y=V2), inherit.aes = F)
    }
  }
  
  
  #cep1<-cep+geom_point(data=mean.ff, aes(x=as.numeric(mean.ff$Rate1), y=as.numeric(mean.ff$Rate2)),size=3, shape=19, stroke=1, inherit.aes = T, show.legend=F)+geom_line(data=mean.ff, aes(x=as.numeric(mean.ff$Rate1), y=as.numeric(mean.ff$Rate2)), inherit.aes = F)+geom_point(data=not.ce, aes(x=not.ce$Rate1, not.ce$Rate2), size=3, shape=10, inherit.aes = T, show.legend=F)+labs(fill='Intervention')+theme_bw()
  
  #)
  
  #assign(cep.name, cep1)
  #cep<-cep1<-NULL
  #setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB', 'coverage', cov.f[cov.var]))
  #ggsave(paste0('CEP',strain.name[strain],cov.f[cov.var],dcount,".png"),plot = cep, width=11, height=8.5, units='in',device='png')
  
  not.ce<<-not.ce
  mean.ff<<-mean.ff
  return(cep1)
} 
      

#test graph
cep.grapher(1,3,3, 15000)
#debug(cep.graph)
#setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB', 'coverage', '0.55'))
library(ggpubr)
wtpl<-20000
all323<-ggarrange(cep.grapher(1,2,3, wtpl)+rremove("xy.title"), cep.grapher(2,2,3, wtpl)+rremove("xy.title")+rremove("xy.title"), cep.grapher(3,2,3, wtpl)+rremove("xy.title")+rremove("xy.title"), common.legend=T, legend='bottom',nrow=1, ncol=3)

all323<-annotate_figure(all323,
                top = text_grob("Cost-Effectiveness Plane"),
                bottom = text_grob("QALY Differences"),
                left = text_grob("Cost in million (₤GBP) /year", rot=90)
)

ggsave(paste0('CEPstrainspecdiff',cov.f[cov.var],version,'3.5',".png"),plot = all323, width=13, height=5.75, units='in',device='png')


###########################
wtpl<-20000
all4<-ggarrange(cep.grapher.all(2,3,wtpl)+rremove("xy.title")+rremove('x.text')+rremove('x.ticks'), 
                cep.grapher(1,2,3, wtpl)+rremove("xy.title")+rremove('y.text')+rremove('y.ticks')+rremove('x.text')+ rremove('x.ticks'), 
                cep.grapher(2,2,3, wtpl)+rremove("xy.title"), 
                cep.grapher(3,2,3, wtpl)+rremove("xy.title")+rremove('y.text')+rremove('y.ticks'),
                common.legend=T, legend='bottom',nrow=2, ncol=2)

all4.2<-annotate_figure(all4,
                        top = text_grob("Cost-Effectiveness Plane"),
                        bottom = text_grob("QALY Differences"),
                        left = text_grob("Cost in million (₤GBP) /year", rot=90)
)

version<-'v1'
ggsave(paste0('CEPstrainx4',cov.f[cov.var],version,'3.5',".png"),plot = all4, width=10, height=8.6, units='in',device='png')










plot(x=dates, y=cov7, add=T)

predict(ff.t, dates2)
plot(ff)

(dates2*7.191e-03)+-1.120e+02

##discount variation
discount3<-
  lapply(1:3, function(dcount) ggarrange(
    cep.grapher(1,1,dcount, wtpl)+rremove("xy.title"), 
    cep.grapher(2,1,dcount, wtpl)+rremove("xy.title"), 
    cep.grapher(3,1,dcount,  wtpl)+rremove("xy.title"), 
    cep.grapher(1,2,dcount,  wtpl)+rremove("xy.title"), 
    cep.grapher(2,2,dcount, wtpl)+rremove("xy.title"), 
    cep.grapher(3,2,dcount, wtpl)+rremove("xy.title"),
    cep.grapher(1,3,dcount,  wtpl)+rremove("xy.title"), 
    cep.grapher(2,3,dcount,  wtpl)+rremove("xy.title"), 
    cep.grapher(3,3,dcount,  wtpl)+rremove("xy.title"),
    ncol=3, nrow=3, common.legend = T, legend='bottom'))

library('gridExtra')
graphs.allpost<-marrangeGrob(discount3,
                             nrow = 1, ncol=1,
                top = text_grob("Cost-Effectiveness Plane"),
                bottom = text_grob("QALY Differences"),
                left = text_grob("Incremental Cost in million ₤/year", rot=90)
)

ggsave(paste0("CEPspecstrainvcoverage",version, ".pdf"),plot=graphs.allpost, width=10, height=8, units='in',device='pdf')


##################
pdf(file='cepintvstratbyage1.pdf',onefile = T, height=8.5, width=11)
for(j in 1:17)
{
  print(ggplot(llg, aes(as.numeric(as.character(qaly)), cost, color=program))+ theme_grey()+
          geom_point(alpha=.5)+stat_ellipse(type='norm', alpha = .9, aes(color = program), level=0.90)+scale_y_continuous(breaks=scales::pretty_breaks(n = 8))+scale_x_continuous(breaks=scales::pretty_breaks(n=8))+geom_hline(aes(yintercept=0))+geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+geom_abline(slope = 15000, intercept=0, color='blue')+geom_abline(slope = 20000, intercept=0, color='light blue')+facet_wrap_paginate(Age~program, ncol = 1, nrow = 3, page=j, scales='free')+xlab('Contact Matrix')+ylab('Cost per 1 QALY Gained')+ labs(fill='Intervention'))
  
}
dev.off()


#######################################################################################
cols<-c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC")

inv.names<-c("I1: 2-4 yrs old"  ,                 "I2: 5-11 yrs old"       ,              "I3: 12-16 yrs old"      ,             "I4: 2-11 yrs old" ,"I7: 2-16 yrs old", "I5: 2-4 & 12-16 yrs old"   ,      "I6: 5-16 yrs old" )
#inv.nam

discount.loader<-function(strain, cov.var, dcount)
{
  #setwd(paste0(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB','coverage', as.character(cov.f[cov.var]))))
  
  #Rearrange tables to correct order of interventions------------------
  cost.tables<-list.files(pattern=glob2rx(paste0('GBCEAcosts', strain.name[strain],cov.var, 'p*',  dcount)))
  
  load(cost.tables[1]); avg.incidence1<-coststab; #status quo
  load(cost.tables[2]); avg.incidence2<-coststab; #preschool
  load(cost.tables[3]); avg.incidence3<-coststab; #primary
  load(cost.tables[4]); avg.incidence4<-coststab; #secondary
  load(cost.tables[5]); avg.incidence5<-coststab; #secondary
  load(cost.tables[6]); avg.incidence6<-coststab; #secondary
  load(cost.tables[7]); avg.incidence7<-coststab; #secondary
  
  
  var.name<-paste0('inv.costs',strain.name[strain], dcount.l[dcount])
  #('WTP','Preschool','Primary','Secondary','Preschool+Primary','PrePrimeSeconday','Preschool+Secondary','Primary+Secondary')
  assign(var.name, do.call( rbind.data.frame, list(avg.incidence1,avg.incidence2,avg.incidence3, avg.incidence4, avg.incidence5,avg.incidence6, avg.incidence7)
  ), envir = .GlobalEnv)
  
}

###############################################
#########Graph barplot for economic outcomes by strain and discount
#coverage not specified here so be careful

barplot<-function(strain, cov.var, dcount)
{
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB','coverage', as.character(cov.f[cov.var])))
  
  #make status quo category
  for(tt in 1:2)
  {
    discount.loader(strain,cov.f[cov.var], dcount)
    inv.costs<-get(paste0('inv.costs',strain.name[strain], dcount.l[dcount]))
    
    inv.costs$variable <- factor(inv.costs$variable, levels=c(paste0('vaccine cost',tt), paste0('Healthcare costs',tt), paste0('GP costs',tt), paste0('hospital costs',tt),paste0('QALYs lost',tt),paste0('QALYs lost case',tt), paste0('QALYs lost death',tt), paste0('QALYs lost hosp', tt)), labels = c('Vaccination Costs','Healthcare Costs','GP Consult Fees','Hospitalization Costs','Years of Life Lost (All Sources)', 'Years of Life Lost (Symptomatic Cases)', 'Years of Life Lost (Deaths)','Years of Life Lost (Hospitalization)'))
    
    inv.costs$strategy<-factor(inv.costs$strategy, levels<-levels(inv.costs$strategy), labels=inv.names)
    
    if(tt==1) inv.costs$strategy<-"SQ: Risk Groups, 65+ yrs old"
    ifelse(tt==1, inv.costs.f<-na.omit(inv.costs)[1:6,], inv.costs.f<-rbind(inv.costs.f, na.omit(inv.costs)))
  }
  #cols<-my.col
  
  q<-ggplot(inv.costs.f, aes(x=strategy, y=value, fill=strategy))+geom_bar(position = "dodge", stat = "identity", aes(fill=strategy))+facet_wrap(~variable, ncol=3, scales = 'free')+scale_y_continuous(breaks=scales::pretty_breaks(n=7))+labs(title = NULL,subtitle=paste('discount=', dcount.l[dcount], '%', ',','coverage=', cov.f[cov.var],',', 'strain=', strain.name[strain]), x = 'Intervention Strategy', y = 'Average Annual Costs(£)', fill = "Economic Impact\n") +geom_errorbar(aes(ymin=CIL, ymax=CIU), width=.2, position=position_dodge(.9))+theme_bw()+theme(
    axis.text.x = element_blank(), axis.ticks = element_blank(), plot.subtitle=element_text(size=10, hjust=0.5, color="black"), legend.position="bottom")+scale_fill_manual(values=cols)
  
  
  ggsave(paste0('costsplot',strain.name[strain],cov.f[cov.var],dcount.l[dcount],".png"),plot = q, width=11, height=7.5, units='in',device='png')
  return(q)
}



library(ggpubr)
library(gridExtra)

for(strain in 1:3)
{
##arrange by discount and coverage level
barplot1<-ggarrange(
  barplot(strain,cov.var=1,dcount=1)+rremove("xy.title"), 
  barplot(strain,2,1)+rremove("xy.title"),
  barplot(strain,3,1)+rremove("xy.title"),
  barplot(strain,1,2)+rremove("xy.title"),
  barplot(strain,2,2)+rremove("xy.title"),
  barplot(strain,3,2)+rremove("xy.title"),
  barplot(strain,1,3)+rremove("xy.title"),
  barplot(strain,2,3)+rremove("xy.title"),
  barplot(strain,3,3)+rremove("xy.title"),
  ncol=1, nrow=1, common.legend = F, legend='bottom')

barplot32<-marrangeGrob(barplot1, nrow=1, ncol=1,
                        top = text_grob(paste("Economic Outcomes Strain", strain.name[strain])),
                        bottom = text_grob("Intervention Strategy"),
                        left = text_grob("Average Annual Costs(£)", rot=90)
)

ggsave(paste0("barplot",strain.name[strain],'u7',".pdf"),plot=barplot32, width=11, height=8.5, units='in',device='pdf')
}
dev.off()


###
  

barplot.compare<-function(cov.var, dcount)
{
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB', 'coverage', cov.f[cov.var]))
  
  discount.loader(1,cov.var, dcount)
  discount.loader(2,cov.var, dcount)
  discount.loader(3,cov.var, dcount)
  
  #inv.costs<-get(paste0('inv.costs',strain.name[strain], dcount.l[dcount]))
  
  #make status quo category
  for(strain in 1:3)
  {
  for(tt in 1:2)
  {
    inv.costs<-get(paste0('inv.costs',strain.name[strain], dcount.l[dcount]))
    
    inv.costs$variable <- factor(inv.costs$variable, levels=c(paste0('vaccine cost',tt), paste0('Healthcare costs',tt), paste0('GP costs',tt), paste0('hospital costs',tt),paste0('QALYs lost',tt),paste0('QALYs lost case',tt), paste0('QALYs lost death',tt), paste0('QALYs lost hosp', tt)), labels = c('Vaccination Costs','Healthcare Costs','GP Consult Fees','Hospitalization Costs','Years of Life Lost (All Sources)', 'Years of Life Lost (Symptomatic Cases)', 'Years of Life Lost (Deaths)','Years of Life Lost (Hospitalization)'))
    
    inv.costs$strategy<-factor(inv.costs$strategy, levels<-levels(inv.costs$strategy), labels=inv.names)
    
    if(tt==1) inv.costs$strategy<-"SQ: Risk Groups, 65+ yrs old"
    ifelse(tt==1, inv.costs.f<-na.omit(inv.costs)[1:6,], inv.costs.f<-rbind(inv.costs.f, na.omit(inv.costs)))
  }
    ifelse(strain==1, inv.costs.g<-inv.costs.f, inv.costs.g<-rbind(inv.costs.g, inv.costs.f))
  }
  #cols<-my.col
  
  ind.scales<-c(c(-5000,NA),c(-5000,NA),c(-5000,NA),c(-5000,NA), c(-3000,NA), c(-100,NA))
  
  q<-ggplot(inv.costs.g, aes(x=strategy, y=value, fill=strategy, group=strain))+geom_bar(position = "dodge", stat = "identity",color='black', aes(fill=strategy, alpha=strain))+facet_wrap(~variable, ncol=3, scales = 'free')+scale_y_continuous(breaks=scales::pretty_breaks(n=7))+labs(title = NULL,subtitle=paste('discount=', dcount.l[dcount], '%', ',','coverage=', cov.f[cov.var]), x = 'Intervention Strategy', y = 'Average Annual Costs(£)', fill = "Economic Impact\n") +geom_errorbar(aes(ymin=CIL, ymax=CIU), width=.2, position=position_dodge(.9))+theme_bw()+theme(
    axis.text.x = element_blank(), axis.ticks = element_blank(), plot.subtitle=element_text(size=10, hjust=0.5, color="black"), legend.position="bottom")+scale_fill_manual(values=cols)+scale_alpha_manual(values=c(0.95,0.65,0.35),
                     name="Strain",
                     breaks=c("H1N1", "H3N2", "B"),         # can use to set order
                     labels=c("H1N1", "H3N2", "B"))
  
  #+geom_text(aes(x=strategy, y=-20, label=strain),position = position_dodge(width=0.9), color = "black",angle = 90, size=2.5)
  
  
  ggsave(paste0('costsplotvsstrain',cov.f[cov.var],dcount.l[dcount],".pdf"),plot = q, width=11, height=7.5, units='in',device='pdf')
  return(q)
}

for(dcount in 1:3) {for(cov.var in 1:3) {barplot.compare(cov.var, dcount)}}


####################################################################################
####Acceptability curve, probability each intervention has optimal cost-effectiveness
###########################################################################################

inv.names2<-c('Status Quo',"Preschool"   ,                "Primary"      ,               "Secondary"  ,                 "Preschool+Primary School"   , "Preschool+Primary+Secondary" ,"Preschool+Secondary"   ,      "Primary+Secondary")

setwd(paste0('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/GB/coverage/', cov.f[cov.var]))

accept.curve<-function(threshold.val, qaly, cost)
{
  qaly2<- Map('*', -1, Map('*',qaly,threshold.val))
  nmb<-Map('-',qaly2, cost)
  
  return(nmb)
}

c.names<-c('program','0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR', 'total')


  
optimal.strat.graph<-function(strain, cov.var, dcount)
{
  qal.all<-list()
  cost.all<-list()
  outerf<-list()
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage',cov.f[cov.var]))
  qal.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.qaly.loss', strain.name[strain],cov.f[cov.var], 'p*', 'd',dcount.l[dcount])))
  cos.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.cost.loss', strain.name[strain],cov.f[cov.var], '*', 'd',dcount.l[dcount])))
  
 for(ii in 1:length(qal.tables)) 
 {
   #optimal strategy compares all strategy outcomes to each other and sums to 1
  load(qal.tables[[ii]])
  load(cos.tables[[ii]])
  
  if(ii==1){qal.all[[ii]]<-tot.QALY.loss[[1]];
            cost.all[[ii]]<-tot.costs.loss[[1]];
  }
  qal.all[[ii+1]]<-tot.QALY.loss[[2]]
  cost.all[[ii+1]]<-tot.costs.loss[[2]]
}
  
  
  ss<-seq(0,30000,by = 1000)
  for(hh in 1:length(ss)) {
    
    AC.all<-lapply(1:length(qal.all), function(oo) accept.curve(ss[hh], qal.all[[oo]], cost.all[[oo]]) )
    
   AC.all2<-lapply(1:length(qal.all), function(pp) data.frame(as.factor(rep(inv.names2[pp],n.samples)), Reduce('+',AC.all[[pp]])/tab.diff, rowSums(Reduce('+',AC.all[[pp]])/tab.diff)))
   
   lapply(AC.all2, setNames, c.names)
   
   new.ac<-array(unlist(AC.all2), dim=c(n.samples, dim(AC.all2[[1]])[2], 8))
   
    for(aa in 2:24){
      e.frame<-data.frame(new.ac[,aa,1:8])
      xx<-apply(e.frame, 1, function(z) {which.max(z)})  
      ex.frame<-matrix(rep(0,dim(e.frame)[1]*dim(e.frame)[2]), ncol=8)
      for(cc in 1:length(xx))
      {ex.frame[cc,xx[cc]]<-1}
      
      
      outp<-data.frame('GB', c.names[aa], ss[hh], t(colMeans(ex.frame)))
      
      ifelse(aa==2&hh==1, outerf<-outp, outerf<-rbind(outerf, outp))
    }
  }
  #ifelse(ii==1, outerf.all<-outerf, outerf.all<-cbind(outerf.all, outerf$X2))


#setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')
colnames(outerf)<-c('country', 'Age','threshold','Status Quo',"Preschool"   ,                "Primary"      ,               "Secondary"  ,                 "Preschool+Primary School"   , "Preschool+Primary+Secondary" ,"Preschool+Secondary"   ,      "Primary+Secondary")


#dev.off()
setwd(paste0('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/GB/coverage/', cov.f[cov.var]))

#for totals only
outer.g<-melt(outerf, id.vars=c('country','Age', 'threshold'))
colnames(outer.g)<-c('country', 'Age', 'Threshold','Intervention','value')
mk<-ggplot(data= outer.g[outer.g$Age=='total',], aes(x=as.numeric(as.character(Threshold)), y=as.numeric(as.character(value)), color=Intervention, linetype=Intervention))

label<-paste(paste0('strain=', strain.name[strain]),paste0('discount=', dcount.l[dcount], '%'), paste0('coverage=', cov.f[cov.var]))

OC3<-mk+labs(title=NULL, subtitle=label, fill='Intervention')+geom_line(aes(color=Intervention))+scale_y_continuous(breaks=scales::pretty_breaks(n=7))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_vline(xintercept = c(15000,20000), linetype=3)+theme(plot.subtitle=element_text(size=10, hjust=0.5, color="black"))+guides(fill=guide_legend(title="Intervention"))

  if(cov.var==2&dcount==3)
    {ggsave(paste0('OC.total',strain.name[strain],cov.f[cov.var],'d',dcount.l[dcount],".png"),plot = OC3, width=11,height=8.5, units='in',
         device='png')};

return(OC3)
}


#arrange Optimal graph output
optimal3<-NULL
optimal3<-lapply(1:3, function(dcount) ggarrange(
    optimal.strat.graph(1,1,dcount)+rremove("xy.title"), optimal.strat.graph(2,1,dcount)+rremove("xy.title"), optimal.strat.graph(3,1,dcount)+rremove("xy.title"), 
    optimal.strat.graph(1,2,dcount)+rremove("xy.title"), optimal.strat.graph(2,2,dcount)+rremove("xy.title"), optimal.strat.graph(3,2,dcount)+rremove("xy.title"),
    optimal.strat.graph(1,3,dcount)+rremove("xy.title"), optimal.strat.graph(2,3,dcount)+rremove("xy.title"), optimal.strat.graph(3,3,dcount)+rremove("xy.title"),
    ncol=3, nrow=3, common.legend = T, legend='bottom'))

graphs.ocpost<-NULL
graphs.ocpost<-marrangeGrob(optimal3, nrow = 1, ncol=1, 
                            top = text_grob('Optimal Acceptability Curve per Strategy'),
                bottom = text_grob('Willingness-to-pay Threshold (£)'),
                left = text_grob('Probability Cost-Effective', rot=90))

ggsave(paste0("OC3.strainvcoverage",".png"),plot=graphs.ocpost, width=12, height=9, units='in',device='png')


#############################################################
#Attack Rate
###############################################################

setwd( '/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/GB')


AR.loader<-function(strain, cov.var, dcount)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB'))
  #Rearrange tables to correct order of interventions------------------
  AR.tables<-list.files(pattern=glob2rx(paste0('GBAttackRate', strain.name[strain],cov.f[cov.var], 'p*', 'd', dcount.l[dcount])))
  
  load(AR.tables[1]); avg.incidence1<-AR; #status quo
  load(AR.tables[2]); avg.incidence2<-AR; #preschool
  load(AR.tables[3]); avg.incidence3<-AR; #primary
  load(AR.tables[4]); avg.incidence4<-AR; #secondary
  load(AR.tables[5]); avg.incidence5<-AR; #secondary
  load(AR.tables[6]); avg.incidence6<-AR; #secondary
  load(AR.tables[7]); avg.incidence7<-AR; #secondary
  
  
  var.name<-paste0('ar',strain.name[strain],cov.f[cov.var], dcount.l[dcount])
  #('WTP','Preschool','Primary','Secondary','Preschool+Primary','PrePrimeSeconday','Preschool+Secondary','Primary+Secondary')
  assign(var.name, list(avg.incidence1[[1]],avg.incidence1[[2]],avg.incidence2[[2]],avg.incidence3[[2]], avg.incidence4[[2]], avg.incidence5[[2]],avg.incidence6[[2]], avg.incidence7[[2]])
  , envir = .GlobalEnv)
  
}

library(viridis)
for(strain in 1:3)
{
for(dcount in 1:3)
{
for(cov.var in 1:3)
{

AR.loader(strain,cov.var,dcount)

  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage',cov.f[cov.var]))
  
ar.vals<-data.frame(cbind(matrix(unlist(lapply(get(paste0('ar', strain.name[strain], cov.f[cov.var], dcount.l[dcount])), colMedians)), ncol=18, byrow = T), ppname))
colnames(ar.vals)<-c('0-1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 y_HR', '2-4 y_HR', '5-11y_HR', '12-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR', 'program')
ar.vals.g<-melt(ar.vals, id.vars='program')

arv<-ggplot(ar.vals.g,aes(x = as.factor(variable),y = as.factor(program), fill = as.numeric(value)*100)) + 
  geom_tile() + theme(legend.position='bottom', axis.text.x = element_text(angle = 90, hjust = 1), plot.title=element_text(size=10, hjust=0.5, color="black"))+
  scale_fill_viridis(name='Attack Rate Intensity (%)', option='plasma', limits=c(0,50))+labs(x='Stratified Age Groups', y='Intervention', title=paste('Attack Rate Stratfied by Age Group', 'strain=', strain.name[strain], 'coverage=', cov.f[cov.var], ',', 'discount=', dcount.l[dcount]))+ guides(fill = guide_colourbar(nbin = 100))

ggsave(paste0("AR",strain.name[strain],cov.f[cov.var],dcount.l[dcount],".pdf"),plot=arv, width=8, height=5, units='in',device='pdf')
#return(arv)
}
}
}

AR.loader(1,2,3)

########################################################
#for specific age groups
#########################################################
outer.g<-melt(outerf.all, id.vars=c('country','Age', 'threshold'))
m<-ggplot(data=outer.g, aes(x=as.numeric(as.character(threshold)), y=as.numeric(as.character(value)), color=variable, linetype=variable))  
pdf(file='OC.3program.pdf',onefile = T, height=11, width=8.75)
for(gg in 1:32){
  print(m+ facet_wrap_paginate(country~Age, ncol = 2, nrow = 4, scales='free', page=gg)+
          geom_line(aes(color=variable, group=variable, linetype=variable))+labs(title='Probability a strategy is Optimal', x='Willing-to-pay Threshold (£)', y='Probability Cost-Effective',color='Strategy') +
          scale_y_continuous(breaks=scales::pretty_breaks(n=8))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_vline(xintercept = c(15000,20000), linetype=3))
}

dev.off()


##########################################################################################################
# ACCEPTIBILITY CURVE for single intervention probability it is cost effective on WTP scale
#########################################################################################################################


c.names<-c('program','0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR', 'total')
qal.diff<-cost.diff<-list()

prob.curve<-function(threshold.val, strain, cov.var, dcount)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage', cov.f[cov.var]))
  cea.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'ICER.strat', strain.name[strain], cov.f[cov.var],'p', '*', 'd',dcount.l[dcount])))
  
  for(ii in 1:length(cea.tables)) 
  {
    #optimal strategy compares all strategy outcomes to each other and sums to 1
    load(cea.tables[[ii]]);
    
    qal.diff[[ii]]<-ICER.set2$`QALY diff`;
    cost.diff[[ii]]<-ICER.set2$`cost diff`;
  }
  
  qal2<-lapply(1:7, function(x) Map('*', qal.diff[[x]], as.numeric(threshold.val)))
  nmb<-lapply(1:7, function(x) Map('-',qal2[[x]], cost.diff[[x]]))
  nmb2<-lapply(1:7, function(x) as.matrix(Reduce('+',nmb[[x]])/tab.diff))
  nmb3<-lapply(1:7, function(x) cbind(nmb2[[x]], rowSums(nmb2[[x]])))
  
  ex.zeros<-list(
    matrix(rep(0,length(nmb3[[1]])), ncol=dim(nmb3[[1]])[2]),
    matrix(rep(0,length(nmb3[[1]])), ncol=dim(nmb3[[1]])[2]),
    matrix(rep(0,length(nmb3[[1]])), ncol=dim(nmb3[[1]])[2]),
    matrix(rep(0,length(nmb3[[1]])), ncol=dim(nmb3[[1]])[2]),
    matrix(rep(0,length(nmb3[[1]])), ncol=dim(nmb3[[1]])[2]),
    matrix(rep(0,length(nmb3[[1]])), ncol=dim(nmb3[[1]])[2]),
    matrix(rep(0,length(nmb3[[1]])), ncol=dim(nmb3[[1]])[2])
  )
    
 
  x<-NULL
  #ex.zeros[nmb3 > 0]<-1
  for(xxx in 1:7) {ex.zeros[[xxx]][nmb3[[xxx]] > 0]<-1}
  nmb4<-lapply(1:7, function(x) colMeans(ex.zeros[[x]]))
  
  nmb5<-do.call(rbind, nmb4)
  
  curve.data.strat<-data.frame(cbind(rep('GB', length(cea.tables)), inv.names, rep(threshold.val,length(cea.tables)), nmb5))

  return(curve.data.strat)
}
  

accept.strat.graph<-function(strain, cov.var, dcount)
{
  ss<-seq(1000,30000,by = 500)
  for(hh in 1:length(ss)) {ifelse(hh==1, accept.out<-prob.curve(ss[hh], strain, cov.var, dcount), 
                                  accept.out<-rbind(accept.out, prob.curve(ss[hh], strain, cov.var, dcount)))}
  
  colnames(accept.out)<-c('country','intervention','threshold' ,'0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR', 'total')

  mj<-ggplot(data=accept.out, aes(x=as.numeric(as.character(threshold)), y=as.numeric(as.character(total)), color=intervention, linetype=intervention))+geom_line(aes(color=intervention, group=intervention))
  
  label<-paste(paste0('strain=', strain.name[strain]),paste0('discount=', dcount.l[dcount], '%'), paste0('coverage=', cov.f[cov.var]))
  
  AC3<-mj+labs(title=NULL, subtitle=label, fill='Intervention')+geom_line(aes(color=intervention))+scale_y_continuous(limits = c(0,1), breaks=scales::pretty_breaks(n=7))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_vline(xintercept = c(15000,20000), linetype=3)+theme(plot.subtitle=element_text(size=10, hjust=0.5, color="black"))+guides(fill=guide_legend(title="Intervention"))
  
return(AC3)
}


#arrange output
accept3<-NULL
accept3<-lapply(1:3, function(dcount) ggarrange(
  accept.strat.graph(1,1,dcount)+rremove("xy.title"), accept.strat.graph(2,1,dcount)+rremove("xy.title"), accept.strat.graph(3,1,dcount)+rremove("xy.title"), 
  accept.strat.graph(1,2,dcount)+rremove("xy.title"), accept.strat.graph(2,2,dcount)+rremove("xy.title"), accept.strat.graph(3,2,dcount)+rremove("xy.title"),
  accept.strat.graph(1,3,dcount)+rremove("xy.title"), accept.strat.graph(2,3,dcount)+rremove("xy.title"), accept.strat.graph(3,3,dcount)+rremove("xy.title"),
  ncol=3, nrow=3, common.legend = T, legend='bottom'))

graphs.acpost<-NULL
graphs.acpost<-marrangeGrob(accept3, nrow = 1, ncol=1, 
                            top = text_grob('Acceptability Curve per Strategy'),
                            bottom = text_grob('Willingness-to-pay Threshold (£)'),
                            left = text_grob('Probability Cost-Effective', rot=90))

ggsave(paste0("AC3strainvcoverage",version,".pdf"),plot=graphs.acpost, width=12, height=9, units='in',device='pdf')

###################################################
#All strains
#########################################################################


#reset directory
setwd(paste0("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/GB/coverage/0.55"))

function(fstrain, new.cov, dcount)
{
alt.program.loader(strain = 1, i.cov = 0.55, dcount = 3.5)
alt.program.loader(strain = 2, i.cov = 0.55, dcount = 3.5)
alt.program.loader(strain = 3, i.cov = 0.55, dcount = 3.5)

for(fstrain in 1:3)
{
  if(fstrain==1) {fforgraph<-strategy.summaryH1N1}
  if(fstrain==2) {fforgraph<-strategy.summaryH3N2}
  if(fstrain==3) {fforgraph<-strategy.summaryB}
  
    choice<-(seq(1,dim(fforgraph)[1],by = 2))
    #Status Quo
    sq.int<-fforgraph[choice,]
    sq<-sq.int[1:12,]
    sq$variable <- factor(sq$variable, levels=c(paste0('cases',1), paste0('deaths',1), paste0('GP consults',1), paste0('hospitalized',1)), labels = c("Incidence", "Mortality",'GP Consults','Hospitalizations'))
    sq$strategy<-'Status Quo'
    sq.clinical<-na.omit(sq)
    
    
   #Interventions
    int<-fforgraph[choice+1,]
    int$variable <- factor(int$variable, levels=c(paste0('cases',2), paste0('deaths',2), paste0('GP consults',2), paste0('hospitalized',2)), labels = c("Incidence", "Mortality",'GP Consults','Hospitalizations'))
    int.clinical<-na.omit(int)
    
    clinical<-rbind(int.clinical, sq.clinical)
    
    #clinical graph
    p <-ggplot(clinical, aes(x=strategy, y=value, fill=strategy))+geom_bar(position = "dodge", stat = "identity", aes(fill=strategy))+facet_wrap(~variable, ncol=2, scales = 'free')+labs(title = paste(strain.name[fstrain],'Clinical Outcomes'), x = 'Strategy', y = 'Annual Incidence', fill = "Intervention\n Strata") +geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2, position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = rev(cols[1:8]))+theme(
      axis.text.x = element_blank(), axis.ticks = element_blank())

    #costs graphs
    sq.costs<-sq.int[1:12,]
    sq.costs$variable <- factor(sq.costs$variable, levels=c(paste0('vaccine cost',1), paste0('Healthcare costs',1), paste0('GP costs',1), paste0('hospital costs',1),paste0('QALYs lost',1),paste0('case QALYs',1), paste0('death QALYs',1), paste0('nondeath QALYs', 1)), labels = c('Vaccination Costs','Healthcare Costs','GP Consult Fees','Hospitalization Costs','Years of Life Lost', 'YLL (cases)', 'YLL (deaths)','YLL (nondeath)'))
    sq.costs$strategy<-'Status Quo'
    sq.costs<-na.omit(sq.costs)
    
    int.costs<-fforgraph[choice+1,]
    int.costs$variable <- factor(int.costs$variable, levels=c(paste0('vaccine cost',2), paste0('Healthcare costs',2), paste0('GP costs',2), paste0('hospital costs',2),paste0('QALYs lost',2),paste0('case QALYs',2), paste0('death QALYs', 2), paste0('nondeath QALYs', 2)), labels = c('Vaccination Costs','Healthcare Costs','GP Consult Fees','Hospitalization Costs','Years of Life Lost', 'YLL (cases)', 'YLL (deaths)','YLL (nondeath)'))
    int.costs<-na.omit(int.costs)
    
    costs<-rbind(int.costs, sq.costs)
    
    q <-ggplot(costs, aes(x=strategy, y=value, fill=strategy))+geom_bar(position = "dodge", stat = "identity", aes(fill=strategy))+facet_wrap(~variable, ncol=3, scales = 'free')+labs(title = paste(strain.name[fstrain], 'Economic Outcomes'), x = 'Strategy', y = 'Costs(£)', fill = "Economic Impact\n") +geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2, position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = rev(cols[1:8]))+theme(
      axis.text.x = element_blank(), axis.ticks = element_blank())

    #sn<-c('SQ','Sec','All3')
    #print(p)
    ggsave(paste0(strain.name[fstrain],new.cov,'cloutcomes',dcount,".png"),plot = p, width=13, units='in',
           device='png');
    ggsave(paste0(strain.name[fstrain],new.cov,'costs',dcount,".png"),plot = q, width=13, units='in',device='png')
    
    ifelse(fstrain==1, all.costs<-costs, all.costs<-rbind(all.costs, costs))
    ifelse(fstrain==1, all.clinical<-clinical, all.clinical<-rbind(all.clinical, clinical))
}


a3<-ggplot(all.costs, aes(fill=strain, y=value, x=strategy)) +
  geom_bar(position="dodge", stat="identity")+facet_wrap(~variable, ncol=3, scales = 'free')+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2, position=position_dodge(.9))+theme(text = element_text(size=10),axis.text.x = element_text(angle=20, hjust=1)) 

ggsave(paste0('allstrains',new.cov,'costs',dcount,".png"),plot = a3, width=13, units='in',
       device='png');

c3<-ggplot(all.clinical, aes(fill=strain, y=value, x=strategy)) +
  geom_bar(position="dodge", stat="identity")+facet_wrap(~variable, ncol=3, scales = 'free')+geom_errorbar(aes(ymin=value-SE, ymax=value+SE), width=.2, position=position_dodge(.9))+theme(text = element_text(size=10),axis.text.x = element_text(angle=20, hjust=1)) 

ggsave(paste0('allstrains',new.cov,'clinical',dcount,".png"),plot = c3, width=13, units='in',
       device='png');
}




#######################################################################################
####### 8 GRAPH JAN Style Plot with age groups ########
####################################################################################################


setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB','coverage', cov.f[cov.var]))

#make status quo category

#cols<-my.col

inv.costs.f$strategy<-factor(inv.costs.f$strategy, levels=c('Preschool','Primary','Secondary','Preschool+Primary','PrePrimeSeconday', 'Preschool+Secondary','Primary+Secondary', 'Status Quo'), labels=c( "S1: 2-4 yrs old"  ,                 "I2: 5-11 yrs old"       ,              "I3: 12-16 yrs old"      ,             "I4: 2-11 yrs old" ,"I7: 2-16 yrs old", "I5: 2-4 & 12-16 yrs old"   ,      "I6: 5-16 yrs old", 'SQ: Risk Groups, 65+ yrs old'))

q<-ggplot(inv.costs.f, aes(x=strategy, y=value, fill=strategy))+geom_bar(position = "dodge", stat = "identity", aes(fill=strategy))+facet_wrap(~variable, ncol=3, scales = 'free')+scale_y_continuous(breaks=scales::pretty_breaks(n=7))+labs(title = NULL,subtitle=paste('discount=', dcount.l[dcount], '%', ',','coverage=', cov.f[cov.var],',', 'strain=all'), x = 'Intervention Strategy', y = '"Average Annual Costs(£)"', fill = "Economic Impact\n") +geom_errorbar(aes(ymin=CIL, ymax=CIU), width=.2, position=position_dodge(.9))+theme_bw()+theme(
  axis.text.x = element_blank(), axis.ticks = element_blank(), plot.subtitle=element_text(size=10, hjust=0.5, color="black"), legend.position="bottom")


ggsave(paste0('costsplot',cov.f[cov.var],dcount,".png"),plot = q, width=11, height=7.5, units='in',device='png')
return(q)


for(tt in 1:2)
{
  discount.loader2(cov.var, dcount)
  inv.costs<-get(paste0('inv.costs',dcount.l[dcount]))
  
  
  inv.costs$variable <- factor(inv.costs$variable, levels=c(paste0('vaccine cost',tt), paste0('Healthcare costs',tt), paste0('GP costs',tt), paste0('hospital costs',tt),paste0('QALYs lost',tt),paste0('QALYs lost case',tt), paste0('QALYs lost death',tt), paste0('QALYs lost hosp', tt)), labels = c('Vaccination Costs','Healthcare Costs','GP Consult Fees','Hospitalization Costs','Years of Life Lost (All Sources)', 'Years of Life Lost (Symptomatic Cases)', 'Years of Life Lost (Deaths)','Years of Life Lost (Hospitalization)'))
  if(tt==1) inv.costs$strategy<-'Status Quo'
  ifelse(tt==1, inv.costs.f<-na.omit(inv.costs)[1:6,], inv.costs.f<-rbind(inv.costs.f, na.omit(inv.costs)))
}


for(fstrain in 1:3)
{
  for(intv in 1:length())
  
    
    #add 7 age groups
    
    choice<-(seq(1+pp,dim(fagetab)[1],by = 3))
    int<-fagetab[choice,]
    #int$strain <- factor(int$strain, levels=c('H1N1','H3N2','B'), labels = c(1:3))
    int<-int[int$strain==fstrain,]
    int$type <- factor(int$type, levels=c(paste0('Healthcare costs',pp+1), paste0('GP costs',pp+1), paste0('hospital costs',pp+1),paste0('QALYs lost',pp+1)), labels = c('Healthcare Costs','GP Consult Fees','Hospitalization Costs','Years of Life Lost'))
    int<-int[order(int$type),]
    
    p <-ggplot(int, aes(x=as.factor(country), y=as.numeric(value), fill=type))+theme_bw()+geom_bar(stat = "identity", aes(fill=variable))+facet_wrap(~type, ncol=2, scales = 'free')+labs(title = paste(strain.name[fstrain],inv.names[pp]), x = 'Coverage', y ='Costs(£)', fill = "Economic Impact\n") +scale_fill_manual(values = cols[1:12])
    
    
    #--------Age statified Outcomes under different interventions
    #choice2<-(seq(1+pp,dim(fageout)[1],by = 3))
    
    discount.loader(strain,cov.var, dcount)
    inv.costs<-get(paste0('inv.costs', strain.name[strain], dcount))
    
    
    outs<-fageout[choice2,]
    outs$strain <- factor(outs$strain, levels=c('H1N1','H3N2','B'), labels = c(1:3))
    outs<-outs[outs$strain==fstrain,]
    outs$ctype <- factor(outs$ctype, levels=c(paste0('cases',pp+1), paste0('deaths',pp+1), paste0('GP consults',pp+1), paste0('hospitalized',pp+1)), labels = c("Incidence", "Mortality",'GP Consults','Hospitalizations'))
    
    q <-ggplot(outs, aes(x=as.factor(country), y=as.numeric(value), fill=ctype))+geom_bar(stat = "identity", aes(fill=variable))+facet_wrap(~ctype, ncol=2, scales = 'free')+labs(title = paste(strain.name[fstrain],interventions[pp+1]), x = 'Countries', y = 'Counts', fill = "Clinical Impact\n") +theme_bw()+scale_fill_manual(values = cols[1:12])
    
    sn<-c('SQ','Sec','All3')
    #print(p)
    ggsave(paste0(strain.name[fstrain],sn[pp+1],'agestratcosts',".png"),plot = p, width=13, units='in',
           device='png');
    ggsave(paste0(strain.name[fstrain],sn[pp+1],'agestratoutcomes',".png"),plot = q, width=13, units='in',device='png')
  }
}

#grouped

########################################################################################
####### 8 GRAPH JAN Style Plot with age groups ########
####################################################################################################

for(fstrain in 1:3)
{
  for(pp in 0:2)
  {
    p<-NULL
    interventions<-c('Status Quo Vaccination','12-16 Year Old Vaccination','2-16 Year Old Vaccination')
    
    choice<-(seq(1+pp,dim(fagetab)[1],by = 3))
    int<-fagetab[choice,]
    int$strain <- factor(int$strain, levels=c('H1N1','H3N2','B'), labels = c(1:3))
    int<-int[int$strain==fstrain,]
    int$type <- factor(int$type, levels=c(paste0('Healthcare costs',pp+1), paste0('GP costs',pp+1), paste0('hospital costs',pp+1),paste0('QALYs lost',pp+1)), labels = c('Healthcare Costs','GP Consult Fees','Hospitalization Costs','Years of Life Lost'))
    int<-int[order(int$type),]
    
    p <-ggplot(int, aes(x=as.factor(country), y=as.numeric(value), fill=type))+theme_bw()+geom_bar(stat = "identity", aes(fill=variable))+facet_wrap(~type, ncol=2, scales = 'free')+labs(title = paste(strain.name[fstrain],interventions[pp+1]), x = 'Countries', y ='Costs(£)', fill = "Economic Impact\n") +scale_fill_manual(values = cols[1:12])
    
    
    #--------Age statified Outcomes under different interventions
    choice2<-(seq(1+pp,dim(fageout)[1],by = 3))
    outs<-fageout[choice2,]
    outs$strain <- factor(outs$strain, levels=c('H1N1','H3N2','B'), labels = c(1:3))
    outs<-outs[outs$strain==fstrain,]
    outs$ctype <- factor(outs$ctype, levels=c(paste0('cases',pp+1), paste0('deaths',pp+1), paste0('GP consults',pp+1), paste0('hospitalized',pp+1)), labels = c("Incidence", "Mortality",'GP Consults','Hospitalizations'))
    
    q <-ggplot(outs, aes(x=as.factor(country), y=as.numeric(value), fill=ctype))+geom_bar(stat = "identity", aes(fill=variable))+facet_wrap(~ctype, ncol=2, scales = 'free')+labs(title = paste(strain.name[fstrain],interventions[pp+1]), x = 'Countries', y = 'Counts', fill = "Clinical Impact\n") +theme_bw()+scale_fill_manual(values = cols[1:12])
    
    sn<-c('SQ','Sec','All3')
    #print(p)
    ggsave(paste0(strain.name[fstrain],sn[pp+1],'agestratcosts',".png"),plot = p, width=13, units='in',
           device='png');
    ggsave(paste0(strain.name[fstrain],sn[pp+1],'agestratoutcomes',".png"),plot = q, width=13, units='in',device='png')
  }
}

########################################################################################
####### 9 GRAPH Age structured outcomes summed across all strain means ########
####################################################################################################

setwd(paste0("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/GB"))


for(pp in 0:2)
{
  r<-NULL
  interventions<-c('Status Quo Vaccination','12-16 Year Old Vaccination','2-16 Year Old Vaccination')
  choice<-(seq(1+pp,dim(allfagetab)[1],by = 3))
  int<-allfagetab[choice,]
  int$type <- factor(int$type, levels=c(paste0('Healthcare costs',pp+1), paste0('GP costs',pp+1), paste0('hospital costs',pp+1),paste0('QALYs lost',pp+1)), labels = c('Healthcare Costs','GP Consult Fees','Hospitalization Costs','Years of Life Lost'))
  int<-int[order(int$type),]
  
  r <-ggplot(int, aes(x=as.factor(country), y=as.numeric(value), fill=type))+theme_bw()+geom_bar(stat = "identity", aes(fill=variable))+facet_wrap(~type, ncol=2, scales = 'free')+labs(title = paste('All Strains',interventions[pp+1]), x = 'Countries', y ='Costs(£)', fill = "Economic Impact\n") +scale_fill_manual(values = cols[1:12])
  
  
  #--------Age statified Outcomes under different interventions
  choice2<-(seq(1+pp,dim(allfageout)[1],by = 3))
  outs<-allfageout[choice2,]
  outs$ctype <- factor(outs$ctype, levels=c(paste0('cases',pp+1), paste0('deaths',pp+1), paste0('GP consults',pp+1), paste0('hospitalized',pp+1)), labels = c("Incidence", "Mortality",'GP Consults','Hospitalizations'))
  
  s <-ggplot(outs, aes(x=as.factor(country), y=as.numeric(value), fill=ctype))+geom_bar(stat = "identity", aes(fill=variable))+facet_wrap(~ctype, ncol=2, scales = 'free')+labs(title = paste('All Strains',interventions[pp+1]), x = 'Countries', y = 'Counts', fill = "Clinical Impact\n") +theme_bw()+scale_fill_manual(values = cols[1:12])
  
  sn<-c('SQ','Sec','All3')
  #print(p)
  ggsave(paste0('Allstrains',sn[pp+1],'agestratcosts',".png"),plot = r, width=13, units='in',
         device='png');
  ggsave(paste0('Allstrains',sn[pp+1],'agestratoutcomes',".png"),plot = s, width=13, units='in',device='png')
  
}


##############################################
#CEP


##########################################
####Acceptability curve by AGE GROUPS
#################################################

#inv.names2<-c('Status Quo','Preschool', 'Primary', 'Secondary')

setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')

accept.curve<-function(threshold.val, qaly, cost)
{
  qaly2<- Map('*', -1, Map('*',qaly,threshold.val))
  nmb<-Map('-',qaly2, cost)
  #nmb3<-Reduce('+', nmb)/tab.diff
  
  return(nmb)
}


c.names<-c('program','0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR', 'total')

outerf<-list()
accept.out1<-accept.out2<-accept.out3<-accept.out4<-list()
for(i.country in 1:length(fname))
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  qal.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.qaly.loss', strain.name[strain],'all.programs', '*', 'd',dcount)))
  cos.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.cost.loss', strain.name[strain],'all.programs', '*', 'd',dcount)))
  
  
  load(qal.tables[[1]])
  load(cos.tables[[1]])
  
  ss<-seq(0,30000,by = 300)
  for(hh in 1:length(ss)) {
    
    
    accept.out1<-accept.curve(ss[hh], tot.QALY.loss[[1]], tot.costs.loss[[1]]) 
    accept.out2<-accept.curve(ss[hh],tot.QALY.loss[[2]], tot.costs.loss[[2]])
    accept.out3<-accept.curve(ss[hh],tot.QALY.loss[[3]], tot.costs.loss[[3]])
    accept.out4<-accept.curve(ss[hh],tot.QALY.loss[[4]], tot.costs.loss[[4]])
    
    pat1<-data.frame(as.factor(rep(inv.names2[1],n.samples)), Reduce('+',accept.out1)/tab.diff, rowSums(Reduce('+',accept.out1)/tab.diff))
    pat2<-data.frame(as.factor(rep(inv.names2[2],n.samples)), Reduce('+',accept.out2)/tab.diff, rowSums(Reduce('+',accept.out2)/tab.diff))
    pat3<-data.frame(as.factor(rep(inv.names2[3],n.samples)), Reduce('+',accept.out3)/tab.diff, rowSums(Reduce('+',accept.out3)/tab.diff))
    pat4<-data.frame(as.factor(rep(inv.names2[4],n.samples)), Reduce('+',accept.out4)/tab.diff, rowSums(Reduce('+',accept.out4)/tab.diff))
    
    colnames(pat1)<-c.names
    
    
    colnames(pat2)<-colnames(pat1)
    colnames(pat3)<-colnames(pat1)
    colnames(pat4)<-colnames(pat1)
    
    #am.1[[hh]]<-pat1
    #am.2[[hh]]<-pat2
    #am.3[[hh]]<-pat3
    #am.4[[hh]]<-pat4
    
    
    for(aa in 2:dim(pat1)[2]){
      e.frame<-data.frame(pat1[,aa], pat2[,aa],pat3[,aa],pat4[,aa])
      xx<-apply(e.frame, 1, function(z) {which.max(z)})  
      ex.frame<-matrix(rep(0,dim(e.frame)[1]*dim(e.frame)[2]), ncol=4)
      for(cc in 1:length(xx))
      {ex.frame[cc,xx[cc]]<-1}
      
      
      outp<-data.frame(fname[i.country], c.names[aa], ss[hh], t(colMeans(ex.frame)))
      
      
      ifelse(aa==2&hh==1&i.country==1, outerf<-outp, outerf<-rbind(outerf, outp))
    }
    
  }
}


setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')
colnames(outerf)<-c('country', 'Age','threshold','status quo','preschool','primary','secondary')


outer.g<-melt(outerf, id.vars=c('country','Age', 'threshold'))
m<-ggplot(data=outer.g, aes(x=as.numeric(as.character(threshold)), y=as.numeric(as.character(value)), color=variable, linetype=variable))  
pdf(file='OC.3program.pdf',onefile = T, height=11, width=8.75)
for(gg in 1:32){
  print(m+ facet_wrap_paginate(country~Age, ncol = 2, nrow = 4, scales='free', page=gg)+
          geom_line(aes(color=variable, group=variable, linetype=variable))+labs(title='Probability a strategy is Optimal', x='Willing-to-pay Threshold (£)', y='Probability Cost-Effective',color='Strategy') +
          scale_y_continuous(breaks=scales::pretty_breaks(n=8))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_vline(xintercept = c(15000,20000), linetype=3))
}
dev.off()


#totals only

outer.g<-melt(outerf, id.vars=c('country','Age', 'threshold'))
mk<-ggplot(data= outer.g[outer.g$Age=='total',], aes(x=as.numeric(as.character(threshold)), y=as.numeric(as.character(value)), color=variable, linetype=variable))  
pdf(file='OC.totals.3program.pdf',onefile = T, height=11, width=8.75)
for(gg in 1:2){
  print(mk+ facet_wrap_paginate(country~Age, ncol = 2, nrow = 3, scales='free', page=gg)+
          geom_line(aes(color=variable, group=variable, linetype=variable))+labs(title='Optimal Acceptability Curve per Strategy', x='Willingness-to-pay Threshold (£)', y='Probability Cost-Effective',color='Strategy') +
          scale_y_continuous(breaks=scales::pretty_breaks(n=8))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_vline(xintercept = c(15000,20000), linetype=3))
}
dev.off()



###GB only


gb.grob<-outer.g[outer.g$Age=='total'&outer.g$country=='GB',]

col<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

oc.gb<-ggplot(data= gb.grob, aes(x=as.numeric(as.character(threshold)), y=as.numeric(as.character(value)), group=variable,color=variable))+scale_color_manual(values=c(col))+
  geom_line()+labs(title='Optimal Acceptability Curve \n per Strategy', x='Willingness-to-pay Threshold (£)', y='Probability Cost-Effective',color='Strategy') +
  scale_y_continuous(breaks=scales::pretty_breaks(n=8))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_vline(xintercept = c(15000,20000), linetype=3)+ theme_bw()+
  annotate("text", x = 29000, y = 0.05, label = paste0('GB'))


dev.off()



##########################################################################################################
# ACCEPTIBILITY CURVE for single intervention probability it is cost effective on WTP scale
#################################################################################################################################
prob.curve<-function(threshold.val, int, i.country)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  cea.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'ICER.strat', strain.name[strain],'p', '*', 'd',dcount)))
  
  
  load(cea.tables[[int-1]]);
  ICER.set.n<-get(paste0('ICER.set', int))
  
  
  qal2<-Map('*', ICER.set.n$`QALY diff`, threshold.val)
  nmb<-Map('-',qal2, ICER.set.n$`cost diff`)
  nmb2<-as.matrix(Reduce('+',nmb)/tab.diff)
  nmb3<-cbind(nmb2, rowSums(nmb2))
  
  #nmb2<-Reduce('+', ICER.set2$`icer`)/tab.diff
  
  ex.zeros<-matrix(rep(0,length(nmb3)), ncol=dim(nmb3)[2])
  
  ex.zeros[nmb3 > 0]<-1  
  
  curve.data.strat<-c(as.character(fname[i.country]), threshold.val, as.numeric(colMeans(ex.zeros)))
  
  return(curve.data.strat)
}

outer<-list()
for(i.country in 1:length(fname))
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])))
  
  ss<-seq(1000,30000,by = 1000)
  for(hh in 1:length(ss)) {ifelse(hh==1, accept.out2<-prob.curve(ss[hh], 2, i.country ), 
                                  accept.out2<-rbind(accept.out2, prob.curve(ss[hh], 2, i.country )))}
  
  for(hh in 1:length(ss)) {ifelse(hh==1, accept.out3<-prob.curve(ss[hh], 3, i.country ), 
                                  accept.out3<-rbind(accept.out3, prob.curve(ss[hh], 3, i.country )))}
  
  for(hh in 1:length(ss)) {ifelse(hh==1, accept.out4<-prob.curve(ss[hh], 4, i.country ), 
                                  accept.out4<-rbind(accept.out4, prob.curve(ss[hh], 4, i.country )))}
  
  
  colnames(accept.out2)<-c('country','threshold' ,'0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR', 'total')
  
  #colnames(accept.out2)<-c('country','threshold' ,'0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR', 'total')
  
  colnames(accept.out2)<-colnames(accept.out3)
  colnames(accept.out2)<-colnames(accept.out4)
  
  outer[[i.country]]<-list(accept.out2, accept.out3, accept.out4)
}

names(outer)<-fname[1:10]

#install.packages("ggthemes") # Install 


for(gg in 1:length(outer))
{
  for(hh in 1:length(outer[[4]]))
  {
    dtm<-outer[[gg]][[hh]]
    
    dtm<-data.frame(dtm)
    colnames(dtm)<-c('country','threshold' ,'0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR', 'total')
    
    
    accept.m<-melt(dtm,id.vars = c('country', 'threshold'))
    program<-rep(inv.names2[hh+1], dim(accept.m)[1])
    
    accept.r<-cbind(accept.m, program)
    
    ifelse(hh==1, accept.out<-accept.r, accept.out<-rbind(accept.out, accept.r))
    #ifelse(gg==1&hh==1, accept.out<-accept.r, accept.out<-rbind(accept.out, accept.r))
  }
  
}


pdf(file='AC.countrybyprogram.pdf',onefile = T, height=11, width=8.75)
for(gg in 1:18){
  m<-ggplot(data=accept.out, aes(x=as.numeric(as.character(threshold)), y=as.numeric(as.character(value)), color=country, linetype=country)) 
  
  print(m+ facet_wrap_paginate(program~variable, ncol = 2, nrow = 3, scales='free', page=gg)+
          geom_line(aes(color=country, group=country))+labs(title='Acceptability Curve per Strategy', x='Willing-to-pay Threshold (£)', y='Probability Cost-Effective',color='Strategy') +
          scale_y_continuous(breaks=scales::pretty_breaks(n=8))+scale_x_continuous(breaks=scales::pretty_breaks(n=10))+geom_vline(xintercept = c(15000,20000), linetype=3))
}
dev.off()

pdf(file='AC.total.3program.pdf',onefile = T, height=11, width=8.75)
m<-ggplot(data=accept.out[accept.out$variable=='total',], aes(x=as.numeric(as.character(threshold)), y=as.numeric(as.character(value)), color=program, linetype=program)) 

for(gg in 1:32){
  
  print(m+facet_wrap_paginate(country~variable, ncol = 2, nrow = 3, scales='free', page=gg)+
          geom_line(aes(color=program, group=program))+labs(title='Acceptability Curve per Strategy', x='Willingness-to-pay (£)', y='Probability Cost-Effective',color='Strategy') +
          scale_y_continuous(breaks=scales::pretty_breaks(n=8))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_vline(xintercept = c(15000,20000), linetype=3))
}
dev.off()


ggsave(paste0('ACsingle',dcount,".png"),plot = AC.single, width=10, units='in',device='png')

####

ac.gb<-ggplot(data=accept.out[accept.out$variable=='total',], aes(x=as.numeric(as.character(threshold)), y=as.numeric(as.character(value)), color=program, group=program))+scale_color_manual(values=col[2:4])+
  geom_line()+labs(title='Acceptability Curve \n per Strategy', x='Willingness-to-pay (£)', y='Probability Cost-Effective',color='Strategy') +
  scale_y_continuous(breaks=scales::pretty_breaks(n=8))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_vline(xintercept = c(15000,20000), linetype=3)+theme_bw()+theme(legend.position="none")+
  annotate("text", x = 29000, y = 0.05, label = paste0('GB'))

nl<-get_legend(oc.gb)

pg<-plot_grid(oc.gb+theme(legend.position="none"), ac.gb,nrow=1, ncol=2, labels = 'AUTO', align = 'h', label_size = 12)
plot_grid(pg, nl, rel_widths = c(2, .3))
ggsave(paste0('GB.CEA.curves','.jpg'),device='jpg',plot=last_plot(), width=21,height=10,units='cm')


#########Calculating Cost-Effectiveness Average for the Table
kerneld <- kde2d(rowSums(check.p[[1]]), rowSums(check.p[[2]]), n=500, lims=c(0,limit,0,100))
#check[[1]]<-check[[1]][check[[1]]>=0]
#check[[2]]<-check[[2]][check[[2]]>=0]
pp <- array()
for (i in 1:500){
  z.x <- max(which(kerneld$x < rowSums(check.p[[1]])[i]))
  z.y <- max(which(kerneld$y < rowSums(check.p[[2]])[i]))
  pp[i] <- kerneld$z[z.x, z.y]
}
confidencebound <- quantile(pp, 0.1, na.rm = TRUE)
cc<-Map('-',Map('*', -1, tot.costs.loss[[2]]),Map('*', -1, tot.costs.loss[[1]]))
for(yyy in 1:length(cc))
{cc[[yyy]][cc[[yyy]]==0]<-1e-5 }


{ICER.set2$`cost diff`[[yyy]][ICER.set2$`cost diff`[[yyy]]==0]<-1e-5}
rowSums(Reduce('+',)/length(ICER.set2$icer))

mean(unlist(lapply(ICER.set2$`cost diff`,rowSums)))/mean(unlist(lapply(ICER.set2$`QALY diff`,rowSums)))
intl.ICER.loader(1,1,3.5)


mean(unlist(lapply(Map('-',tot.costs.loss[[4]], tot.costs.loss[[2]]), rowSums)))/mean(unlist(lapply(Map('-',tot.QALY.loss[[4]], tot.QALY.loss[[2]]), rowSums)))

mean(unlist(lapply(Map('-',tot.costs.loss[[2]], tot.costs.loss[[1]]), rowSums)))/mean(unlist(lapply(Map('-',tot.QALY.loss[[2]], tot.QALY.loss[[1]]), rowSums)))

mean(unlist(lapply(Map('-',tot.costs.loss[[3]], tot.costs.loss[[1]]), rowSums)))/mean(unlist(lapply(Map('-',tot.QALY.loss[[3]], tot.QALY.loss[[1]]), rowSums)))

mean(unlist(lapply(Map('-',tot.costs.loss[[4]], tot.costs.loss[[1]]), rowSums)))/mean(unlist(lapply(Map('-',tot.QALY.loss[[4]], tot.QALY.loss[[1]]), rowSums)))


#average total cost lost due to healthcare and vaccine spent costs
mean(unlist(lapply(lapply(tot.costs.loss[[1]], rowSums), mean)))
mean(unlist(lapply(lapply(tot.costs.loss[[2]], rowSums), mean)))
mean(unlist(lapply(lapply(tot.costs.loss[[3]], rowSums), mean)))
mean(unlist(lapply(lapply(tot.costs.loss[[4]], rowSums), mean)))

#average total qalys lost
mean(unlist(lapply(lapply(tot.QALY.loss[[1]], rowSums), mean)))
mean(unlist(lapply(lapply(tot.QALY.loss[[2]], rowSums), mean)))
mean(unlist(lapply(lapply(tot.QALY.loss[[3]], rowSums), mean)))
mean(unlist(lapply(lapply(tot.QALY.loss[[4]], rowSums), mean)))



q.gained<-Map('-',Map('*', -1, tot.QALY.loss[[2]]),Map('*', -1, tot.QALY.loss[[1]]))

mean(lapply(rowSums, function(x) Map('/',q.gained,cc)[[x]]))

Map('-',tot.costs.loss[[2]],tot.costs.loss[[1]])
