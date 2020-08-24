
#install.packages("odin")
library(devtools)
#install_github("MJomaba/flu-evidence-synthesis", dependencies = TRUE, build_vignettes = F)
#
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
library( ggpubr )
rm(list = ls())

#grew<-data.frame(cbind(seq(40,90,10), c(20690, 9270, 3461,1039,205,14), c(12528, 5627, 2128, 614, 69, -53)))
#names(grew)<-c('Age','Male','Female')
#grew2<-melt(grew, id.vars='Age')

#ggplot(grew2, aes(x=Age, y=value, color=variable))+geom_point()+geom_line()+ylab('Incremental Net Benefit (£)')+labs(color='Gender')

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
#source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/INTL_SENS files/OUTPUT_ICER_conversion.R')
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/UK Files/FUNC_odin_UK.R')
source('/Users/Natasha/Dropbox/UKfluworkGIT/Outputfunctions/UK Files/OUTPUT_icer_conversion_uk.R')

#remainder items
y.name<<-c(1995:2014)

setwd('/Users/Natasha/Dropbox/UKfluworkGIT/DavidsCode/RSave/')
death.risk.tables<-list.files(pattern=glob2rx('tab_risk_death_*')); #loads in the order B,H1,H3
hosp.risk.tables<-list.files(pattern=glob2rx('tab_risk_hosp*'));
GP.risk.tables<-list.files(pattern=glob2rx('tab_risk_GP*'))

setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputfunctions', 'UK Files'));
fname<-c('BE', 'DE', 'FI', 'GB', 'IT', 'LU', 'NL', 'PL', 'PE', 'FR', 'ZI')
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
  load(sens.tables[1]); avg.incidence5<-total.keep; #preprime
  load(sens.tables[5]); avg.incidence6<-total.keep; #preprimesec
  load(sens.tables[6]); avg.incidence7<-total.keep; #secondary
  load(sens.tables[7]); avg.incidence8<-total.keep; #secondary
  
  
  var.name<-paste0('uk.table',strain.name[strain])
  #('WTP','Preschool','Primary','Secondary','Preschool+Primary','PrePrimeSeconday','Preschool+Secondary','Primary+Secondary')
  assign(var.name, list(avg.incidence1,avg.incidence2,avg.incidence3, avg.incidence4, avg.incidence5,
                        avg.incidence6, avg.incidence7,avg.incidence8
  ), envir = .GlobalEnv)
}


###########################################################################################
#Major Runs
################################################################################################


####### All Strains ICER calculation ########
#cov.level<-list(thirty.cov, fiftyfive.cov, seventy.cov)
cov.level<-c(mL3, mL55, mL7)
dcount.l<-c(0,1.5,3.5)
cov.f<-c(0.3, 0.55, 0.7)
#cov.var<-1

  GP.costs.l<-Hosp.costs.l<-sum.costs.l<-aged.cases.l<-QALY.tot.l<-QALY.death.l<-cumi.l<-add.dose<-avt.dose<-list()
  
  for(dd in 1:length(dcount.l))
  {
    for(new.strat in 2:8)
    {
    #new.strat<-
      cov.sens.loader(strain=2, end.cov=cov.f[2], version='v1')
      cov.sens.loader(strain=3, end.cov=cov.f[2], version='v1')
      cov.sens.loader(strain=1, end.cov=cov.f[2], version='v1')
      
      uk.CEA.allstrains(Dataset1H1N1 = uk.tableH1N1[[1]], Dataset2H1N1 = uk.tableH1N1[[new.strat]],
                        Dataset1H3N2 = uk.tableH3N2[[1]], Dataset2H3N2 = uk.tableH3N2[[new.strat]],
                        Dataset1B = uk.tableB[[1]], Dataset2B = uk.tableB[[new.strat]],
                        dcount = dd, program1=1, program2=new.strat, i.cov=2, new.cov=mL55, version='v1')
    }
  }
  
  
  
  GP.costs.l<-Hosp.costs.l<-sum.costs.l<-aged.cases.l<-QALY.tot.l<-QALY.death.l<-cumi.l<-add.dose<-avt.dose<-list()
  
  for(dd in 1:length(dcount.l))
  {
    for(new.strat in 2:8)
    {
      #new.strat<-
      cov.sens.loader(strain=2, end.cov=cov.f[3], version='v1')
      cov.sens.loader(strain=3, end.cov=cov.f[3], version='v1')
      cov.sens.loader(strain=1, end.cov=cov.f[3], version='v1')
      
      uk.CEA.allstrains(Dataset1H1N1 = uk.tableH1N1[[1]], Dataset2H1N1 = uk.tableH1N1[[new.strat]],
                        Dataset1H3N2 = uk.tableH3N2[[1]], Dataset2H3N2 = uk.tableH3N2[[new.strat]],
                        Dataset1B = uk.tableB[[1]], Dataset2B = uk.tableB[[new.strat]],
                        dcount = dd, program1=1, program2=new.strat, i.cov=3, new.cov=mL7, version='v1')
    }
  }
  

  
  GP.costs.l<-Hosp.costs.l<-sum.costs.l<-aged.cases.l<-QALY.tot.l<-QALY.death.l<-cumi.l<-add.dose<-avt.dose<-list()
  
  for(dd in 1:length(dcount.l))
  {
    for(new.strat in 2:8)
    {
      #new.strat<-
      cov.sens.loader(strain=2, end.cov=cov.f[1], version='v1')
      cov.sens.loader(strain=3, end.cov=cov.f[1], version='v1')
      cov.sens.loader(strain=1, end.cov=cov.f[1], version='v1')
      
      uk.CEA.allstrains(Dataset1H1N1 = uk.tableH1N1[[1]], Dataset2H1N1 = uk.tableH1N1[[new.strat]],
                        Dataset1H3N2 = uk.tableH3N2[[1]], Dataset2H3N2 = uk.tableH3N2[[new.strat]],
                        Dataset1B = uk.tableB[[1]], Dataset2B = uk.tableB[[new.strat]],
                        dcount = dd, program1=1, program2=new.strat, i.cov=1, new.cov=mL3, version='v1')
    }
  }




####################################################################################################################



#n.samples<-100
for(strainpull in 1:length(strain.name))
{
  for(strategy in 1:4)
  {
    if(strategy==2) next
    singleSENS.samp.uk(strain=strainpull, i.country = 4, program=strategy,num.samp = n.samples, short = T, new.cov=seventy.cov, version='u7');
  }
}


#########################################################################################################
#Compare to Edwin
##########################################################################################################
n.samples<-2000
vss<-'v1'
  uk.tableB<-uk.tableH1N1<-uk.tableH3N2<-NULL
  pp<-ss<-list()
  
cov.sens.loader(2,cov.f[3], version = vss)
cov.sens.loader(1,cov.f[3], version = vss)
cov.sens.loader(3,cov.f[3], version = vss)

uk.prime<-Map('+', Map('+', uk.tableH3N2[[3]], uk.tableH1N1[[3]]), uk.tableB[[3]])
uk.sq<-Map('+', Map('+', uk.tableH3N2[[1]], uk.tableH1N1[[1]]), uk.tableB[[1]])
uk.sec<-Map('+', Map('+', uk.tableH3N2[[4]], uk.tableH1N1[[4]]), uk.tableB[[4]])


#uk.prime<-Map('+', uk.tableH3N2[[3]], uk.tableH1N1[[3]])
#uk.sq<-Map('+', uk.tableH3N2[[1]], uk.tableH1N1[[1]])
#uk.sec<-Map('+', uk.tableH3N2[[4]], uk.tableH1N1[[4]])

###############
#unlist(rowSums(uk.sq[[19]]-uk.prime[[19]])/sum(popv[1:22]))/unlist(rowSums(uk.sq[[19]])/sum(popv[1:22]))

pp<-Map('-',1, Map('/',lapply(Map('-',uk.sq,uk.prime),rowSums), lapply(uk.sq, rowSums)))
ss<-Map('-', 1,Map('/',lapply(Map('-',uk.sq,uk.sec),rowSums), lapply(uk.sq, rowSums)))


pp<-Map('/',lapply(Map('-',uk.sq,uk.prime),rowSums), lapply(uk.sq, rowSums))
ss<-Map('/',lapply(Map('-',uk.sq,uk.sec),rowSums), lapply(uk.sq, rowSums))


ppH3<-Map('/',lapply(Map('-',uk.tableH3N2[[1]], uk.tableH3N2[[3]]), rowSums),  lapply(uk.tableH3N2[[1]], rowSums))
ppH1<-Map('/',lapply(Map('-',uk.tableH1N1[[1]], uk.tableH1N1[[3]]), rowSums),  lapply(uk.tableH1N1[[1]], rowSums))
ppB<-Map('/',lapply(Map('-',uk.tableB[[1]], uk.tableB[[3]]), rowSums),  lapply(uk.tableB[[1]], rowSums))


#ppH3<-Map('/',lapply(Map('-',uk.tableH3N2[[3]], uk.tableH3N2[[1]]), rowSums),  lapply(uk.tableH3N2[[1]], rowSums))
#ppH1<-Map('/',lapply(Map('-',uk.tableH1N1[[3]], uk.tableH1N1[[1]]), rowSums),  lapply(uk.tableH1N1[[1]], rowSums))
#ppB<-Map('/',lapply(Map('-',uk.tableB[[3]], uk.tableB[[1]]), rowSums),  lapply(uk.tableB[[1]], rowSums))



primary<- Map('/', Map('+', Map('+',ppB, ppH1), ppH3), 3)


#ss<-Map('/',lapply(Map('-',uk.sq,uk.sec),rowSums), lapply(uk.sq, rowSums))


ssH3<-Map('/',lapply(Map('-',uk.tableH3N2[[1]], uk.tableH3N2[[4]]), rowSums),  lapply(uk.tableH3N2[[1]], rowSums))
ssH1<-Map('/',lapply(Map('-',uk.tableH1N1[[1]], uk.tableH1N1[[4]]), rowSums),  lapply(uk.tableH1N1[[1]], rowSums))
ssB<-Map('/',lapply(Map('-',uk.tableB[[1]], uk.tableB[[4]]), rowSums),  lapply(uk.tableB[[1]], rowSums))

#ssH3<-Map('/',lapply(Map('-',uk.tableH3N2[[4]], uk.tableH3N2[[1]]), rowSums),  lapply(uk.tableH3N2[[1]], rowSums))
#ssH1<-Map('/',lapply(Map('-',uk.tableH1N1[[4]], uk.tableH1N1[[1]]), rowSums),  lapply(uk.tableH1N1[[1]], rowSums))
#ssB<-Map('/',lapply(Map('-',uk.tableB[[4]], uk.tableB[[1]]), rowSums),  lapply(uk.tableB[[1]], rowSums))

secondary<-Map('/', Map('+', Map('+',ssB, ssH1), ssH3), 3)

#pp<-Map('/',lapply(Map('-',uk.sq,uk.prime),rowSums), lapply(uk.sq, rowSums))
#ss<-Map('/',lapply(Map('-',uk.sq,uk.sec),rowSums), lapply(uk.sq, rowSums))
#assign('p1',pp, envir=.GlobalEnv)
#assign('s1',ss, envir=.GlobalEnv)
#return(list(pp,ss))


for(ii in 14:length(primary)) {
  m.out.prime<-data.frame(cbind(c(ppH1[[ii]], ppH3[[ii]], ppB[[ii]]), rep('primary',n.samples*3), rep(paste0(y.name[[ii]],'/', y.name[ii+1]))))
  m.out.sec<-data.frame(cbind(c(ssH1[[ii]], ssH3[[ii]], ssB[[ii]]), rep('secondary',n.samples*3), rep(paste0(y.name[[ii]],'/', y.name[ii+1]))))
  
  m.prime<-rbind(m.prime, m.out.prime)
  m.sec<-rbind(m.sec, m.out.sec)
  
}



#1-(lapply(uk.prime,rowSums)[[19]]/sum(popv[1:22]))/(lapply(uk.sq,rowSums)[[19]]/sum(popv[1:22]))

#primary<-pp
#secondary<-ss

for(hhj in 1:19)
{
primary[[hhj]]<-c(unlist(tt1[[hhj]]), unlist(tt2[[hhj]]), unlist(tt3[[hhj]]))
secondary[[hhj]]<-c(unlist(yy1[[hhj]]), unlist(yy2[[hhj]]), unlist(yy3[[hhj]]))
}

n.samples<-2000

m.prime<-c()
m.sec<-c()

for(ii in 14:length(primary)) {
  m.out.prime<-data.frame(cbind(unlist(primary[[ii]]), rep('primary',n.samples), rep(paste0(y.name[[ii]],'/', y.name[ii+1]))))
  m.out.sec<-data.frame(cbind(unlist(secondary[[ii]]), rep('secondary',n.samples), rep(paste0(y.name[[ii]],'/', y.name[ii+1]))))
  
  m.prime<-rbind(m.prime, m.out.prime)
  m.sec<-rbind(m.sec, m.out.sec)
  
}

new.set<-rbind(m.prime, m.sec)
ggplot(new.set, aes(x=X3, y=as.numeric(as.character(X1)), fill=as.factor(X2))) +scale_y_continuous(limits = c(0,1))+
  geom_violin(scale='width',position=position_dodge(1))+scale_fill_manual(values=c("#E69F00", "#56B4E9"))+theme_bw()

#unlist(rowSums(uk.sq[[19]]-uk.prime[[19]])/sum(popv[1:22]))/unlist(rowSums(uk.sq[[19]])/sum(popv[1:22]))

xx1<-Map('/',Map('/',lapply(Map('-',uk.tableH1N1[[1]],uk.tableH1N1[[3]]),rowSums), sum(popv[1:22])), Map('/',lapply(uk.tableH1N1[[1]], rowSums), sum(popv[1:22])))



tt1<-Map('/', Map('/', lapply(uk.tableH1N1[[3]],rowSums), sum(popv[1:22])), Map('/',lapply(uk.tableH1N1[[1]], rowSums), sum(popv[1:22])))
tt2<-Map('/', Map('/', lapply(uk.tableH3N2[[3]],rowSums), sum(popv[1:22])), Map('/',lapply(uk.tableH3N2[[1]], rowSums), sum(popv[1:22])))
tt3<-Map('/', Map('/', lapply(uk.tableB[[3]],rowSums), sum(popv[1:22])), Map('/',lapply(uk.tableB[[1]], rowSums), sum(popv[1:22])))

primary<- Map('/', Map('+', Map('+',tt1, tt2), tt3), 3)



yy1<-Map('/', Map('/', lapply(uk.tableH1N1[[4]],rowSums), sum(popv[1:22])), Map('/',lapply(uk.tableH1N1[[1]], rowSums), sum(popv[1:22])))
yy2<-Map('/', Map('/', lapply(uk.tableH3N2[[4]],rowSums), sum(popv[1:22])), Map('/',lapply(uk.tableH3N2[[1]], rowSums), sum(popv[1:22])))
yy3<-Map('/', Map('/', lapply(uk.tableB[[4]],rowSums), sum(popv[1:22])), Map('/',lapply(uk.tableB[[1]], rowSums), sum(popv[1:22])))

secondary<- Map('/', Map('+', Map('+',yy1, yy2), yy3), 3)


#primary<-xx3
xx2<-Map('/',Map('/',lapply(uk.tableH3N2[[1]],uk.tableH3N2[[3]]),rowSums), sum(popv[1:22])), Map('/',lapply(uk.tableH3N2[[1]], rowSums), sum(popv[1:22])))
xx3<-Map('/',Map('/',lapply(Map('-',uk.tableB[[1]],uk.tableB[[3]]),rowSums), sum(popv[1:22])), Map('/',lapply(uk.tableB[[1]], rowSums), sum(popv[1:22])))

primary<- Map('/', Map('+', Map('+',xx1, xx2), xx3), 3)

yy1<-Map('/',Map('/',lapply(Map('-',uk.tableH1N1[[1]],uk.tableH1N1[[4]]),rowSums), sum(popv[1:22])), Map('/',lapply(uk.tableH1N1[[1]], rowSums), sum(popv[1:22])))
#secondary<-yy3
yy2<-Map('/',Map('/',lapply(Map('-',uk.tableH3N2[[1]],uk.tableH3N2[[4]]),rowSums), sum(popv[1:22])), Map('/',lapply(uk.tableH3N2[[1]], rowSums), sum(popv[1:22])))
yy3<-Map('/',Map('/',lapply(Map('-',uk.tableB[[1]],uk.tableB[[4]]),rowSums), sum(popv[1:22])), Map('/',lapply(uk.tableB[[1]], rowSums), sum(popv[1:22])))

secondary<- Map('/', Map('+', Map('+',yy1, yy2), yy3), 3)



lapply(Map('/',lapply(Map('-',uk.tableH1N1[[1]],uk.tableH1N1[[4]]),rowSums), lapply(uk.tableH1N1[[1]], rowSums)), median)
lapply(Map('/',lapply(Map('-',uk.tableH3N2[[1]],uk.tableH3N2[[4]]),rowSums), lapply(uk.tableH3N2[[1]], rowSums)), median)


xx<-Map('/',lapply(uk.tableH1N1[[1]], rowSums), sum(popv[1:22]))
yy<-Map('/',lapply(uk.tableH1N1[[3]], rowSums), sum(popv[1:22]))

primary<-Map('/',Map('-', xx, yy),xx) 


xx<-Map('/',lapply(uk.tableH1N1[[1]], rowSums), sum(popv[1:22]))
yy<-Map('/',lapply(uk.tableH1N1[[4]], rowSums), sum(popv[1:22]))

primary<-Map('/',Map('-', xx, yy),xx) 

yy1<-Map('/',lapply(Map('-',uk.tableH1N1[[1]],uk.tableH1N1[[3]]),rowSums),lapply(uk.tableH1N1[[1]], rowSums))
primary<-yy1

yy2<-Map('/',lapply(Map('-',uk.tableH1N1[[1]],uk.tableH1N1[[4]]),rowSums),lapply(uk.tableH1N1[[1]], rowSums))
secondary<-yy2

yy1<-Map('/',lapply(Map('-',uk.tableH3N2[[1]],uk.tableH3N2[[3]]),rowSums),lapply(uk.tableH3N2[[1]], rowSums))
primary<-yy1
yy2<-Map('/',lapply(Map('-',uk.tableH3N2[[1]],uk.tableH3N2[[4]]),rowSums),lapply(uk.tableH3N2[[1]], rowSums))
secondary<-yy2

yy3<-Map('/',lapply(Map('-',uk.tableB[[1]],uk.tableB[[3]]),rowSums),lapply(uk.tableB[[1]], rowSums))




mark<-1
library(coda)
mcmc.diagnotic<-function(strain.choice,vern, mark,pull)
{
 strain.cap<-c('H1','H3','B')
 if(strain.choice==3){flu.tables<-list.files(pattern=glob2rx(paste0('GBflu.', vern, '*B')));}
 if(strain.choice==1){flu.tables<-list.files(pattern=glob2rx(paste0('GBflu.', vern,'*H1')));}
 if(strain.choice==2){flu.tables<-list.files(pattern=glob2rx(paste0('GBflu.',vern, '*H3')));}

 setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/GB")
 year.length<-length(flu.tables)

  gelman.output<-list()
  last.pull<<-list()
  for(hh in 1:year.length)
  {
    
    load(flu.tables[[hh]]); chain1<-as.mcmc(mcmc.result$batch[(mark:dim(mcmc.result$batch)[1]),])
    last.pull[[hh]]<<- mcmc.result$batch[(dim(mcmc.result$batch)[1]),]
    #load(flu.tablesB); chain3<- as.mcmc(mcmc.result$batch[(mark:dim(post.sample)[1]),])
    
    if(pull==F)
    {
      load(flu.tables.b[[hh]]); chain2<- as.mcmc(mcmc.result$batch[(mark:dim(mcmc.result$batch)[1]),])
      xx<-mcmc.list(chain1, chain2)
      f<-gelman.diag(xx, confidence = 0.95, transform=T, autoburnin=F, multivariate = T)
      
      gelman.output[[hh]]<- f
      gelman.plot(xx)
    }
    
  }
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata','GB'));
  save(last.pull, file=paste0('GB',strain.cap[strain.choice],vern, 'last.pull'))
}

mcmc.diagnotic(1,vern='u9',50000, pull=T)
mcmc.diagnotic(2,vern='u9',50000, pull=T)
mcmc.diagnotic(3,vern='u9',50000, pull=T)

load('GBH1N1last.pull')


######################################################################################################################
#CEP Plane Graphs
#########################################################################################################

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 8
cols <- gg_color_hue(n)
cols<-c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC")

inv.names<-c("V2-4y"  ,                 "V5-11y"       ,              "V12-16y"      ,             "V2-11y" ,"V2-16y", "V2-4y/V12-16y"   ,      "V5-16y" )
#inv.names<-c("Preschool","Primary","Secondary","Preschool+Primary School","Preschool+Primary+Secondary","Preschool+Secondary","Primary+Secondary")
library(scales)

cep.graph2 <-function(cov.var, dcount){  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB', 'coverage',cov.f[cov.var]))
  cea.tables<-list.files(pattern=glob2rx(paste0('GBICER.all.strains', cov.f[cov.var], 'p', '*', 'd',dcount)))
  
  cos.f<-qds.f<-qds.m<-cos.m<-cep<-cep1<-contour_95<-mean.ff<-mean.po<-not.ce<-NULL
  
  for(ff in 1:length(inv.names))
  {#ff=interventions. Strains are already chosen
    contour_95<-NULL
    load(cea.tables[[ff]]) #H1N1
    
    cer<-ICER.set2 #grab ICERs
    
    qa<-Reduce('+', cer$`net QALY`)/tab.diff
    co<-Reduce('+', cer$`net cost`)/tab.diff
    
    #m.qa<-rowSums(qa)
    #m.co<-rowSums(co)
    
    
    qds<-data.frame(cbind(rep(fname[i.country], length(qa)), rep(inv.names[[ff]], length(qa)), rep(dcount.l[dcount], length(qa)), qa))
    cos<-data.frame(cbind(rep(fname[i.country], length(qa)) ,rep(inv.names[[ff]], length(co)), rep(dcount.l[dcount], length(co)), co))
    
    colnames(qds)<-c('country','Intervention','discount','total')
    
    colnames(qds)<-colnames(cos)
    
    ifelse(ff==1, qds.f<-qds, qds.f<-rbind(qds.f, qds))
    ifelse(ff==1, cos.f<-cos, cos.f<-rbind(cos.f, cos))
    
    xy.points<-data.frame(qa, co)
    names(xy.points)<-c('x', 'y')
    kd <- ks::kde(xy.points, compute.cont=TRUE)
    contour_95 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                        z=estimate, levels=cont["10%"])[[1]])
    contour_95 <- data.frame(contour_95, rep(inv.names[[ff]], length(contour_95$x)))
    colnames(contour_95)<-c('level','x','y','Intervention')
    
    p1<-geom_path(data=contour_95, aes(x, y, color=Intervention))
    contour.name<-paste0('cont.all.strains', cov.var, dcount.l[dcount], ff)
    assign(contour.name,  p1)
  }
  
  colnames(qds.f)<-c('country','Intervention','discount','total')
  qds.i<-melt(qds.f, id.vars = c('country', 'Intervention','discount'))
  colnames(cos.f)<-colnames(qds.f)
  cos.i<-melt(cos.f, id.vars = c('country', 'Intervention','discount'))
  #setwd('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/IntlSens')
  
  llg<-cbind(qds.i[,1:4],as.numeric(as.character(qds.i$value)), as.numeric(as.character(cos.i$value))) 
  colnames(llg)<-c('country','Intervention','discount', 'Age','qaly','cost')
  llg<-na.omit(llg)
  
  totals<-llg[llg$Age=='total',]
  #library('plyr')
  mean.p<-ddply(totals, .(Intervention), summarize, Rate1=median(qaly), Rate2=median(cost))
  mean.po<-mean.p[order(mean.p$Rate2, decreasing = T),]
  ##   continent n_uniq_countries
  
  
  
  ######CE evaluation
  not.ce<-not.ce1<-not.ce2<-not.ce3<-not.ce4<-counter<-NULL
  
  rudi.icer<-diff(rbind(c(0,0),cbind(rev(mean.po[,2]), rev(mean.po[,3]))))
  mean.f<-mean.po[which((rudi.icer[,2]/rudi.icer[,1]) >0&(rudi.icer[,2]/rudi.icer[,1] <20000)),]
  levels(mean.f$Intervention)= c("I1: 2-4 yrs old"  ,                 "I2: 5-11 yrs old"       ,              "I3: 12-16 yrs old"      ,             "I4: 2-11 yrs old" ,"S7: 2-16 yrs old", "I5: 2-4 & 12-16 yrs old"   ,      "I6: 5-16 yrs old" , 'SQ: Risk Groups, 65+ yrs old')
  mean.ff<-data.frame(rbind(c('SQ: Risk Groups, 65+ yrs old',as.numeric(0),as.numeric(0)), mean.f[order(mean.f$Rate2, decreasing=F),]))
  
  not.ce<-mean.po[which(rudi.icer[,2]/rudi.icer[,1] <0 | rudi.icer[,2]/rudi.icer[,1] >20000),]
  
  if(length(not.ce)!=0)
  {
    counter<-2
    ##second icer eval
    rudi.icer2<-diff(rbind(c(0,0),cbind(rev(mean.f[,2]), rev(mean.f[,3]))))
    mean.f2<-mean.f[which(rudi.icer2[,2]/rudi.icer2[,1] >0&rudi.icer2[,2]/rudi.icer2[,1] <20000),]
    levels(mean.f2$Intervention)= c("I1: 2-4 yrs old"  ,                 "I2: 5-11 yrs old"       ,              "I3: 12-16 yrs old"      ,             "I4: 2-11 yrs old" ,"I7: 2-16 yrs old", "I5: 2-4 & 12-16 yrs old"   ,      "I6: 5-16 yrs old" , 'SQ: Risk Groups, 65+ yrs old')
    mean.ff2<-data.frame(rbind(c('SQ: Risk Groups, 65+ yrs old',as.numeric(0),as.numeric(0)), mean.f2[order(mean.f2$Rate2, decreasing=F),]))
    
    not.ce2<-rbind(not.ce, mean.f[which(rudi.icer2[,2]/rudi.icer2[,1] <0 | rudi.icer2[,2]/rudi.icer2[,1] >20000),])
  }
  
  if(any(dim(not.ce2)!=dim(not.ce)))
  {
    counter<-counter+1
    ##3rd CE check
    #mean.po<-mean.p[order(mean.p$Rate2, decreasing = T),]
    rudi.icer3<-diff(rbind(c(0,0),cbind(rev(mean.f2[,2]), rev(mean.f2[,3]))))
    mean.f3<-mean.f2[which(rudi.icer3[,2]/rudi.icer3[,1] >0&rudi.icer3[,2]/rudi.icer3[,1] <20000),]
    levels(mean.f3$Intervention)= c("I1: 2-4 yrs old"  ,                 "I2: 5-11 yrs old"       ,              "I3: 12-16 yrs old"      ,             "I4: 2-11 yrs old" ,"I7: 2-16 yrs old", "I5: 2-4 & 12-16 yrs old"   ,      "I6: 5-16 yrs old" , 'SQ: Risk Groups, 65+ yrs old')
    mean.ff3<-data.frame(rbind(c('SQ: Risk Groups, 65+ yrs old',as.numeric(0),as.numeric(0)), mean.f3[order(mean.f3$Rate2, decreasing=F),]))
    
    not.ce3<-rbind(not.ce2, mean.f[which(rudi.icer3[,2]/rudi.icer3[,1] <0 | rudi.icer3[,2]/rudi.icer3[,1] >20000),])
  }
  
  if(any(dim(not.ce3)!=dim(not.ce2)))
  {
    counter<-counter+1
    ##4th CE check
    #mean.po<-mean.p[order(mean.p$Rate2, decreasing = T),]
    rudi.icer4<-diff(rbind(c(0,0),cbind(rev(mean.f3[,2]), rev(mean.f3[,3]))))
    mean.f4<-mean.f3[which(rudi.icer4[,2]/rudi.icer4[,1] >0&rudi.icer4[,2]/rudi.icer4[,1] <20000),]
    levels(mean.f4$Intervention)= c("I1: 2-4 yrs old"  ,                 "I2: 5-11 yrs old"       ,              "I3: 12-16 yrs old"      ,             "I4: 2-11 yrs old" ,"I7: 2-16 yrs old", "I5: 2-4 & 12-16 yrs old"   ,      "I6: 5-16 yrs old" , 'SQ: Risk Groups, 65+ yrs old')
    mean.ff4<-data.frame(rbind(c('SQ: Risk Groups, 65+ yrs old',as.numeric(0),as.numeric(0)), mean.f4[order(mean.f4$Rate2, decreasing=F),]))
    
    not.ce4<-rbind(not.ce3, mean.f[which(rudi.icer4[,2]/rudi.icer4[,1] <0 | rudi.icer4[,2]/rudi.icer4[,1] >20000),])
  }
  
  if(counter>1)
  {
    mean.ff<-get(paste0('mean.ff', counter))
    not.ce<-get(paste0('not.ce', counter))
  }
  
  
  ####
  cep<-cep1<-NULL
  #llg[llg$cost>1e6|llg$cost < (-1e6),]<-NA
  #cep<-ggplot(totals, aes(qaly, cost, color=Intervention))+
   # scale_y_continuous(breaks=scales::pretty_breaks(n = 6))+scale_x_continuous(breaks=scales::pretty_breaks(n=6))+geom_hline(aes(yintercept=0))+geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+geom_abline(slope = 15000, intercept=0, linetype='dotdash')+geom_abline(slope = 20000, intercept=0, linetype='longdash')+labs(subtitle=paste(paste0('discount=',dcount.l[dcount],'%'), paste0('strain=all'), paste0('coverage=',cov.f[cov.var])), fill='Intervention')
  
  cep<-ggplot(totals, aes(qaly, cost, color=Intervention))+
    scale_y_continuous(breaks=scales::pretty_breaks(n = 6))+scale_x_continuous(breaks=scales::pretty_breaks(n=6))+geom_hline(aes(yintercept=0))+geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+geom_abline(slope = 15000, intercept=0, linetype='dotdash')+geom_abline(slope = 20000, intercept=0, linetype='longdash')+labs(fill='Intervention')
  
  
  #+ylab('Incremental Cost in million ₤/year')+xlab('QALY Difference')+labs(subtitle=paste(paste0('discount=',dcount.l[dcount],'%'), paste0('strain=',strain.name[strain]), paste0('coverage=',cov.f[cov.var])), fill='Intervention')
  
  for(fff in 1:length(inv.names))
  {add.on<-get(paste0('cont.all.strains', cov.var, dcount.l[dcount], fff))
  cep<-cep+add.on}
  #cep.name<-paste0('cep',strain, cov.var, dcount.l[dcount])
  
  #assign(cep.name, 
  
  cep1<-cep+geom_point(data=mean.ff, aes(x=as.numeric(mean.ff$Rate1), y=as.numeric(mean.ff$Rate2)),size=3, shape=19, stroke=1, inherit.aes = T, show.legend=F)+geom_line(data=mean.ff, aes(x=as.numeric(mean.ff$Rate1), y=as.numeric(mean.ff$Rate2)), inherit.aes = F)+geom_point(data=not.ce, aes(x=not.ce$Rate1, not.ce$Rate2), size=3, shape=10, inherit.aes = T, show.legend=F)+labs(fill='Intervention', x='QALY Differences', y="Incremental Cost in £/year")+theme_bw()
  #)
  
  #assign(cep.name, cep1)
  #cep<-cep1<-NULL
  #setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB', 'coverage', cov.f[cov.var]))
  #ggsave(paste0('CEP',strain.name[strain],cov.f[cov.var],dcount,".png"),plot = cep, width=11, height=8.5, units='in',device='png')
  return(cep1)
  
}


cep.grapher.all<-function(cov.var, dcount, wtp)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB', 'coverage',cov.f[cov.var]))
  cea.tables<-list.files(pattern=glob2rx(paste0('GBICER.all.strains', cov.f[cov.var], 'p', '*', 'd',dcount)))
  
  cos.f<-qds.f<-qds.m<-cos.m<-cep<-cep1<-contour_95<-mean.ff<-mean.po<-not.ce<-NULL
  
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
    afford<-sum(unlist(ICER.set2$icer)<10000)/(tab.diff*length(ICER.set2$icer[[1]]))
    ifelse(ff==1, afford.out<<-c(inv.names[ff], afford), afford.out<<-rbind(afford.out, c(inv.names[ff], afford)))
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
    contour.name<-paste0('cont',cov.var,dcount.l[dcount], ff)
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
  
  
  ######CE evaluation
  #not.ce<-not.ce1<-not.ce2<-not.ce3<-not.ce4<-counter<-NULL
  #counter<-1
  #rudi.icer<-diff(rbind(rep(0,6),cbind(rev(mean.po[,3]), rev(mean.po[,4]))))
  #mean.f<-mean.po[which(rudi.icer[,2]/rudi.icer[,1] >0&rudi.icer[,2]/rudi.icer[,1] <wtp),]
  #levels(mean.f$Intervention)= c("Preschool"  ,"Primary"  ,  "Secondary", 'Status Quo')
  #mean.ff1<-data.frame(rbind(c(fname[i.country],'Status Quo',as.numeric(0), as.numeric(0)), mean.f[order(as.numeric(as.character(mean.f$Rate2, decreasing=F))),]))
  
  #not.ce1<-mean.po[which(rudi.icer[,2]/rudi.icer[,1] <0 | rudi.icer[,2]/rudi.icer[,1] >wtp),]
  
  #if(length(not.ce)!=0 & dim(mean.f)[1]!=0)
  #{
  #counter<-counter+1
  ##second icer eval
  #rudi.icer2<-diff(rbind(c(0,0),cbind(rev(mean.f[,3]), rev(mean.f[,4]))))
  #mean.f2<-mean.f[which(rudi.icer2[,2]/rudi.icer2[,1] >0&rudi.icer2[,2]/rudi.icer2[,1] <wtp),]
  #levels(mean.f2$Intervention)= c("Preschool"  , "Primary" , "Secondary"  ,  'Status Quo')
  #mean.ff2<-data.frame(rbind(c(fname[i.country],'Status Quo',as.numeric(0),as.numeric(0)), mean.f2[order(mean.f2$Rate2, decreasing=F),]))
  
  # not.ce2<-rbind(not.ce, mean.f[which(rudi.icer2[,2]/rudi.icer2[,1] <0 | rudi.icer2[,2]/rudi.icer2[,1] >wtp),])
  
  
  #if(any(dim(not.ce2)!=dim(not.ce))&length(mean.ff2)!=4&dim(mean.f2)[1]!=0)
  #{
  # counter<-counter+1
  ##3rd CE check
  #mean.po<-mean.p[order(mean.p$Rate2, decreasing = T),]
  #rudi.icer3<-diff(rbind(c(0,0),cbind(rev(mean.f2[,3]), rev(mean.f2[,4]))))
  #mean.f3<-mean.f2[which(rudi.icer3[,2]/rudi.icer3[,1] >0&rudi.icer3[,2]/rudi.icer3[,1] <wtp),]
  #levels(mean.f3$Intervention)= c("Preschool"  ,"Primary",  "Secondary"  , 'Status Quo')
  #mean.ff3<-data.frame(rbind(c(fname[i.country], 'Status Quo',as.numeric(0),as.numeric(0)), mean.f3[order(mean.f3$Rate2, decreasing=F),]))
  
  #not.ce3<-rbind(not.ce2, mean.f[which(rudi.icer3[,2]/rudi.icer3[,1] <0 | rudi.icer3[,2]/rudi.icer3[,1] >wtp),])
  #}
  #}
  
  
  #mean.ff<-get(paste0('mean.ff', counter))
  #not.ce<-get(paste0('not.ce', counter))
  #colnames(mean.ff)<-c('country','Intervention', 'QALYS', 'Costs')
  #colnames(not.ce)<-c('country','Intervention', 'QALYS', 'Costs')
  
  
  
  ####
  cep<-cep1<-NULL
  cols<-c("#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC")
  
  #llg[llg$cost>1e6|llg$cost < (-1e6),]<-NA
  cep<-ggplot(cost.utility, aes(as.numeric(as.character(QALYS)), as.numeric(as.character(Costs)), color=Intervention))+scale_y_continuous(breaks=scales::pretty_breaks(n = 7), limits=c(0,1e8))+scale_x_continuous(breaks=scales::pretty_breaks(n=6), limits=c(-1000,71000))+geom_hline(aes(yintercept=0))+geom_hline(aes(yintercept=0))+geom_vline(aes(xintercept=0))+geom_abline(slope = 15000, intercept=0, linetype='dotdash')+geom_abline(slope = wtp, intercept=0, linetype='longdash')+labs(subtitle=paste(paste0('discount=',dcount.l[dcount],'%'), paste0('strain=All'), paste0('coverage=',cov.f[cov.var])), fill='Intervention')+theme_bw()+scale_colour_manual(values=cols)+theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12))
  
  
  #+ylab('Incremental Cost in million ₤/year')+xlab('QALY Difference')+labs(subtitle=paste(paste0('discount=',dcount.l[dcount],'%'), paste0('strain=',strain.name[strain]), paste0('coverage=',cov.f[cov.var])), fill='Intervention')
  
  for(fff in 1:length(inv.names))
  {add.on<-get(paste0('cont',cov.var,dcount.l[dcount], fff))
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
      cep1<-cep+geom_point(data=mean.ff, aes(x=as.numeric(as.character(Rate1)), y=as.numeric(as.character(Rate2))),size=3, shape=19, stroke=1, inherit.aes = T, show.legend=F)+geom_path(data=mean.ff, aes(x=as.numeric(Rate1), y=as.numeric(Rate2)), inherit.aes = F)+geom_point(data=not.ce, aes(x=as.numeric(as.character(not.ce$Rate1)), y=as.numeric(as.character(not.ce$Rate2))), size=3, shape=10, inherit.aes = T, show.legend=F)+labs(fill='Intervention')+theme_bw(base_size = 14)+geom_line(data=as.data.frame(points), aes(x=V1, y=V2), inherit.aes = F)
    }else{
      colnames(mean.ff)<-c('Intervention','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB','cost.sd','qal.sd', 'icer', 'icerUB', 'icerLB', 'icersdLB', 'icersdUB')
      cep1<-cep+geom_point(data=mean.ff, aes(x=as.numeric(as.character(Rate1)), y=as.numeric(as.character(Rate2))),size=3, shape=19, stroke=1, inherit.aes = T, show.legend=F)+theme_bw(base_size = 14)+geom_line(data=as.data.frame(points), aes(x=V1, y=V2), inherit.aes = F)
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
tab.diff<-19
cep0.55<-cep.grapher.all(cov.var = 2,dcount = 3,wtp = 20000)+rremove('xy.title')
cep0.552<-annotate_figure(cep0.55,
                top = text_grob("Cost-Effectiveness Plane"),
                bottom = text_grob("QALY Differences"),
                left = text_grob("Incremental Cost (£GBP) /year", rot=90))
ggsave(paste0('CEPallstrainsfig',cov.f[2],'3.5', version,".png"),plot = cep0.552, width=9, height=5.75, units='in',device='png')

 

wtpl<-20000
####ggarrange for all strains, discount=3.5%, different coverages
all32<-ggarrange(cep.grapher.all(1,3, wtpl)+rremove("xy.title"), cep.grapher.all(2,3, wtpl)+rremove("xy.title"), cep.grapher.all(3,3, wtpl)+rremove("xy.title")+rremove("xy.title"), common.legend=T, legend='bottom',nrow=1, ncol=3)

all32<-annotate_figure(all32,
                       top = text_grob("Cost-Effectiveness Plane"),
                       bottom = text_grob("QALY Differences"),
                       left = text_grob("Incremental Cost (£GBP) /year", rot=90)
)

ggsave(paste0('CEPdiffcoverage','3.5',version,".png"),plot = all32, width=13, height=5.75, units='in',device='png')







###########################CEP for all strains######################################################
discount3<-ggarrange(
    cep.grapher.all(1,1, wtpl)+rremove("xy.title"), 
    cep.grapher.all(2,1, wtpl)+rremove("xy.title"),
    cep.grapher.all(3,1, wtpl)+rremove("xy.title"),
    cep.grapher.all(1,2, wtpl)+rremove("xy.title"), 
    cep.grapher.all(2,2, wtpl)+rremove("xy.title"),
    cep.grapher.all(3,2, wtpl)+rremove("xy.title"),
    cep.grapher.all(1,3, wtpl)+rremove("xy.title"), 
    cep.grapher.all(2,3, wtpl)+rremove("xy.title"),
    cep.grapher.all(3,3, wtpl)+rremove("xy.title"),
    ncol=3, nrow=3, common.legend = T, legend='bottom')

discount32<-annotate_figure(discount3,
                top = text_grob("Cost-Effectiveness Plane (All Strains)"),
                bottom = text_grob("QALY Differences"),
                left = text_grob("Incremental Cost (GBP) /year", rot=90)
)

ggsave(paste0("CEPallstrainvcoverage",version,".png"),plot=discount32, width=11, height=8.5, units='in',device='png')



######################################################coverage 0.55, discount 3.5, all strains
all33<-ggarrange(cep.graph2(2,3)+rremove("xy.title"), common.legend=T, legend='bottom',nrow=1, ncol=1)

all33<-annotate_figure(all33,
                       top = text_grob("Cost-Effectiveness Plane"),
                       bottom = text_grob("QALY Differences"),
                       left = text_grob("Incremental Cost in ₤/year", rot=90)
)
ggsave(paste0("CEP05535",version,".pdf"),plot=all33, width=7, height=9, units='in',device='pdf')



##############################################################################################
#########Graph barplot for economic outcomes by strain and discount
##############################################################################################

###loader
discount.loader2<-function(cov.var, dcount)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB','coverage', cov.f[cov.var]))
  
  #Rearrange tables to correct order of interventions------------------
  cost.tables<-list.files(pattern=glob2rx(paste0('GBCEAcosts.allstrain',cov.f[cov.var], 'p*', 'd', dcount)))
  
  load(cost.tables[1]); avg.incidence1<-coststab; #status quo
  load(cost.tables[2]); avg.incidence2<-coststab; #preschool
  load(cost.tables[3]); avg.incidence3<-coststab; #primary
  load(cost.tables[4]); avg.incidence4<-coststab; #secondary
  load(cost.tables[5]); avg.incidence5<-coststab; #secondary
  load(cost.tables[6]); avg.incidence6<-coststab; #secondary
  load(cost.tables[7]); avg.incidence7<-coststab; #secondary
  
  
  var.name<-paste0('inv.costs', dcount.l[dcount])
  #('WTP','Preschool','Primary','Secondary','Preschool+Primary','PrePrimeSeconday','Preschool+Secondary','Primary+Secondary')
  assign(var.name, do.call( rbind.data.frame, list(avg.incidence1,avg.incidence2,avg.incidence3, avg.incidence4, avg.incidence5,avg.incidence6, avg.incidence7)
  ), envir = .GlobalEnv)
  
}


###########################################################################################################
#Loader for outcomes
############################################################################################################
outcomes.loader2<-function(cov.var, dcount)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB','coverage', cov.f[cov.var]))
  
  #Rearrange tables to correct order of interventions------------------
  ost.tables<-list.files(pattern=glob2rx(paste0('GBCEAoutcomes.allstrain',cov.f[cov.var], 'p*', 'd', dcount)))
  
  load(ost.tables[1]); avg.incidence1<-avg.outcomes; #status quo
  load(ost.tables[2]); avg.incidence2<-avg.outcomes; #preschool
  load(ost.tables[3]); avg.incidence3<-avg.outcomes; #primary
  load(ost.tables[4]); avg.incidence4<-avg.outcomes; #secondary
  load(ost.tables[5]); avg.incidence5<-avg.outcomes; #secondary
  load(ost.tables[6]); avg.incidence6<-avg.outcomes; #secondary
  load(ost.tables[7]); avg.incidence7<-avg.outcomes; #secondary
  
  
  var.name<-paste0('inv.outcomes', dcount.l[dcount])
  #('WTP','Preschool','Primary','Secondary','Preschool+Primary','PrePrimeSeconday','Preschool+Secondary','Primary+Secondary')
  assign(var.name, do.call( rbind.data.frame, list(avg.incidence1,avg.incidence2,avg.incidence3, avg.incidence4, avg.incidence5,avg.incidence6, avg.incidence7)
  ), envir = .GlobalEnv)
  
}

outcomes.loader2(2, 3)
discount.loader2(2,3)
############################################################################################################
###Barplot graph for all strains
############################################################################################################

library(ggpubr)
library(gridExtra)

barplot2<-function(cov.var, dcount)
{
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB','coverage', cov.f[cov.var]))
  
  #make status quo category
  for(tt in 1:2)
  {
    discount.loader2(cov.var, dcount)
    inv.costs<-get(paste0('inv.costs',dcount.l[dcount]))
  
    
    inv.costs$variable <- factor(inv.costs$variable, levels=c(paste0('vaccine cost',tt), paste0('Healthcare costs',tt), paste0('GP costs',tt), paste0('hospital costs',tt),paste0('QALYs lost',tt),paste0('QALYs lost case',tt), paste0('QALYs lost death',tt), paste0('QALYs lost hosp', tt)), labels = c('Vaccination Costs','Healthcare Costs','GP Consult Fees','Hospitalization Costs','Years of Life Lost (All Sources)', 'Years of Life Lost (Symptomatic Cases)', 'Years of Life Lost (Deaths)','Years of Life Lost (Hospitalization)'))
    if(tt==1) inv.costs$strategy<-'Status Quo'
    ifelse(tt==1, inv.costs.f<-na.omit(inv.costs)[1:6,], inv.costs.f<-rbind(inv.costs.f, na.omit(inv.costs)))
  }
  #cols<-my.col
  
  inv.costs.f$strategy<-factor(inv.costs.f$strategy, levels=c('Preschool','Primary','Secondary','Preschool+Primary','PrePrimeSeconday', 'Preschool+Secondary','Primary+Secondary', 'Status Quo'), labels=c( "I1: 2-4 yrs old"  ,                 "I2: 5-11 yrs old"       ,              "I3: 12-16 yrs old"      ,             "I4: 2-11 yrs old" ,"I7: 2-16 yrs old", "I5: 2-4 & 12-16 yrs old"   ,      "I6: 5-16 yrs old", 'SQ: Risk Groups, 65+ yrs old'))
  
  q<-ggplot(inv.costs.f, aes(x=strategy, y=value, fill=strategy))+geom_bar(position = "dodge", stat = "identity", aes(fill=strategy))+facet_wrap(~variable, ncol=3, scales = 'free')+scale_y_continuous(breaks=scales::pretty_breaks(n=7))+labs(title = NULL,subtitle=paste('discount=', dcount.l[dcount], '%', ',','coverage=', cov.f[cov.var],',', 'strain=all'), x = 'Intervention Strategy', y = '"Average Annual Costs(£)"', fill = "Economic Impact\n") +geom_errorbar(aes(ymin=CIL, ymax=CIU), width=.2, position=position_dodge(.9))+theme_bw()+theme(
    axis.text.x = element_blank(), axis.ticks = element_blank(), plot.subtitle=element_text(size=10, hjust=0.5, color="black"), legend.position="bottom")
  
  
  ggsave(paste0('costsplot',cov.f[cov.var],dcount,".png"),plot = q, width=11, height=7.5, units='in',device='png')
  return(q)
}

##arrange by discount and coverage level
barplot3<-ggarrange(
  barplot2(1,1)+rremove("xy.title"), 
  barplot2(2,1)+rremove("xy.title"),
  barplot2(3,1)+rremove("xy.title"),
  barplot2(1,2)+rremove("xy.title"),
  barplot2(2,2)+rremove("xy.title"),
  barplot2(3,2)+rremove("xy.title"),
  barplot2(1,3)+rremove("xy.title"),
  barplot2(2,3)+rremove("xy.title"),
  barplot2(3,3)+rremove("xy.title"),
  ncol=1, nrow=1, common.legend = F, legend='bottom')

barplot32<-marrangeGrob(barplot3, nrow=1, ncol=1,
                            top = text_grob("Economic Outcomes (All Strains)"),
                            bottom = text_grob("Intervention Strategy"),
                            left = text_grob("Average Annual Costs(GBP£)", rot=90)
)

ggsave(paste0("barplotallstrain",version,".pdf"),plot=barplot32, width=11, height=8.5, units='in',device='pdf')

dev.off()


ggsave(paste0("barplot35055",version,".pdf"),plot=barplot2(2,3), width=11, height=8.5, units='in',device='pdf')



####################################################################################
#### Optimal Acceptability curve, probability each intervention has optimal cost-effectiveness
###########################################################################################

#inv.names2<-c('Status Quo',"Preschool"   ,                "Primary"      ,               "Secondary"  ,                 "Preschool+Primary School"   , "Preschool+Primary+Secondary" ,"Preschool+Secondary"   ,      "Primary+Secondary")


inv.names2<-c('SQ: Risk Groups, 65+ yrs old', "I1: 2-4 yrs old"  ,                 "I2: 5-11 yrs old"       ,              "I3: 12-16 yrs old"      ,             "I4: 2-11 yrs old" ,"I7: 2-16 yrs old", "I5: 2-4 & 12-16 yrs old"   ,      "I6: 5-16 yrs old")

inv.names3<-c('SQ: Risk Groups, 65+ yrs old', "I1: 2-4 yrs old"  ,                 "I2: 5-11 yrs old"       ,              "I3: 12-16 yrs old"      ,             "I4: 2-11 yrs old" , "I5: 2-4 & 12-16 yrs old"   ,      "I6: 5-16 yrs old", "I7: 2-16 yrs old")



setwd(paste0('/Users/Natasha/Dropbox/UKfluworkGIT/Outputdata/GB/coverage/', cov.f[cov.var]))


###calculate optimal outcomes from total qalys and costs
accept.curve<-function(threshold.val, qaly, cost)
{
  qaly2<- Map('*', -1, Map('*',qaly,threshold.val))
  nmb<-Map('-',qaly2, cost)
  
  return(nmb)
}

c.names<-c('program','0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR', 'total')


###overall function
optimal.strat.graph2<-function(cov.var, dcount)
{
  qal.all<-list()
  cost.all<-list()
  outerf<-list()
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage',cov.f[cov.var]))
  qal.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.qaly.lossall.strains',cov.f[cov.var], 'p*', 'd',dcount)))
  cos.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.cost.lossall.strains',cov.f[cov.var], '*', 'd',dcount)))
  
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
  
  
  ss<-seq(1000,30000,by = 500)
  for(hh in 1:length(ss)) {
    
    AC.all<-lapply(1:length(qal.all), function(oo) accept.curve(ss[hh], qal.all[[oo]], cost.all[[oo]]) )
    
    AC.all2<-lapply(1:length(qal.all), function(pp) data.frame(as.factor(rep(inv.names2,n.samples)), Reduce('+',AC.all[[pp]])/tab.diff, rowSums(Reduce('+',AC.all[[pp]])/tab.diff)))
    
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
  colnames(outerf)<-c('country', 'Age','threshold','SQ: Risk Groups, 65+ yrs old', "I1: 2-4 yrs old"  ,                 "I2: 5-11 yrs old"       ,              "I3: 12-16 yrs old"      ,             "I4: 2-11 yrs old" ,"I7: 2-16 yrs old", "I5: 2-4 & 12-16 yrs old"   ,      "I6: 5-16 yrs old")
  
    outer.g<-melt(outerf, id.vars=c('country','Age', 'threshold'))
  colnames(outer.g)<-c('country', 'Age', 'Threshold','Intervention','value')
  mk<-ggplot(data= outer.g[outer.g$Age=='total',], aes(x=as.numeric(as.character(Threshold)), y=as.numeric(as.character(value)), color=Intervention, linetype=Intervention))
  
  label<-paste('strain=all',paste0('discount=', dcount.l[dcount], '%'), paste0('coverage=', cov.f[cov.var]))
  
  #final graph
  OC3<-mk+labs(title=NULL, subtitle=label, fill='Intervention')+geom_line(aes(color=Intervention))+scale_y_continuous(breaks=scales::pretty_breaks(n=7), limits = c(0,1))+scale_x_continuous(breaks=scales::pretty_breaks(n=7))+geom_vline(xintercept = c(15000,20000), linetype=3)+theme(plot.subtitle=element_text(size=10, hjust=0.5, color="black"))+guides(fill=guide_legend(title="Intervention"))
  if(cov.var==2&dcount==3)
  {ggsave(paste0('OC.total.all.strains',cov.f[cov.var],'d',dcount.l[dcount],".png"),plot = OC3, width=11,height=8.5, units='in',device='png')};
  
  return(OC3)
}

cos.tables[[cov.var]]<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.cost.lossall.strains',cov.f[cov.var], '*', 'd','3')))[[2]]

###################################################
#QALY Gains
qal.tables<-cos.tables<-c()
setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage',cov.f[cov.var]))
qal.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.qaly.lossall.strains',cov.f[cov.var], 'p*', 'd',dcount)))

load(qal.tables[1])
Map('-',tot.QALY.loss[[1]],tot.QALY.loss[[2]])




##

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 8
cols <- gg_color_hue(n)
cols<-c("#FF61CC", "#F8766D","#CD9600","#7CAE00", "#00BE67","#C77CFF","#00BFC4","#00A9FF" )
#,'SQ: Risk Groups, 65+ yrs old', "S1: 2-4 yrs old"  ,                 "S2: 5-11 yrs old"       ,              "S3: 12-16 yrs old"      ,             "S4: 2-11 yrs old" ,"S7: 2-16 yrs old", "S5: 2-4 & 12-16 yrs old"   ,      "S6: 5-16 yrs old")
optimal.strat.new<-function(cov.var, dcount)
{
  qal.all<-list()
  cost.all<-list()
  outerf<-list()
  
   qal.tables<-cos.tables<-c()
    setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage',cov.f[cov.var]))
    qal.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.qaly.lossall.strains',cov.f[cov.var], 'p*', 'd',dcount)))
    cos.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.cost.lossall.strains',cov.f[cov.var], '*', 'd',dcount)))
  
  for(ii in 1:length(qal.tables)) 
  {
  #optimal strategy compares all strategy outcomes to each other and sums to 1
  #setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage',cov.f[ii]))
    load(qal.tables[[ii]])
    load(cos.tables[[ii]])
    
    if(ii==1){
      qal.all[[ii]]<-tot.QALY.loss[[1]];
      cost.all[[ii]]<-tot.costs.loss[[1]];
    }
    qal.all[[ii+1]]<-tot.QALY.loss[[2]]
    cost.all[[ii+1]]<-tot.costs.loss[[2]]
    
  }
  
  tab.diff<-length(qal.all[[1]])
  intv.number<-length(qal.all)
  qal.m<-lapply(lapply(1:length(qal.all), function(oo) Reduce('+', qal.all[[oo]])/tab.diff), rowSums)
  cost.m<-lapply(1:length(qal.all), function(oo) Reduce('+', cost.all[[oo]])/tab.diff)
  
  ss<-seq(1,25000,by = 250)
  for(hh in 1:length(ss)) {
    
    #cov.name3<-c('0%','30%','55%','70%')
    #AC.all<-lapply(1:length(qal.all), function(oo) accept.curve(ss[hh], qal.all[[oo]], cost.all[[oo]]) )
    #AC.all2<-lapply(1:length(qal.all), function(pp) data.frame(rep(cov.name2[pp],n.samples), Reduce('+',AC.all[[pp]])/tab.diff, rowSums(Reduce('+',AC.all[[pp]])/tab.diff)))
    
    ttest0<-(-qal.m[[1]]+0)-(-cost.m[[1]]-0)/ss[hh]
    ttesti<-list()
    #QALYS gained-additional cost/WTP
    for(oo in 2:intv.number) {ttesti[[oo]]<-(-qal.m[[oo]]+qal.m[[1]])-(-cost.m[[oo]]+cost.m[[1]])/ss[hh]}
    #for(oo in 2:intv.number) {ttesti[[oo]]<-(qal.m)-(cost.m)/ss[hh]}
    
    OC.all<-list(ttest0, ttesti[[2]], ttesti[[3]], ttesti[[4]], ttesti[[5]], ttesti[[6]], ttesti[[7]], ttesti[[8]])
    
    #OC.all2<- lapply(OC.all, rowSums)
    n.samples<-length(OC.all[[1]])
    OC.all3<-lapply(lapply(1:length(OC.all), function(pp) data.frame(rep(as.character(inv.names2[pp],n.samples)), OC.all[[pp]])), unname)
    
    
    new.ac<-array(unlist(OC.all3), dim=c(n.samples, dim(OC.all3[[1]])[2], intv.number))
    
    #for(aa in 2:24){
    e.frame<-data.frame(new.ac[,2,1:intv.number])
    xx<-apply(e.frame, 1, function(z) {which.max(z)})  
    ex.frame<-matrix(rep(0,dim(e.frame)[1]*dim(e.frame)[2]), ncol=intv.number)
    for(cc in 1:length(xx))
    {ex.frame[cc,xx[cc]]<-1}
    outp<-data.frame('GB', ss[hh], t(colMeans(ex.frame)))
    
    ifelse(hh==1, outerf<-outp, outerf<-rbind(outerf, outp))
    #}
  }
  
  colnames(outerf)<-c('country', 'threshold',inv.names2)
  #colnames(outerf)<-c('country','threshold','0%','30%','55%','70%')
  
  #outer.g<-melt(outerf, id.vars=c('country','Age', 'threshold'))
  outer.g<-melt(outerf, id.vars=c('country', 'threshold'))
  
  #colnames(outer.g)<-c('country', 'Age', 'Threshold','Coverage','value')
  colnames(outer.g)<-c('country', 'Threshold','Intervention','value')
  #mk<-ggplot(data= outer.g[outer.g$Age=='total',], aes(x=as.numeric(as.character(Threshold)), y=as.numeric(as.character(value)), color=Coverage, linetype=Coverage))
  
  mk<-ggplot(data= outer.g, aes(x=as.numeric(as.character(Threshold)), y=as.numeric(as.character(value))))
  
  label<-paste('strain=all',paste0(', discount=', dcount.l[dcount], '%'), ', coverage=',cov.f[cov.var])
  
  #final graph
  OC3<-mk+labs(subtitle=label, fill='Coverage')+geom_line(aes(color=Intervention))+scale_color_manual(breaks=inv.names3, values=cols)+scale_y_continuous(breaks=scales::pretty_breaks(n=7), limits = c(0,1))+scale_x_continuous(breaks=scales::pretty_breaks(n=6))+geom_vline(xintercept = c(15000,20000), linetype=3)+theme(plot.subtitle=element_text(size=12, hjust=0.5, color="black"), axis.text=element_text(size=11), axis.title=element_text(size=14))+guides(fill=guide_legend(title="Coverage"))
  
  #OC3<-OC3+scale_color_manual(values=cols[1:8])
  return(OC3)
}

#arrange output by coverage and discount
optimal3<-NULL
optimal3<-ggarrange(
  optimal.strat.new(1,1)+rremove("xy.title"), optimal.strat.new(2,1)+rremove("xy.title"), optimal.strat.new(3,1)+rremove("xy.title"), 
  optimal.strat.new(1,2)+rremove("xy.title"), optimal.strat.new(2,2)+rremove("xy.title"), optimal.strat.new(3,2)+rremove("xy.title"),
  optimal.strat.new(1,3)+rremove("xy.title"), optimal.strat.new(2,3)+rremove("xy.title"), optimal.strat.new(3,3)+rremove("xy.title"),
  ncol=3, nrow=3, common.legend = T, legend='bottom')

graphs.ocpost2<-NULL
graphs.ocpost2<-annotate_figure(optimal3,
                            top = text_grob('Optimal Acceptability Curve per Strategy'),
                            bottom = text_grob('Willingness-to-pay Threshold (£)'),
                            left = text_grob('Probability Strategy is Optimal', rot=90))

ggsave(paste0("OC3.allstraincoveragevdcount", version,".png"),plot=graphs.ocpost2, width=12, height=7.5, units='in',device='png')



just55<-optimal.strat.new(2,3)+rremove('xy.title')
graphsjust5<-annotate_figure(just55,
                                top = text_grob('Optimal Acceptability Curve per Strategy'),
                                bottom = text_grob('Willingness-to-pay Threshold (£GBP)'),
                                left = text_grob('Probability Strategy is Optimal', rot=90))


ggsave(paste0("OC3.allstrain05535", version,".png"),plot=graphsjust5, width=7, height=4, units='in',device='png')


##########################################################################################################
# ACCEPTIBILITY CURVE for single intervention probability it is cost effective on WTP scale
##############################################################################################################################
c.names<-c('program','0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR', 'total')

cols2<-c("#F8766D","#CD9600","#7CAE00", "#00BE67","#C77CFF","#00BFC4","#00A9FF","#FF61CC" )
####calcualte probability function
prob.curve2<-function(threshold.val, cov.var, dcount)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage', cov.f[cov.var]))
  cea.tables<-list.files(pattern=glob2rx(paste0('GBICER.all.strains', cov.f[cov.var],'p', '*', 'd',dcount)))
  qal.diff<-cost.diff<-list()
  
  for(ii in 1:length(cea.tables)) 
  {
    #optimal strategy compares all strategy outcomes to each other and sums to 1
    load(cea.tables[[ii]]);
    
    qal.diff[[ii]]<-ICER.set2$`net QALY`;
    cost.diff[[ii]]<-ICER.set2$`net cost`;
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

#########main acceptability per program function
accept.strat.graph2<-function(cov.var, dcount)
{
  ss<-c(seq(0,10000,by = 150), seq(10000, 25000, by=250))
  for(hh in 1:length(ss)) {ifelse(hh==1, accept.out<-prob.curve2(ss[hh], cov.var, dcount), 
                                  accept.out<-rbind(accept.out, prob.curve2(ss[hh], cov.var, dcount)))}
  
  colnames(accept.out)<-c('country','Intervention','threshold' , 'total')
  mj<-ggplot(data=accept.out, aes(x=as.numeric(as.character(threshold)), y=as.numeric(as.character(total)),color=Intervention))+geom_line()+scale_color_manual(breaks=inv.names3, values = cols[2:8])
  
  label<-paste(paste0('strain=all'),',',paste0('discount=', dcount.l[dcount], '%'), ',',paste0('coverage=', cov.f[cov.var]))
  
  AC3<-mj+labs(title=NULL, subtitle=label, fill='Intervention')+scale_y_continuous(limits = c(0.9,1), breaks=scales::pretty_breaks(n=7))+scale_x_continuous(breaks=scales::pretty_breaks(n=5))+geom_vline(xintercept = c(15000,20000), linetype=3)+geom_hline(yintercept = c(1), linetype=3)+theme(plot.subtitle=element_text(size=12, hjust=0.5, color="black"), axis.text=element_text(size=12))+guides(fill=guide_legend(title="Intervention"))
  
  return(AC3)
}

accept.strat.graph2(2,3)

#arrange output
accept32<-NULL
accept32<-ggarrange(
  accept.strat.graph2(1,1)+rremove("xy.title"), accept.strat.graph2(2,1)+rremove("xy.title"), accept.strat.graph2(3,1)+rremove("xy.title"), 
  accept.strat.graph2(1,2)+rremove("xy.title"), accept.strat.graph2(2,2)+rremove("xy.title"), accept.strat.graph2(3,2)+rremove("xy.title"),
  accept.strat.graph2(1,3)+rremove("xy.title"), accept.strat.graph2(2,3)+rremove("xy.title"), accept.strat.graph2(3,3)+rremove("xy.title"),
  ncol=3, nrow=3, common.legend = T, legend='bottom')

graphs.acpost<-NULL
graphs.acpost<-annotate_figure(accept32, 
                            top = text_grob('Acceptability Curve per Strategy (All Strains)'),
                            bottom = text_grob('Willingness-to-pay Threshold (£)'),
                            left = text_grob('Probability Cost-Effective', rot=90))
ggsave(paste0("AC3.allstrainvcoverage",version,".png"),plot=graphs.acpost, width=12, height=9, units='in',device='png')




just55a<-accept.strat.graph2(2,3)+rremove('xy.title')
graphsjust5a<-annotate_figure(just55a,
                             top = text_grob('Acceptability Curve per Strategy'),
                             bottom = text_grob('Willingness-to-pay Threshold (£GBP)'),
                             left = text_grob('Probability Strategy is Cost-Effective', rot=90))


ggsave(paste0("AC3allstrain05535", version,".png"),plot=graphsjust5, width=7, height=4, units='in',device='png')



#####Optimal and general Acceptibility Curve

just55ao<-ggarrange(optimal.strat.new(2,3)+rremove('xy.title'), accept.strat.graph2(2,3)+rremove('xy.title'), nrow=1, ncol=2, common.legend = T, legend='bottom', labels=c('A', 'B'))

just55ao2<-annotate_figure(just55ao,
                             bottom = text_grob('Willingness-to-pay Threshold (£GBP)'),
                             left = text_grob('Probability', rot=90))


ggsave(paste0("OCAC3.allstrain05535", version,".png"),plot=just55ao2, width=9, height=5, units='in',device='png')


####################################################################################################################################
#QALY Gains Stratified by Age
######################################################################################################################
cov.var<-2 #55%
dcount<-3#3.5%
qal.all<-list()
cost.all<-list()
outerf<-list()

qaly.gains.all<-function(cov.var, dcount)
{

setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage',cov.f[cov.var]))
qal.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.qaly.lossall.strains',cov.f[cov.var], 'p*', 'd',dcount)))
cos.tables<-list.files(pattern=glob2rx(paste0(fname[i.country], 'tot.cost.lossall.strains',cov.f[cov.var], '*', 'd',dcount)))

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

for(i.program in 1:length(inv.names))
{
  
  #qal.all2<-qal.all[[i.program]]
  age.small<-list()
  risk.age<-list()
  prop.small<-list()


#for(iter in 1:2)
#{

  #redis<-get(paste0('aged.', i.vout), envir = .GlobalEnv)[[iter]]
  redis<-Reduce('+', Map('-',qal.all[[1]], qal.all[[i.program+1]]))/tab.diff
  
  #redis<-aged.cases[[iter]]
  
  age.small[[iter]]<-cbind(redis[,1]+redis[,2]+redis[,12]+redis[,13], #0-2, 
                           redis[,3]+redis[,14], #2-4,
                           redis[,4]+redis[,15], #5-11
                           redis[,5]+redis[,6]+redis[,16]+redis[,17],
                           redis[,7]+redis[,18],
                           redis[,8]+redis[,19],
                           redis[9]+redis[,20],#45-64
                           redis[10]+redis[,21]+redis[11]+redis[,22]#75+
  )
 
 colnames(age.small[[iter]])<-c('0-1 yr olds','2-4 yr olds','5-11 yr olds','12-16 yr olds','17-24 yr olds', '25+ year olds')

test<-melt(age.small)
test$L1<-factor(test$L1)
levels(test$L1)<-c(paste0('I',i.program))

#levels(test$L1)<-c('SQ',paste0('S',i.program))
#test2<-cbind(test, rep(inv.n,dim(test)[1]))
ifelse(i.program==1, agg.hosp<-test, agg.hosp<-rbind(agg.hosp, test))
}

colnames(agg.hosp)<-c('X1','Age','value','Strategy')

hosp.m2<-ddply(agg.hosp[,2:4], .(Age, Strategy), summarise, m1=mean(value), LB=quantile(value, 0.025, na.rm=F), UB=quantile(value, 0.975, na.rm=F))
#hosp.m2$CM<-factor(hosp.m2$CM, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
hosp.m2$Strategy<-factor(hosp.m2$Strategy, levels=c('S1','S2','S3','S4','S5','S6','S7'), labels = c(inv.names))
#relevel
hosp.m2$Strategy<-factor(hosp.m2$Strategy, levels=c(inv.names3 [2:8]), labels=c("I1: 2-4 yrs","I2: 5-11 yrs",        "I3: 12-16 yrs", "I4: 2-11 yrs",   "I5: 2-4 & 12-16 yrs", "I6: 5-16 yrs",  "I7: 2-16 yrs"))
#hosp.m2$Age<-factor(hosp.m2$Age, levels=c('0-1 yr olds', '2-4 yr olds', '5-11 yr olds', '12-16 yr olds', '17-24 yr olds', '25-44 yr olds', '45-64 yr olds', '65+ yr olds'), labels=c('0-1','2-4','5-11','12-16','17-24', '25-44', '45-64', '65+'))

hosp.m2$Age<-factor(hosp.m2$Age, levels=c('0-1 yr olds', '2-4 yr olds', '5-11 yr olds', '12-16 yr olds', '17-24 yr olds', '25-44 yr olds', '45-64 yr olds', '65+ yr olds'), labels=c('0-1 yr olds', '2-4 yr olds', '5-11 yr olds', '12-16 yr olds', '17-24 yr olds', '25-44 yr olds', '45-64 yr olds', '65+ yr olds'))


    epio<-ggplot(data = hosp.m2, aes(x=Strategy, y=m1/1000, color=Strategy, group=Strategy)) + geom_pointrange(aes(ymax=UB/1000, ymin=LB/1000), position=position_dodge(0.85))+scale_y_continuous(breaks=scales::pretty_breaks(n=7))+theme_bw()+facet_wrap(~Age, ncol=3, scales = 'free_y')+scale_color_manual(values=cols2[2:8])+labs(y=paste('QALYs Gained (Thousands)'), x='Strategies', title='QALYs Gained per Strategy', subtitle = paste('coverage=', paste0(cov.f[cov.var]*100,'%')))+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position = 'None')+geom_hline(yintercept = 0)
    
    #epio<-ggplot(data = hosp.m2, aes(x=Age, y=m1/1000, color=Strategy, group=Strategy)) + geom_pointrange(aes(ymax=UB/1000, ymin=LB/1000), position=position_dodge(0.85))+theme_bw()+facet_wrap(~Strategy, ncol=3)+scale_color_manual(values=cols2)+labs(y=paste('QALYs Gained (Thousands)'), x='Age Groups')+theme_bw()+theme(legend.position = 'None')+scale_y_continuous(breaks=c(seq(-30,30,10)))+theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_hline(yintercept = 0)
    

return(epio)
    }

qaly.gains.all(2,3)
ggsave(paste0("QALYgains0553.5",version,".png"),plot=qaly.gains.all(2,3), width=12, height=7, units='in',device='png')
ggsave(paste0("QALYgains303.5",version,".png"),plot=qaly.gains.all(1,3), width=12, height=7, units='in',device='png')
ggsave(paste0("QALYgains703.5",version,".png"),plot=qaly.gains.all(3,3), width=12, height=7, units='in',device='png')


qalydots<-NULL
qalydots<-lapply(1:3, function(x) qaly.gains.all(x,3))
# qaly.gains.all(1,3)+rremove("xy.title"), qaly.gains.all(2,3)+rremove("xy.title"), qaly.gains.all(3,3)+rremove("xy.title"))
m1<-marrangeGrob(qalydots, nrow=1, ncol=1)
ggsave(paste0("epidots.QALYsallstrain",version,".png"),plot=m1, width=8.5, height=11, units='in',device='png')

######################################################################################################################
#QALY Gains Stratified by Age
######################################################################################################################

###QALYS
n.samples<-1000
avg.peragestrat<-lapply(1:length(qal.all), function(pp) data.frame(as.factor(rep(inv.names2[pp],n.samples)), Reduce('+',qal.all[[pp]])/tab.diff, rowSums(Reduce('+',qal.all[[pp]])/tab.diff)))

avg.allages<-lapply(1:length(qal.all), function(pp) rowSums(Reduce('+',qal.all[[pp]])/tab.diff))

for(rug in 1:8) {
grapple.m<-round(median(as.numeric(avg.allages[[rug]])))
grapple.sd<-round(quantile(as.numeric(avg.allages[[rug]]), c(0.025, 0.975)))

grapple<-cbind(inv.names2[rug],grapple.m, unname(grapple.sd[1]), unname(grapple.sd[2]))
ifelse(rug==1, m.qaly<-grapple, m.qaly<-data.frame(rbind(m.qaly, grapple)))}

colnames(m.qaly)<-c('intervention','median QALYS','LCI', 'UCI')


##########################################################################################
#Epi Pulls for Age-stratified Deaths, Cases, GP Appointments, and Hospitalizations
#####################################################################################################################
cols2<-c("#FF61CC","#F8766D", "#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF" )


######loader
age.strat.loader<-function(vout, version, cov.var, dcount, strategy)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage', cov.f[cov.var]))
  #strategy.summary<-NULL
  age.tables<-list.files(pattern=glob2rx(paste0('GBagestrat.',vout,'.allstrainp*',cov.f[cov.var],dcount.l[dcount],version)))
  load(age.tables[strategy], envir = .GlobalEnv)
} 

#########v.out=variable (e.g. death)

epi.pulls.allstrains<-function(i.vout, i.version, cov.var, dcount)
{
  agg.hosp<-c()
for(i.program in 1:length(inv.names))
{
    age.strat.loader(vout=i.vout, version=i.version,cov.var, dcount, strategy=i.program)
    
    age.small<-list()
    risk.age<-list()
    prop.small<-list()
    
    for(iter in 1:2)
    {
      redis<-get(paste0('aged.', i.vout), envir = .GlobalEnv)[[iter]]
      
      #redis<-aged.cases[[iter]]
      
      age.small[[iter]]<-cbind(redis[,1]+redis[,2]+redis[,12]+redis[,13], #0-2, 
                               redis[,3]+redis[,14], #2-4,
                               redis[,4]+redis[,15], #5-11
                               redis[,5]+redis[,6]+redis[,16]+redis[,17],
                               redis[,7]+redis[,18],
                               redis[,8]+redis[,19]+redis[9]+redis[,20]+redis[10]+redis[,21]+redis[11]+redis[,22]#65+
      )
      
      
      
      #prop.small[[iter]]<-age.small[[iter]]/rowSums(age.small[[iter]])
      
      risk.age[[iter]]<-cbind(redis[,1]+redis[,12],
                              redis[,2]+redis[,13], 
                              redis[,3]+redis[,14],
                              redis[,4]+redis[,15],
                              redis[,5]+redis[,16],
                              redis[,6]+redis[,17],
                              redis[,7]+redis[,18],
                              redis[,8]+redis[,19],
                              redis[,9]+redis[,20],
                              redis[,10]+redis[,21],
                              redis[,11]+redis[,22]
      )
      
      colnames(age.small[[iter]])<-c('0-1 yr olds','2-4 yr olds','5-11 yr olds','12-16 yr olds','17-24 yr olds', '25+ year olds')
      
      #colnames(prop.small[[iter]])<-c('0-1 yr olds','2-4 yr olds','5-11 yr olds','12-16 yr olds','17-24 yr olds', '25-44 yr olds', '45-64 yr olds', '65-74 yr olds', '75+ year olds')
    }
    
    test<-melt(age.small)
    test$L1<-factor(test$L1)
    levels(test$L1)<-c('SQ',paste0('I',i.program))
    #test2<-cbind(test, rep(inv.n,dim(test)[1]))
    ifelse(i.program==1, agg.hosp<-test, agg.hosp<-rbind(agg.hosp, test))
  }
  
  colnames(agg.hosp)<-c('X1','Age','value','Strategy')
  
  hosp.m2<-ddply(agg.hosp[,2:4], .(Age, Strategy), summarise, m1=mean(value), LB=quantile(value, 0.025, na.rm=F), UB=quantile(value, 0.975, na.rm=F))
  #hosp.m2$CM<-factor(hosp.m2$CM, levels=c('GB', 'BE', 'DE', 'FI','FR', "IT", "LU", "NL" ,"PL", "PE" ))
  hosp.m2$Strategy<-factor(hosp.m2$Strategy, levels=c('SQ', 'I1','I2','I3','I4','I5','I6','I7'), labels = c(inv.names2))
  #relevel
  hosp.m2$Strategy<-factor(hosp.m2$Strategy, levels=c(inv.names3))
  hosp.m2$Age<-factor(hosp.m2$Age, levels=c('0-1 yr olds','2-4 yr olds','5-11 yr olds','12-16 yr olds','17-24 yr olds', '25+ year olds'))
  
  
  #graph  
  
  if(i.vout %in% c('hosp', 'GP', 'cases'))
  {
    if(i.vout=='hosp')
    {
      epio<-ggplot(data = hosp.m2, aes(x=Strategy, y=m1, color=Strategy, group=Strategy)) + geom_pointrange(aes(ymax=UB, ymin=LB), position=position_dodge(0.85))+theme_bw()+facet_wrap(~Age, ncol=2, scales='free_y')+scale_color_manual(values=cols2)+labs(y=paste('Average Number of \n Hospitalizations'), x='Intervention Strategy')+theme_bw()+theme(legend.position = 'None')+theme(axis.text.x = element_text(angle = 90, hjust = 1))
      
      epistack<-ggplot(data=hosp.m2)+geom_col(aes(x=Strategy,y=m1, fill=Age))+theme_bw()+scale_fill_brewer()+
        labs(title=paste('Average Number of Hospitalizations'), y='Average Count', x='Intervention Strategy', subtitle=paste('Coverage=', paste0(cov.f[cov.var],'%'), 'Discount=', dcount.l[dcount]), fill = 'Age Strata')+theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1, size=14), axis.text.y = element_text( size=14),axis.title.y = element_text(size=14), legend.position = 'right')
    }
      if(i.vout=='GP')
      {
        epio<-ggplot(data = hosp.m2, aes(x=Strategy, y=m1/1000, color=Strategy, group=Strategy)) + geom_pointrange(aes(ymax=UB/1000, ymin=LB/1000), position=position_dodge(0.85))+theme_bw()+facet_wrap(~Age, ncol=2, scales='free_y')+scale_color_manual(values=cols2)+labs(y=paste('Average Number of \n GP Visits (Thousands)'), x='Intervention Strategy')+theme_bw()+theme(legend.position = 'None')+theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        epistack<-ggplot(data=hosp.m2)+geom_col(aes(x=Strategy,y=m1/1000, fill=Age))+theme_bw()+scale_fill_brewer()+labs(title=paste('Average Number of GP Visits (Thousands)'), y='Average Count', x='Intervention Strategy', subtitle=paste('Coverage=', paste0(cov.f[cov.var],'%'), 'Discount=', dcount.l[dcount]), fill = 'Age Stratum')+theme_bw()+theme(axis.text.x = element_text(angle = 60, hjust = 1, size=14), axis.text.y = element_text( size=14),axis.title.y = element_text( size=14), legend.position = 'right')
      }
  
      if(i.vout=='cases')
      {
      epio<-ggplot(data = hosp.m2, aes(x=Strategy, y=m1/1000, color=Strategy, group=Strategy))+ 
        geom_pointrange(aes(ymax=UB/1000, ymin=LB/1000), position=position_dodge(0.85)) +theme_bw()+facet_wrap(~Age, ncol=2, scales='free_y')+scale_color_manual(values = cols2)+labs(y=paste('Average Number of \n', i.vout, '(Thousands)'), x='Intervention Strategy')+theme_bw()+theme(legend.position = 'None')+theme(axis.text.x = element_text(angle = 90, hjust = 1))

      epistack<-ggplot(data=hosp.m2)+geom_col(aes(x=Strategy,y=m1/1000, fill=Age))+theme_bw()+scale_fill_brewer()+labs(title=paste('Annual Number of', i.vout, '(Thousands)'), y='Average Count', x='Intervention Strategy', subtitle=paste('Coverage=', paste0(cov.f[cov.var],'%'), 'Discount=', dcount.l[dcount]), fill = 'Age Stratum')+theme(axis.text.x = element_text(angle = 60, hjust = 1, size=14), axis.text.y = element_text(size=14),axis.title.y = element_text(size=14), legend.position = 'right')
      }
      }else{
        epio<-ggplot(data = hosp.m2, aes(x=Strategy, y=m1, color=Strategy, group=Strategy))+ geom_pointrange(aes(ymax=UB, ymin=LB), position=position_dodge(0.85)) +theme_bw()+facet_wrap(~Age, ncol=2, scales='free_y')+scale_color_manual(values=cols2)+labs(y=paste('Average Number of Influenza Associated', paste0(i.vout,'s')), x='Intervention Strategy')+theme_bw()+theme(legend.position = 'None')+theme(axis.text.x = element_text(angle = 90, hjust = 1))
        
        epistack<-ggplot(data=hosp.m2)+geom_col(aes(x=Strategy,y=m1, fill=Age))+theme_bw()+scale_fill_brewer()+labs(title=paste('Annual Number of Influenza Associated', paste0(i.vout, 's')), y='Average Count', x='Intervention Strategy', subtitle=paste('Coverage=', paste0(cov.f[cov.var],'%'), 'Discount=', dcount.l[dcount]), fill = 'Age Strata')+theme(axis.text.x = element_text(angle = 60, hjust = 1, size=14), axis.text.y = element_text(size=14), axis.title.y = element_text(size=14), legend.position = 'right')
      }
  
  ##savers
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage', cov.f[cov.var]))  
  
  ggsave(paste0('epidots', i.vout,cov.f[cov.var],dcount.l[dcount],i.version,".png"),plot =epio , width=11, height=8.5, units='in',device='png')
  
  ggsave(paste0('epistack', i.vout,cov.f[cov.var],dcount.l[dcount],i.version,".png"),plot =epistack , width=11, height=8.5, units='in',device='png')
}  


#epi.pulls.allstrains(i.vout='death', i.version='u7', cov.var=2, dcount=3)
#epi.pulls.allstrains(i.vout='cases', i.version='u7', cov.var=2, dcount=3)
#epi.pulls.allstrains(i.vout='hosp', i.version='u7', cov.var=2, dcount=3)
#epi.pulls.allstrains(i.vout='GP', i.version='u7', cov.var=2, dcount=3)



for(vv in 1:length(cov.f))
{
  for(dd in 1:length(dcount.l))
  {
    epi.pulls.allstrains(i.vout='death', i.version='v1', cov.var = vv, dcount = dd)
    epi.pulls.allstrains('cases',i.version='v1', cov.var = vv, dcount = dd)
    epi.pulls.allstrains('hosp', i.version='v1', cov.var = vv, dcount = dd)
    epi.pulls.allstrains('GP',i.version='v1', cov.var = vv, dcount = dd)
  }
}


######################################################################################################################
####Incidence
########################################################################################################
version<-'v1'
cov.var<-2
cov.sens.loader(strain=2, end.cov=cov.f[cov.var], version='v1')
cov.sens.loader(strain=3, end.cov=cov.f[cov.var], version='v1')
cov.sens.loader(strain=1, end.cov=cov.f[cov.var], version='v1')


avg.incidenceH1N1<-lapply(1:length(uk.tableH1N1), function(pp) rowSums(Reduce('+',uk.tableH1N1[[pp]])/tab.diff))
avg.incidenceH3N2<-lapply(1:length(uk.tableH3N2), function(pp) rowSums(Reduce('+',uk.tableH3N2[[pp]])/tab.diff))
avg.incidenceB<-lapply(1:length(uk.tableB), function(pp) rowSums(Reduce('+',uk.tableB[[pp]])/tab.diff))

avg.incidence.all<-Map('+', Map('+', avg.incidenceH1N1, avg.incidenceH3N2), avg.incidenceB)

for(rug in 1:8) {
  inc.m<-median(as.numeric(avg.incidence.all[[rug]]))
  inc.u<-mean(as.numeric(avg.incidence.all[[rug]]))
  
  inc.sd<-quantile(as.numeric(avg.incidence.all[[rug]]), c(0.025, 0.975))
  
  inccr<-cbind(inv.names2[rug],inc.m, inc.u, unname(inc.sd[1]), unname(inc.sd[2]))
  ifelse(rug==1, m.inc<-inccr, m.inc<-data.frame(rbind(m.inc, inccr)))}

colnames(m.inc)<-c('intervention','median incidence','mean incidence','LCI', 'UCI')

######################################################################
###########Total COSTS
#####################################################################
avg.costallages<-lapply(1:length(cost.all), function(pp) Reduce('+',cost.all[[pp]])/tab.diff)


for(rug in 1:8) {
  grapple.m<-round(median(as.numeric(avg.costallages[[rug]]))/1e6)
  grapple.sd<-round(quantile(as.numeric(avg.costallages[[rug]]), c(0.025, 0.975))/1e6)
  
  grapple<-cbind(inv.names2[rug],grapple.m, unname(grapple.sd[1]), unname(grapple.sd[2]))
  ifelse(rug==1, m.cost<-grapple, m.cost<-data.frame(rbind(m.cost, grapple)))}

colnames(m.cost)<-c('intervention','median Cost','LCI', 'UCI')


######Net Costs######################################################################
setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage',cov.f[cov.var]))
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




##########ICER vs SQ######################################################################
#rug<-5

for(rug in 2:8) {
una.m<-median(unlist(Map('/',
          lapply(c(1:7,9:19), function(x) rowSums(Map('-',cost.all[[rug]], cost.all[[1]])[[x]])),
          Map('*', -1,lapply(c(1:7,9:19), function(x) rowSums(Map('-',qal.all[[rug]], qal.all[[1]])[[x]]))))) #qalys gained
)
una.sd<-quantile(unlist(Map('/',
                         lapply(c(1:7,9:19), function(x) rowSums(Map('-',cost.all[[rug]], cost.all[[1]])[[x]])),
                         Map('*', -1,lapply(c(1:7,9:19), function(x) rowSums(Map('-',qal.all[[rug]], qal.all[[1]])[[x]]))
))), c(0.025, 0.975))

frad<-cbind(inv.names2[rug],una.m, unname(una.sd[1]), unname(una.sd[2]))
ifelse(rug==2, m.icer2<-frad, m.icer2<-data.frame(rbind(m.icer2, frad)))
}

for(rug in 2:8) {
  frad.m<-median((avg.costallages[[1]]-avg.costallages[[rug]])/(-avg.allages[[1]]+avg.allages[[rug]]))
  frad.sd<-quantile(((avg.costallages[[rug]]-avg.costallages[[1]])/(-avg.allages[[rug]]+avg.allages[[1]])), c(0.025, 0.975))
  frad<-cbind(inv.names2[rug],frad.m, unname(frad.sd[1]), unname(frad.sd[2]))
  ifelse(rug==2, m.icer<-frad, m.icer<-data.frame(rbind(m.icer, frad)))
  }
  


#########################################################################################################
####### Establish Dominance
#########################################################################################################
cov.var<-2; dcount<-3
n.samples<-2000
Qalygain12 <- Map('+',Map('+',QALYdifferences(Dataset1H1N1, Dataset2H1N1, strain=1,discount=dcount.l[dcount]),
                          QALYdifferences(Dataset1H3N2, Dataset2H3N2, strain=2,discount=dcount.l[dcount])), 
                  QALYdifferences(Dataset1B, Dataset2B, strain=3,discount=dcount.l[dcount]))


######CE evaluation
  
icers.benefit<-function(dcount, cov.var)
{
setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', 'GB', 'coverage',cov.f[cov.var]))
  cea.tables<-list.files(pattern=glob2rx(paste0('GBICER.all.strains', cov.f[cov.var], 'p', '*', 'd',dcount)))
  
  qdscos.f<-qds.m<-cos.m<-cep<-cep1<-contour_95<-mean.ff<-mean.po<-not.ce<-ics<-ics.f<-NULL
  
  for(ff in 1:length(inv.names))
  {#ff=interventions. Strains are already chosen
    load(cea.tables[[ff]]) #H1N1
    
    cer<-ICER.set2 #grab ICERs
    
    qa<-Reduce('+', cer$`net QALY`)/tab.diff
    co<-Reduce('+', cer$`net cost`)/tab.diff
    ic<-round(t(c(mean(co/qa), unname(quantile(co/qa, c(0.025, 0.975))), sd(co/qa))))
    
    #m.qa<-rowSums(qa)
    #m.co<-rowSums(co)
    
    par<-length(qa)
    qdscos<-data.frame(cbind(rep(fname[i.country], par), rep(inv.names[[ff]], par), rep(dcount.l[dcount], par), qa, co))
    #cos<-data.frame(cbind(rep(fname[i.country], par), rep(inv.names[[ff]], par), rep(dcount.l[dcount], par), co))
    ics<-data.frame(cbind(fname[i.country], inv.names[[ff]], dcount.l[dcount], ic))
    
    colnames(qdscos)<-c('country','program','discount','qalys','costs')
    colnames(ics)<-c('country','program','discount','icer mu', 'icer LB', 'icerUB', 'icerSD')
    #colnames(qds)<-colnames(cos)
    
    ifelse(ff==1, qdscos.f<-qdscos, qdscos.f<-rbind(qdscos.f, qdscos))
    #ifelse(ff==1, cos.f<-cos, cos.f<-rbind(cos.f, cos))
    ifelse(ff==1, ics.f<-ics, ics.f<-rbind(ics.f, ics))
  } 
    #colnames(qds.f)<-c('country','program','discount','0-1 m_LR', '1 y_LR', '2-4 y_LR', '5-11 y_LR', '12-14 y_LR', '15-16 y_LR', '17-24 y_LR', '25-44 y_LR', '45-64 y_LR', '65-74 y_LR', '75+ y_LR', '0-1 m_HR',  '1 y_HR', '2-4 y_HR', '5-11y_HR', '12-14 y_HR', '15-16 y_HR', '17-24 y_HR', '25-44 y_HR','45-64 y_HR', '65-74 y_HR', '75+ y_HR','total')
    print(ics.f)
}

icers.benefit(3,2)
icers.benefit(2,2)
icers.benefit(1,2)
icers.benefit(3,3)
icers.benefit(3,1)




    colnames(qdscos.f)<-colnames(qdscos)
    qds.i<-melt(qdscos.f, id.vars = c('country', 'program','discount', 'qalys','costs'))
    
    llg<-data.frame(cbind(as.character(qds.i[,2]),as.numeric(as.character(qds.i$qalys)), as.numeric(as.character(qds.i$costs)))) 
    colnames(llg)<-c('program', 'qaly','cost')
    

mean.p<-ddply(as.data.frame(llg), .(program), summarize, Rate1=median(as.numeric(as.character(qaly))), Rate2=median(as.numeric(as.character(cost))), costLB=quantile(as.numeric(as.character(cost)), 0.025), costUB=quantile(as.numeric(as.character(cost)), 0.975), qalyLB=quantile(as.numeric(as.character(qaly)), 0.025), qalyUB=quantile(as.numeric(as.character(qaly)), 0.975)   )

mean.sd<-ddply(as.data.frame(llg), .(program), summarize, cost.sd=sd(as.numeric(as.character(cost))), qal.sd=sd(as.numeric(as.character(qaly))))

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
levels(mean.f$program)= c(levels(mean.f$program), 'Status Quo')

sq<-matrix(c('Status Quo',rep(0,9)), nrow=1)
mean.ff<-data.frame(rbind(sq, as.matrix(mean.f[order(mean.f$Rate2, decreasing=F),])))

not.ce<-mean.po[which(div <0 | div >20000),]
colnames(mean.ff)<-colnames(mean.f)<-c('program','Rate1','Rate2','costLB', 'costUB','qalyLB', 'qalyUB', 'icer', 'icerUB', 'icerLB')

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
  
  mean.f2<-data.frame(cbind(apply(mean.f[which(div2 >0&div2 <20000),],2,rev)[,1:7], round(div2[which(div2 >0&div2 <20000)]), round(div.ub2[which(div2 >0&div2 <20000)]), round(div.lb2[which(div2 >0&div2 <20000)])))

  levels(mean.f2$program)= c(levels(mean.f2$program), 'Status Quo')
  mean.ff2<-data.frame(rbind(sq, as.matrix(mean.f2[order(mean.f2$Rate2, decreasing=F),])))
  
  not.ce2<-rbind(not.ce, mean.f[which(rudi.icer2[,2]/rudi.icer2[,1] <0 | rudi.icer2[,2]/rudi.icer2[,1] >20000),])
  colnames(mean.ff2)<-c('program','Rate1','Rate2','icer')
  
}

if(any(dim(not.ce2)!=dim(not.ce)))
{
  counter<-counter+1
  ##3rd CE check
  #mean.po<-mean.p[order(mean.p$Rate2, decreasing = T),]
  rudi.icer3<-diff(cbind(rev(as.numeric(as.character(mean.ff2[,2]))), rev(as.numeric(as.character(mean.ff2[,3])))))
  div3<-rudi.icer3[,2]/rudi.icer3[,1]
  mean.f3<-cbind(mean.f2[which(div3 >0&div3 <20000),][,1:3], div3[which(div3 >0&div3 <20000)])
  levels(mean.f3$program)= c(levels(mean.f3$program), 'Status Quo')
  mean.ff3<-data.frame(rbind(sq, as.matrix(mean.f3[order(mean.f3$Rate2, decreasing=F),])))
  
  not.ce3<-rbind(not.ce2, mean.f[which(rudi.icer3[,2]/rudi.icer3[,1] <0 | rudi.icer3[,2]/rudi.icer3[,1] >20000),])
  colnames(mean.ff3)<-c('program','Rate1','Rate2','icer')
  
}

if(any(dim(not.ce3)!=dim(not.ce2)))
{
  counter<-counter+1
  ##4th CE check
  #mean.po<-mean.p[order(mean.p$Rate2, decreasing = T),]
  rudi.icer4<-diff(cbind(rev(as.numeric(as.character(mean.ff3[,2]))), rev(as.numeric(as.character(mean.ff3[,3])))))
  div3<-rudi.icer4[,2]/rudi.icer4[,1]
  mean.f4<-cbind(mean.f3[which(div4>0&div4 <20000),][,1:3], div4[which(div4 >0&div4 <20000)])
  levels(mean.f4$program)= c(levels(mean.f4$program), 'Status Quo')
  mean.ff4<-data.frame(rbind(sq, as.matrix(mean.f4[order(mean.f4$Rate2, decreasing=F),])))
  
  not.ce4<-rbind(not.ce3, mean.f[which(rudi.icer4[,2]/rudi.icer4[,1] <0 | rudi.icer4[,2]/rudi.icer4[,1] >20000),])
  colnames(mean.ff4)<-c('program','Rate1','Rate2','icer')
  
}

if(counter>1)
{
  mean.ff<-get(paste0('mean.ff', counter))
  not.ce<-get(paste0('not.ce', counter))
}

  
  
  
  
  #############  #############  #############  #############  #############  #############  #############
  #AVG net benefit all strains
  ################  #############  #############  #############  #############  #############  #############
  #returns number of vaccine doses per season

i.cov<-mL55
for(program in 1:8)
{
  status.quo.program<-alt.program.coverage(new.coverage = i.cov, strategy = 1, strain = 1)#status quo
  new.program<-alt.program.coverage(new.coverage = i.cov,strategy = program, strain = 1)
  
  #return
  vaccines.doses.new2<-t(t(new.program)*as.vector(popv[1:22]))
  vaccines.doses.statquo<-t(t(status.quo.program)*as.vector(popv[1:22]))
  
  ########## 7.1.2 #calculates cost per dose########
  
  vaccine.cost.statquo<-vaccine.cost<-vcost.gp<-vcost.school<-vcost.pharm<-list()
  
  for(tty in 1:dim(vaccines.doses.new2)[1])
  {
    rcost.LAIV.school<-matrix(rtriangle(dim(vaccines.doses.statquo)[2]*n.samples,a=17,b=25,c=20.14), nrow=n.samples) #cost of vaccine per dose per season
    #L7.24+L2.80+L7.64 (reimbursement, service, dose cost) avg nurse 2.80 pounds/100 doses at rate of 36lbs/hr and 20 doses per hr works 7.5 hours a day, 2.5 travel and set up time and lunch
    rcost.GP<-matrix(rtriangle(dim(vaccines.doses.statquo)[2]*n.samples,a=17,b=25,c=19.66) , nrow=n.samples)
    #L7.24+L2.25+L7.64 (reimbursement, service, dose cost)
    rcost.pharmacy<-matrix(rtriangle(dim(vaccines.doses.statquo)[2]*n.samples, a=14, b=22, c=17.29), nrow=n.samples)
    #Katie paper BMJ pharmacy
    
    prop65up.pharm<-0.065 
    prop2to4.pharm<-0.0
    propadults.LR.pharm<-0.055
    propadults.HRpharm<-0.056
    
    #vaccine cost statquo 
    vaccine.cost.statquo[[tty]] <-cbind(t(vaccines.doses.statquo[tty,1:3]*t(rcost.GP[,1:3])), 
                                        t(vaccines.doses.statquo[tty,4:6]*t(rcost.LAIV.school[,4:6])),
                                        t(t(rcost.GP[,7:9])*(1-propadults.LR.pharm)*vaccines.doses.statquo[tty,7:9]+
                                            t(rcost.pharmacy[,7:9])*(propadults.LR.pharm)*vaccines.doses.statquo[tty,7:9]), 
                                        t(t(rcost.GP[,10:11])*(1-prop65up.pharm)*vaccines.doses.statquo[tty,10:11]+
                                            t(rcost.pharmacy[,10:11])*(prop65up.pharm)*vaccines.doses.statquo[tty,10:11]),
                                        t(vaccines.doses.statquo[tty,12:22]*(1-propadults.HRpharm)*t(rcost.GP[,12:22]))+t(vaccines.doses.statquo[tty,12:22]*propadults.HRpharm*t(rcost.pharmacy[,12:22])))
    
    #vaccine cost new 
    vaccine.cost[[tty]] <-cbind(t(vaccines.doses.new2[tty,1:3]*t(rcost.GP[,1:3])), 
                                t(vaccines.doses.new2[tty,4:6]*t(rcost.LAIV.school[,4:6])),
                                t(t(rcost.GP[,7:9])*(1-propadults.LR.pharm)*vaccines.doses.new2[tty,7:9]+
                                    t(rcost.pharmacy[,7:9])*(propadults.LR.pharm)*vaccines.doses.new2[tty,7:9]), 
                                t(t(rcost.GP[,10:11])*(1-prop65up.pharm)*vaccines.doses.new2[tty,10:11]+
                                    t(rcost.pharmacy[,10:11])*(prop65up.pharm)*vaccines.doses.new2[tty,10:11]),
                                t(vaccines.doses.new2[tty,12:22]*(1-propadults.HRpharm)*t(rcost.GP[,12:22]))+
                                  t(vaccines.doses.new2[tty,12:22]*propadults.HRpharm*t(rcost.pharmacy[,12:22])))
    
    
    vcost.gp[[tty]]<-rowSums(vaccines.doses.new2[tty,1:3]*rcost.GP[,1:3])+
                      rowSums(rcost.GP[,7:9]*(1-propadults.LR.pharm)*vaccines.doses.new2[tty,7:9])+
                      rowSums(rcost.GP[,10:11]*(1-prop65up.pharm)*vaccines.doses.new2[tty,10:11])+
                      rowSums(vaccines.doses.new2[tty,12:22]*(1-propadults.HRpharm)*rcost.GP[,12:22])
    
    vcost.pharm[[tty]]<-rowSums(rcost.pharmacy[,10:11]*(prop65up.pharm)*vaccines.doses.new2[tty,10:11])+
                        rowSums(rcost.pharmacy[,7:9]*(propadults.LR.pharm)*vaccines.doses.new2[tty,7:9])+
                        rowSums(vaccines.doses.new2[tty,12:22]*propadults.HRpharm*rcost.pharmacy[,12:22])

    vcost.school[[tty]]<-rowSums(vaccines.doses.new2[tty,4:6]*rcost.LAIV.school[,4:6])
  }
  
  
  
  m.vaccine.cost<-mean(rowSums(Reduce('+',vaccine.cost)/tab.diff))
  sd.vaccine.cost<-as.vector(quantile(rowSums(Reduce('+',vaccine.cost)/tab.diff), c(0.025, 0.975)))
  
  m.vcost.pharm<-mean(Reduce('+',vcost.pharm)/tab.diff)
  sd.vcost.pharm<-as.vector(quantile(Reduce('+',vcost.pharm)/tab.diff, c(0.025, 0.975)))
  
  m.vcost.school<-mean(Reduce('+',vcost.school)/tab.diff)
  sd.vcost.school<-as.vector(quantile(Reduce('+',vcost.school)/tab.diff, c(0.025, 0.975)))
  
  m.vcost.gp<-mean(Reduce('+', vcost.gp)/tab.diff)
  sd.vcost.gp<-as.vector(quantile(Reduce('+',vcost.gp)/tab.diff, c(0.025, 0.975)))
  
  ifelse(program==1, 
         vacc.out<-cbind(inv.names2[program], round(m.vaccine.cost), round(sd.vaccine.cost[1]), round(sd.vaccine.cost[2]),mean(rowSums(vaccines.doses.statquo)), mean(rowSums(vaccines.doses.new2)), round(m.vcost.gp), round(sd.vcost.gp)[1],round(sd.vcost.gp)[2], round(m.vcost.pharm), round(sd.vcost.pharm)[1],round(sd.vcost.pharm)[2], round(m.vcost.school), round(sd.vcost.school)[1], round(sd.vcost.school)[2]), 
         vacc.out<-rbind(vacc.out, cbind(inv.names2[program], round(m.vaccine.cost), round(sd.vaccine.cost[1]), round(sd.vaccine.cost[2]), mean(rowSums(vaccines.doses.statquo)), mean(rowSums(vaccines.doses.new2)), round(m.vcost.gp), round(sd.vcost.gp)[1],round(sd.vcost.gp)[2], round(m.vcost.pharm), round(sd.vcost.pharm)[1],round(sd.vcost.pharm)[2], round(m.vcost.school), round(sd.vcost.school)[1], round(sd.vcost.school)[2])))
  
  colnames(vacc.out)<-c('strategy', 'mean vaccine cost', 'sd vaccine cost SQ', 'sd.vaccine.cost.I', 'doses SQ',
                        'doses new', 'mean gp vacc cost', 'sd gp vacc cost L','sd gp vacc cost U',
                        'mean pharm vacc cost', 'sd pharm vacc cost L', 'sd pharm vacc cost U',
                        'mean school vacc cost', 'sd school vacc cost L', 'sd school vacc cost U')
    }
  


#########################################################################################################
#Net Benefit
#########################################################################################################

  
nb.curve<-function(threshold.val, cov.var, dcount)
{
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage', cov.f[cov.var]))
  cea.tables<-list.files(pattern=glob2rx(paste0('GBICER.all.strains', cov.f[cov.var],'p', '*', 'd',dcount)))
  qal.diff<-cost.diff<-list()
  
  for(ii in 1:length(cea.tables)) 
  {
    #optimal strategy compares all strategy outcomes to each other and sums to 1
    load(cea.tables[[ii]]);
    
    qal.diff[[ii]]<-ICER.set2$`net QALY`;
    cost.diff[[ii]]<-ICER.set2$`net cost`;
    
    #afford<-sum(unlist(ICER.set2$icer)<10000)/(tab.diff*length(ICER.set2$icer[[1]]))
    #ifelse(ii==1, afford.out<<-c(inv.names[ii], afford), afford.out<<-rbind(afford.out, c(inv.names[ii], afford)))
  }
  
  inv.names2<-c('Status Quo',"Preschool"   ,                "Primary"      ,               "Secondary"  ,                 "Preschool+Primary School"   , "Preschool+Primary+Secondary" ,"Preschool+Secondary"   ,      "Primary+Secondary")
  
  qal2<-lapply(1:7, function(x) Map('*', qal.diff[[x]], as.numeric(threshold.val)))
  nmb<-lapply(1:7, function(x) Map('-',qal2[[x]], cost.diff[[x]]))
  nmb2<-lapply(1:7, function(x) as.matrix(Reduce('+',nmb[[x]])/tab.diff))
  #nmb3<-lapply(1:7, function(x) cbind(nmb2[[x]], rowSums(nmb2[[x]])))
  
  nmb4<-lapply(1:7, function(x) mean(nmb2[[x]]))
  nmb.sd<-lapply(1:7, function(x) quantile(nmb2[[x]], c(0.025, 0.975)))
  nmb5<-cbind(do.call(rbind, nmb4), do.call(rbind, nmb.sd))
  
  curve.data.strat<-data.frame(cbind(rep('GB', length(cea.tables)), inv.names, rep(threshold.val,length(cea.tables)), nmb5))
  
  return(curve.data.strat)
}

  #changed discount

nb.curve(15000,2,3)


  nb.curve(20000,2,3)
  nb.curve(20000,2,2)  
  nb.curve(20000,2,1) 
  #changed coverage
  nb.curve(20000,3,3)  
  nb.curve(20000,1,3)  
  
 
  
  
  ss<-seq(1000,25000,by = 2000)
  for(hh in 1:length(ss)) {ifelse(hh==1, accept.out<-prob.curve2(ss[hh], cov.var, dcount), 
                                  accept.out<-rbind(accept.out, prob.curve2(ss[hh], cov.var, dcount)))}
  
  colnames(accept.out)<-c('country','Intervention','threshold' , 'total')
  
  
   
  icer.curve<-function(cov.var, dcount)
  {
    setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country]), 'coverage', cov.f[cov.var]))
    cea.tables<-list.files(pattern=glob2rx(paste0('GBICER.all.strains', cov.f[cov.var],'p', '*', 'd',dcount)))
    qal.diff<-cost.diff<-icey<-list()
    
    for(ii in 1:length(cea.tables)) 
    {
      #optimal strategy compares all strategy outcomes to each other and sums to 1
      load(cea.tables[[ii]]);
      
      icey[[ii]]<-lapply(1:19, function(x) ICER.set2$icer[[x]])
      
      qal.diff[[ii]]<-lapply(1:19, function(x) ICER.set2$`net QALY`[[x]])
      cost.diff[[ii]]<-lapply(1:19, function(x) ICER.set2$`net cost`[[x]])
    }
    nmb2<-lapply(1:7, function(x) mean(unlist(cost.diff[[x]]))/mean(unlist(qal.diff[[x]])))
    
    icer.sd<-matrix(unlist(lapply(1:7, function(x) quantile(mean(unlist(cost.diff[[x]])), c(0.025, 0.975))))/unlist(lapply(1:7, function(x) quantile(mean(unlist(qal.diff[[x]])), c(0.025, 0.975)))), ncol=2)
    icer2<-cbind(inv.names,data.frame(unname(unlist(lapply(1:7, function(x) mean(unlist(icey[[x]])))))), icer.sd)
    
    colnames(icer2)<-c('strategy', 'icer','2.5%', '97.5%')
    
    return(icer2)
  }

  cbind(unlist(icer.curve(1,3)), unlist(icer.curve(2,3)),unlist(icer.curve(3,3)))
  
 