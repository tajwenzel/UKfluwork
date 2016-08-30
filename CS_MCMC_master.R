#Based on Marc Baguelin's 2013 paper, written in R by Dr. Edwin van Leevan, modified by NWenzel, University of Washington School of Public Health, Department of Epidemiology. August 2016

#Cost Effectivness
#Should take output from ODE and calculate cost. FUNC_vstrategy will also change for different strategies. 
#Different strategies being examined

####0) STATUS QUO is vaccination of elderly (65+) and high risk proportion (ages 0-85)

####1) Preschool Vaccination only; ages 2-4
####2) Primary School Vaccination only; ages 5-11
####3) Secondary school vaccination only; ages 12-16
#---------------------------------------------------------
####4) Preschool+Primary School Vaccination; ages 2-11
####5) Preschool+Primary+Secondary school vaccination; ages 2-16
####6) Preschool+Secondary school vaccination; ages 2-4, 12-16
####7) Primary+Secondary school vaccination; ages 4-16


###Analysis Questions
#year by year cost-effectiveness, versus cost-effectiveness over time because of strain protection and shifting disease up to other age groups (e.g. Greece/Rubella example). Long term simulation is necessary too. 
#
#######################################################################################################
library(devtools)
#install.packages('flu-evidence-synthesis',type='source',repos='https://api.github.com/repos/MJomaba/flu-evidence-synthesis/zipball/master')
#install_github("MJomaba/flu-evidence-synthesis", force=TRUE)
###load package
library(fluEvidenceSynthesis)
library(plyr)
library(pander)
#install.packages('data.table', type = 'source',repos = 'http://Rdatatable.github.io/data.table')
library(data.table)
#library(parallel)
library(gdata)
#library(Rcpp)
#Sys.setenv("PKG_CXXFLAGS"="-fopenmp -std=c++11")
#library(pbapply)

##FILE PATH & CLEAN UP
rm(list = ls())
setwd("/Users/Natasha/Dropbox/UKfluworkGIT")

#agegrouplimits<-age.group.limits;
#agegroupsizes<-age.group.sizes;
#riskratios<-risk.ratios.null
#######################################################################
#IMPORT DATA
#####################################################################

#read in contact matricies for different countries


#polymod<-allpolymod[[1]]
#i.top=1
#for(i.top in 8:length(allpolymod))
#{


  
  ###################################################
  
  ####0) STATUS QUO is vaccination of elderly (65+) and high risk proportion (ages 0-85)
  
  ####1) Preschool Vaccination only; ages 2-4
  ####2) Primary School Vaccination only; ages 5-11
  ####3) Secondary school vaccination only; ages 12-16
  #---------------------------------------------------------
  ####4) Preschool+Primary School Vaccination; ages 2-11
  ####5) Preschool+Primary+Secondary school vaccination; ages 2-16
  ####6) Preschool+Secondary school vaccination; ages 2-4, 12-16
  ####7) Primary+Secondary school vaccination; ages 4-16
programme<-c(1:8) #see above for key
master.mcmc<-function(bb)
{

#for (bb in 1:length(programme))
#{
#####################################################################################
# INPUT
#####################################################################################

  data("age_sizes")
  
  initial.parameters <- dget(file="INPUT_UKInitial.R") #unknowns and known
  vstrategy<-dget(file='FUNC_cov_strategy.R')
  
  #respecify groups
  age.group.limits<-c(1,5,12,15,16,25,45,65,75) #upper limits
  
  
  risk.ratios.ce<-matrix(c(0.021,0.055,0.098,0.098,0.098,0.087,0.092,0.183,0.45,0.45,
                           0,0,0,0,0,0,0,0,0,0),ncol=length(age.group.limits)+1, byrow=TRUE)  # Fraction of each age group classified as high risk. Additional risk groups can be added here, as additional rows in our risk.ratios matrix)
  
age.group.sizes<<-stratify_by_age(age_sizes$V1,age.group.limits)

#-----------incidence for likelihood--------------------------------------------------
ILI<-fread(input='/Users/Natasha/Dropbox/UKfluworkGIT/ILIincidence.csv',sep2=',')
confirmed.samples<-ILI[-(1:3),]
confirmed.samples[[3]]<-NULL;
confirmed.samples[[2]]<-NULL;
confirmed.samples[[1]]<-NULL;
colnames(confirmed.samples)=c('V1','V2','V3','V4','V5','V6','V7','V8')

ili.array<-array(as.numeric(unlist(confirmed.samples)), c(52,8,19))
#index by [week,age-group,year(i)-1994]

#-----------reconfigure polymod--------------------------------------------------
polymod<-fread(input='/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/GBtable10.csv',sep = 'auto')
#here 65-74, and 75+ are combined
#repeat column V11 as we assume 65-74 and 75+ people have the same contact rates
polymod<<-cbind(polymod,polymod$V11)
colnames(polymod)=c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12')
current.contact.ids <<- seq(1,nrow(polymod))
proposed.contact.ids <<- current.contact.ids
########################################################################################
# Seasonal vaccination plan
########################################################################################
load('/Users/Natasha/Dropbox/UKfluworkGIT/coverageH1'); coverageH1<-cov.eff;
load('/Users/Natasha/Dropbox/UKfluworkGIT/coverageH3'); coverageH3<-cov.eff;
load('/Users/Natasha/Dropbox/UKfluworkGIT/coverageB'); coverageB<-cov.eff;
cov.eff<-NULL

strain<-2
scenario=1
season=4

  cov.eff.data<-list(coverageH1,coverageH3,coverageB)
  cov.eff.in<-cov.eff.data[[strain]]
  
for (season in 1:length(cov.eff.in))
{
vcalendar<-vstrategy(risk.ratios.ce, scenario,cov.eff.in,season) #seasons and strains are saved 3 separate datasets which have been manipulated into a list format by Nwenzel and saved under cov.effH1, cov.effH3, cov.effB.


####################################################################################################
#Likelihood constant
####################################################################################################
buildLL<-dget(file='FUNC_LL.R',keep.source=TRUE)
llikelihood<-buildLL()

  llprior <- function(pars) {
    if (any(pars[1:8] < 0) || any(pars[1:4] > 1) || any(pars[6:8] > 1)
        || pars[9] < log(0.00001) || pars[9] > log(10) )
      return(-Inf)
    
    lprob <- dnorm(pars[5], 0.1653183, 0.02773053, 1) 
    + dlnorm(pars[1], -4.493789, 0.2860455, 1) 
    + dlnorm(pars[2], -4.117028, 0.4751615, 1) 
    + dlnorm(pars[3], -2.977965, 1.331832, 1);
    
    return(lprob)
  }
  
  ####################################################################################################
  #MCMC run iteration 
  burnin<-30; #potatoes
  out<-150; #meat
  saveiteration<-1; #thin size
  
  contact.ids<-list()
  # Run adaptive.mcmc
  ptm <- proc.time() 
  mcmc.result <- adaptive.mcmc(lprior = llprior, llikelihood=llikelihood, 
                               nburn=burnin,
                               initial = initial.parameters,
                               nbatch = out, blen = saveiteration,
                               outfun= function() {contact.ids[[length(contact.ids)+1]]<<-current.contact.ids},
                               acceptfun= function() {current.contact.ids <<- proposed.contact.ids}, agegrouplimits=age.group.limits, agegroupsizes=age.group.sizes, riskratios=risk.ratios.ce,polymod=polymod)
  proc.time() - ptm
  #prog.name<-c(0:7)
  
  strain.name<-c('H1','H3','B')
  date.labels<-format(as.Date(cov.eff.in[[season]]$V1, origin="1970-01-01"), "%Y")
  
  save(contact.ids,file=paste0('flu',date.labels[[season]],strain.name[strain]))
  save(mcmc.result,file=paste0('flu',date.labels[[season]], strain.name[strain]))
    }
  }
}



sapply(programme, master.mcmc)
