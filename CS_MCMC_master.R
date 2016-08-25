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
  vstrategy<-dget(file='FUNC_vstrategy8.R')
  
  #respecify groups
  age.group.limits <- c(1,5,12,17,26,46,65) # 8 age groups
  #risk.ratios.null<-matrix(c(rep(0,7)), ncol = 7, byrow = T)
  
  risk.ratios.ce<-matrix(c(0.021, 0.021, 0.055, 0.098, 0.087, 0.092, 0.183, 0.45,0,0,0,0,0,0,0,0),ncol=(length(age.group.limits)+1),byrow=TRUE)  # Fraction of each age group classified as high risk. Additional risk groups can be added here, as additional rows in our risk.ratios matrix)
  
  vcalendar<-vstrategy(risk.ratios.ce, scenario)
  
  age.group.sizes<<-stratify_by_age(age_sizes$V1,age.group.limits)
  confirmed.samples<-dget(file='/Users/Natasha/Dropbox/UKfluwork/INPUT_confirmedsamp8.R')
  
  
  load("/Users/Natasha/Dropbox/UKfluworkGIT/MainCode/vaccine.calendars.Rda")
#Calculate the 8 age group polymod
#div8polymod<-dget(file='INPUT_polymod_pull8.R',keep.source=TRUE)
#div8polymod(age.group.sizes,touchtype=1)

#temp2 = list.files(pattern="*table8.csv")
#allpolymod = lapply(temp2, fread)
polymod<<-fread(input='GBtable8.csv',sep = 'auto')
current.contact.ids <<- seq(1,nrow(polymod))
proposed.contact.ids <<- current.contact.ids
########################################################################################
# Seasonal vaccination plan
########################################################################################
load('cov.effH1'); cov.effH1<-vset;
load('cov.effH3'); cov.effH3<-vset;
load('cov.effB'); cov.effB<-vset;
vset<-NULL

for (strain in 1:3)} #end


  for (season in 1:length(cov.effB))
    
vcalendar(risk.ratios.ce,scenario,season) #seasons and strains are saved 3 separate datasets which have been manipulated into a list format by Nwenzel and saved under cov.effH1, cov.effH3, cov.effB.

  






####################################################################################################
#Likelihood constant
####################################################################################################
buildLL<-dget(file='FUNC_LL8.R',keep.source=TRUE)
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
  
 2013/2014 live attenuated vaccine leave out 
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
  proc.time() - ptm}
    
  prog.name<-c(0:7)
  save(contact.ids,file=paste0('ctmatprog2',prog.name[bb]))
  save(mcmc.result,file=paste0('GBprog2', prog.name[bb]))
}



sapply(programme, master.mcmc)
