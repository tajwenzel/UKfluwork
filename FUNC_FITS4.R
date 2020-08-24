#Based on Marc Baguelin's 2013 paper, written in R by Dr. Edwin van Leevan, modified by NWenzel, University of Washington School of Public Health, Department of Epidemiology. August 2016

#FITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFITTINGFIT

# have a go at replace  year missing virological data with older years 2007-2009.


#######################################################################################################
#.libPaths('/home/nwenzel/UKflu')
#library('devtools')
#install.packages('foreach',lib='~/UKflu')
#install.packages('doParallel',lib='~/UKflu')
#install.packages('parallel',lib='~/UKflu')
#install.packages('iterators',lib='~/UKflu')
#install.packages('doSNOW',lib='~/UKflu')
#install.packages('data.table',lib='~/UKflu')
#install.packages('devtools',lib='~/UKflu')
#install.packages('Rcpp',type='source')
#install.packages('fluEvidenceSynthesis',lib='~/UKflu')

library('fluEvidenceSynthesis')
library(doParallel)
library(foreach)
#library(stringi)
library(dplyr)
library(iterators)
library(pander)
library(data.table)
library('Rcpp')
library('base')
#library('snow')
#library(doSNOW)

rm(list = ls())


if(Sys.info()[7]=='Natasha') {clustertf<-0} else {(clustertf<-1)}

##FILE PATH & CLEAN UP
if(clustertf==0)
    {setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata")}else{
        numCores <- 8
        cl <- makeCluster(numCores)
        registerDoParallel(cl)
        .libPaths('~/UKflu')
        setwd('~/UKflu')
        set.seed(round(runif(n=1, min=1,max=500)))
          }


          load('ili.counts.rda',.GlobalEnv)
          load('virological.rda',.GlobalEnv)
  
  
#load('/home/nwenzel/UKflu/ili.counts.rda',.GlobalEnv)
#load('/home/nwenzel/UKflu/virological.rda',.GlobalEnv)

  #####################################################################################
  # INPUT
  #####################################################################################

data("age_sizes")
age.group.limits<-c(1,2,5,12,15,17,25,45,65,75) #upper limits
  
risk.ratios.ce<-matrix(c(0.021,0.055,0.055,0.098,0.098,0.098,0.087,0.092,0.183,0.45,0.45,0,0,0,0,0,0,0,0,0,0,0),ncol=length(age.group.limits)+1, byrow=TRUE)  
  
age.group.sizes<-stratify_by_age(age_sizes$V1,age.group.limits)

country.name<-('GB')
polymod<-fread(input=paste0(country.name,'table11.csv'),sep = 'auto');
polymod<<-cbind(polymod,polymod$V12)

if(clustertf==0)
  {polymod<-fread(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Inputdata', paste0(country.name,'table11.csv')))
polymod<<-cbind(polymod,polymod$V12)}

colnames(polymod)=c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13')
  current.contact.ids <<- seq(1,nrow(polymod))
  proposed.contact.ids <<- current.contact.ids
  
  if(clustertf==0){vstrategy<-dget(file='/Users/Natasha/Dropbox/UKfluworkGIT/Inputfunctions/FUNC_cov_strategy4.R')}  
  vstrategy<-dget(file='FUNC_cov_strategy4.R')  
  load('coverageH1',.GlobalEnv);
  coverageH1<-cov.eff.new;
  load('coverageH3',.GlobalEnv); 
  coverageH3<-cov.eff.new;
  load('coverageB',.GlobalEnv); 
  coverageB<-cov.eff.new;
  cov.eff.new<-NULL
  cov.eff.data<-list(coverageH1,coverageH3,coverageB)
  
  ########################################################################################
# Seasonal vaccination plan, using function strain.choice
 ########################################################################################
  #setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata")
  #season.choice<-dget('/Users/Natasha/Dropbox/UKfluworkGIT/Inputfunctions/FUNC_seas_fit.R')
  setwd("~/UKflu")
  season.choice<-dget('FUNC_seas_fit.R')

strains<-2L
seasons<-15
jam<-expand.grid(seasons,strains)

mapply(FUN=season.choice,jam[[1]],jam[[2]],program=1,country.name)
