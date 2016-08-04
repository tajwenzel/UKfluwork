#Based on Marc Baguelin's 2013 paper, written in R by Dr. Edwin van Leevan, modified by NWenzel, University of Washington School of Public Health, Department of Epidemiology. July 2016

#####Load in package from github
#library(devtools)
#install_github("MJomaba/flu-evidence-synthesis", force=TRUE)
###load package
library(fluEvidenceSynthesis)
library(plyr)
library(pander)


#####serological cases < ili as they are a subset. create warning and stop evaluation
##FILE PATH & CLEAN UP
rm(list = ls())
setwd("/Users/Natasha/Dropbox/UKfluwork")
data("age_sizes")

#country levels BE DE FI GB IT LU NL PL
temp = list.files(pattern="*table.csv")
allpolymod = lapply(temp, read.csv)
polymod<-allpolymod[[4]];
#Load data
#data("age_sizes")
#data("polymod_uk")


#ilitest<-dget(file='INPUT_ilisample.R',keep.source=TRUE)
#confirmed.samples<-dget(file='INPUT_confirmedsamp.R',keep.source=TRUE)

#rawpolymod<-read.csv('polymod_contacts.csv',header=TRUE,sep=',') #which dataset
current.contact.ids <- seq(1,nrow(polymod))
proposed.contact.ids <- current.contact.ids
initial.parameters <- dget(file="INPUT_UKinitial.R", keep.source=TRUE) #unnkowns and known
vstrategy<-dget(file='FUNC_vstrategy.R',keep.source=TRUE)


###SOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVER
#dply by participant, age group and sum 
###MCMC run iteration 
burnin<-100; #potatoes
out<-200; #meat
saveiteration<-1; #keeper

###################################################ALTERNATE
age.group.limits <- c(1,5,15,25,45,65)
risk.ratios.null<-matrix(c(rep(0,7)), ncol = 7, byrow = T)
vcalendar<-vstrategy(risk.ratios.null)


age.group.sizes <- stratify_by_age(age_sizes$V1,c(1,5,15,25,45,65))
# Fraction of each age group classified as high risk. Additional risk groups can be added here, as additional rows in our risk.ratios matrix)
buildLL<-dget(file='FUNC_LL.R',keep.source=TRUE)
llikelihood<-buildLL()

#agegrouplimits<-age.group.limits;
#agegroupsizes<-age.group.sizes;
#riskratios<-risk.ratios.null

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


contact.ids<-list()
# Run adaptive.mcmc
ptm <- proc.time() 
mcmc.result <- adaptive.mcmc(lprior = llprior, llikelihood, 
                             nburn=burnin,
                             initial = initial.parameters,
                             nbatch = out, blen = saveiteration,
                             outfun= function() {contact.ids[[length(contact.ids)+1]]<<-current.contact.ids},
                             acceptfun= function() {contact.ids[[length(contact.ids)+1]]<<-current.contact.ids},
                             pars=initial.parameters, agegrouplimits=age.group.limits, agegroupsizes=age.group.sizes, riskratios=risk.ratios.null)
proc.time() - ptm

