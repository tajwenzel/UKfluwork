#Based on Marc Baguelin's 2013 paper, written in R by Dr. Edwin van Leevan, modified by NWenzel, University of Washington School of Public Health, Department of Epidemiology. July 2016

#####Load in package from github
#library(devtools)
#install_github("MJomaba/flu-evidence-synthesis", force=TRUE)
###load package
library(fluEvidenceSynthesis)
library(plyr)
library(pander)
library(data.table)

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
temp = list.files(pattern="*table.csv")
allpolymod = lapply(temp, fread)

#polymod<-allpolymod[[1]]
#i.top=1
for(i.top in 8:length(allpolymod))
{
  data("age_sizes")
polymod<-allpolymod[[i.top]]; #country levels BE DE FI GB IT LU NL PL
#polymod<-allpolymod[[8]]
#ilitest<-dget(file='INPUT_ilisample.R',keep.source=TRUE)

#rawpolymod<-read.csv('polymod_contacts.csv',header=TRUE,sep=',') #which dataset
current.contact.ids <- seq(1,nrow(polymod))
proposed.contact.ids <- current.contact.ids
initial.parameters <- dget(file="INPUT_UKinitial.R", keep.source=TRUE) #unnkowns and known
vstrategy<-dget(file='FUNC_vstrategy.R',keep.source=TRUE)

###################################################
age.group.limits <- c(1,5,15,25,45,65)
risk.ratios.null<-matrix(c(rep(0,7)), ncol = 7, byrow = T)
vcalendar<-vstrategy(risk.ratios.null)

age.group.sizes<<-stratify_by_age(age_sizes$V1,c(1,5,15,25,45,65))
confirmed.samples<-dget(file='/Users/Natasha/Dropbox/UKfluwork/INPUT_confirmedsamp.R',keep.source=TRUE)

# Fraction of each age group classified as high risk. Additional risk groups can be added here, as additional rows in our risk.ratios matrix)
buildLL<-dget(file='FUNC_LL.R',keep.source=TRUE)
llikelihood<-buildLL()

####################################################################################################
#Likelihood constant
####################################################################################################

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
burnin<-2000; #potatoes
out<-10000; #meat
saveiteration<-1; #thin size

contact.ids<-list()
# Run adaptive.mcmc
ptm <- proc.time() 
mcmc.result <- adaptive.mcmc(lprior = llprior, llikelihood=llikelihood, 
                             nburn=burnin,
                             initial = initial.parameters,
                             nbatch = out, blen = saveiteration,
                             outfun= function() {contact.ids[[length(contact.ids)+1]]<<-current.contact.ids},
                             acceptfun= function() {current.contact.ids <<- proposed.contact.ids}, agegrouplimits=age.group.limits, agegroupsizes=age.group.sizes, riskratios=risk.ratios.null)
proc.time() - ptm
  
sname<-c('BEid', 'DEid', 'FIid', 'GBid', 'ITid', 'LUid', 'NLid', 'PLid')
fwrite(contact.ids,paste0(sname[i.top],".csv"),row.names = FALSE)

}
