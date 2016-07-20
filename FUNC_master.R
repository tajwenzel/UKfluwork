#Based on Marc Baguelin's 2013 paper, written in R by Dr. Edwin van Leevan, modified by NWenzel, University of Washington School of Public Health, Department of Epidemiology. July 2016

#####Load in package from github
library(devtools)
install_github("MJomaba/flu-evidence-synthesis", force=TRUE)
###load package
library(fluEvidenceSynthesis)
library(plyr)
library(pander)


#####serological cases < ili as they are a subset. create warning and stop evaluation
##FILE PATH & CLEAN UP
rm(list = ls())
setwd("/Users/Natasha/Dropbox/UKfluwork")

#Load data
data("age_sizes")
data("polymod_uk")
data("ili")
data("confirmed.samples")

#ilitest<-dget(file='INPUT_ilisample.R',keep.source=TRUE)
#confirmed.samples<-dget(file='INPUT_confirmedsamp.R',keep.source=TRUE)

#set initial parameters
initial.parameters <- dget(file="UKStart.R", keep.source=TRUE) #unnkowns and known
vstrategy<-dget(file='FUNC_vstrategy.R',keep.source=TRUE)

polymod<-polymod_uk
#rawpolymod<-read.csv('polymod_contacts.csv',header=TRUE,sep=',') #which dataset
current.contact.ids <- seq(1,nrow(polymod))
proposed.contact.ids <- current.contact.ids
###Now that vaccination schemes have been input, we will run the model to look at different scenarios. 

###SOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVERSOLVER

###MCMC run iteration 
burnin<-100; #potatoes
out<-100; #meat
saveiteration<-1; #keeper

###################################################ALTERNATE
age.group.limits <- c(1,5,15,25,45,65)
risk.ratios.null<-matrix(c(rep(0,7)), ncol = 7, byrow = T)
vaccine_calendar<-vstrategy(risk.ratios.null)

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
                             outfun<<-function() {contact.ids[[length(contact.ids)+1]]<<-current.contact.ids},
                             acceptfun<<-function() {
                               current.contact.ids <<-proposed.contact.ids},
                             nburn=burnin,
                             initial = initial.parameters,
                             nbatch = out, blen = saveiteration,...pars=initial.parameters,age.group.limits,age.group.sizes,risk.ratios.null)
proc.time() - ptm


###############Cut and Sort into outbreak size, age group, risk group

colnames(mcmc.result$batch) <- c("eps1", "eps2", "eps3", "psi", "q",
                                 "susc1", "susc2", "susc3", "I0")

ggplot(data=melt(mcmc.result$batch)) + facet_wrap( ~ Var2, ncol=3, scales="free" ) + geom_histogram(aes(x=value,factor=Var2), bins=25)


cim<-credible.interval.model(ode.results, mcmc.result$batch, intervals = c(0,0.5,0.975))
cim$row.ID<-as.Date(as.character(cim$row.ID))
#####################################Plot
library(ggplot2)
cases <- sapply(seq(1,nrow(mcmc.result$batch)), function(i)
{
  sum(vaccinationScenario( age_sizes=age_sizes[,1],
                           vaccine_calendar=vaccine_calendar,
                           polymod_data=as.matrix(polymod_uk),
                           contact_ids=mcmc.result$contact.ids[i,],
                           parameters=mcmc.result$batch[i,]
  ))
})
ggplot(data=data.frame("cases"=cases)) + geom_histogram(aes(x=cases),bins=25)
