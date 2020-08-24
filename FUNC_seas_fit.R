season.choice<-function(season, strain,program,country.name)
{
cov.eff.in<-cov.eff.data[[strain]]
vcalendar1<<-vstrategy(risk.ratios.ce,program,cov.eff.in,season)
master.name<-paste0(country.name,'flu.u2')
strain.name<-c('H1','H3','B')
date.labels<-format(as.Date(cov.eff.in[[season]]$V1, origin="1970-01-01"), "%Y")
ctname<-paste0(country.name,'ct.u2',date.labels[1],last(date.labels),strain.name[strain])
rname<-paste0(master.name,date.labels[1],last(date.labels),strain.name[strain])

####################################################################################################
resort<-list(virological$pos.by.strain$H1,virological$pos.by.strain$H3,virological$pos.by.strain$B)
positive<<-as.matrix(resort[[strain]])
total.sampled<<-as.matrix(virological$no.samples)

if(clustertf==0) {setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Inputfunctions/cluster/exact")}
#build.parameters<-dget(file='INPUT_UKInitial.R')
build.parameters <- dget(file="INPUT_INTLInitial.R")
initial.parameters<-build.parameters(country.name, clustertf, season, strain)
#initial.parameters<-build.parameters()
s.index<<-(((season-1)*52)+1):(season*52)
if(clustertf==0) {setwd("/Users/Natasha/Dropbox/UKfluworkGIT/Inputfunctions/cluster/exact")}
  #Likelihood constant
  ####################################################################################################
  
  if( season>14) {llikelihood<-dget(file='FUNC_LL_extended.R')}else{buildLL<-dget(file='FUNC_LL.R');
    llikelihood<-buildLL( )}
  
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
 
  burnin<-50000
  out<-200000
  saveiteration<-1
  
  contact.ids<-list()
  
  mcmc.result <- adaptive.mcmc(lprior = llprior, llikelihood=llikelihood, nburn=burnin,initial = initial.parameters,nbatch = out, blen = saveiteration,outfun= function() {contact.ids[[length(contact.ids)+1]]<<-current.contact.ids}, acceptfun= function() {current.contact.ids <<-proposed.contact.ids}, agegrouplimits=age.group.limits, agegroupsizes=age.group.sizes, riskratios=risk.ratios.ce, season=season, strain=strain, polymod=polymod)

  save(contact.ids,file=ctname)
  save(mcmc.result,file=rname)

 }

