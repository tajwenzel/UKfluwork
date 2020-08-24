season.choice<-function(sea,ii,program,country.name)
{  
cov.eff.in<-cov.eff.data[[ii]]
vcalendar<<-vstrategy(risk.ratios.ce,program,cov.eff.in,sea)
master.name<-paste0(country.name,'flu.p9')
strain.name<-c('H1','H3','B')  
date.labels<-format(as.Date(cov.eff.in[[sea]]$V1, origin="1970-01-01"), "%Y") 
ctname<-paste0(country.name,'ct.p9',date.labels[1],last(date.labels),strain.name[ii])
rname<-paste0(master.name,date.labels[1],last(date.labels),strain.name[ii])  

#############################################################
resort<-list(virological$pos.by.strain$H1,virological$pos.by.strain$H3,virological$pos.by.strain$B)
positive<<-as.matrix(resort[[ii]])
total.sampled<<-as.matrix(virological$no.samples)
 
#total.sampled<<-rbind(total.sampled, total.sampled[987,]) #week 34 repeat
s.index<<-(((sea-1)*52)+1):(sea*52)
#####################################################################################
# Posterior from previous year
#---must make allowance for year 1, 1995
strain.pull<-c(glob2rx(paste0(country.name,'flu.p9*H1')),glob2rx(paste0(country.name,'flu.p9*H3')),glob2rx(paste0(country.name,'flu.p9*B')))
list.files(pattern=strain.pull[ii])

#setwd('/Users/Natasha/Dropbox/UKfluworkGIT/cluster')
######################################################################################
build.parameters <- dget(file="INPUT_UKInitial.R")
initial.parameters<-build.parameters(sea)

#####################################################################################
#match to closests prior on virological data
chase<-c()
year.skips<-function(season)
  {s.index<-(((season-1)*52)+1):(season*52);
   chase<-rbind(chase,sum(positive[s.index,1:5]))}

skip.vec<-sapply(c(1:length(seasons)),year.skips)

skip.index.year<-(sea)-which.min(abs(skip.vec[1:(sea-1)] - skip.vec[sea])) 
#####################################################################################

if(sea > 2)
  {load(paste0(master.name,(as.numeric(date.labels[1])-skip.index.year),(as.numeric(last(date.labels))-skip.index.year),strain.name[ii]));  

sub.year<-(as.numeric(date.labels[1])-skip.index.year); 
print('going back 1 year');
print(sub.year);
 
#remove burn-in form pervious years before calculating new priors. Average burn-in is 60000
      rm.burnin<-round((0.2*dim(mcmc.result$batch)[1]))
  
      med.f<-function(cat) {median(mcmc.result$batch[rm.burnin:dim(mcmc.result$batch)[1],cat])}
      initial.parameters<<-sapply(c(1:dim(mcmc.result$batch)[2]),med.f)
      #med.pull<-sapply(c(1:dim(mcmc.result$batch)[2]),med.f)
   sd.f<-function(dog) {sd(mcmc.result$batch[rm.burnin:dim(mcmc.result$batch)[1],dog])} 
  sd.pull<-sapply(c(1:dim(mcmc.result$batch)[2]),sd.f)
       
   mean.f<-function(mouse) {mean(mcmc.result$batch[rm.burnin:dim(mcmc.result$batch)[1],mouse])} 
  mean.pull<-sapply(c(1:dim(mcmc.result$batch)[2]),mean.f)
  
      mcmc.result<-NULL
  }

####################################################################################################
  #Likelihood constant
  ####################################################################################################
  buildLL<-dget(file='FUNC_LL.R')
  llikelihood<-buildLL()
 
  #build.prior<-dget(file='FUNC_prior.R')
  #llprior<-build.prior()
  if(sea<=2)
  {
  llprior<-function(pars)
    {if (any(pars[1:8] < 0) || any(pars[1:4] > 1) || any(pars[6:8] > 1)
        || pars[9] < log(0.00001) || pars[9] > log(10))
      return(-Inf)
    
    lprob <- dnorm(pars[5], 0.1653183, 0.02773053, 1)
    + dlnorm(pars[1], -4.493789, 0.2860455, 1)
    + dlnorm(pars[2], -4.117028, 0.4751615, 1)
    + dlnorm(pars[3], -2.977965, 1.331832, 1)
    return(lprob)
    }
  } else {
    llprior<-function(pars)
    {
    options(warn=-1);
    if (any(pars[1:8] < 0) || any(pars[1:4] > 1) || any(pars[6:8] > 1)
        || pars[9] < log(0.00001) || pars[9] > log(10) )
      return(-Inf)
    
    parm5<-dnorm(pars[5], mean.pull[5], sd.pull[5], 1) 
    parm1<-dlnorm(pars[1], -4.493789, 0.2860455, 1)
    parm2<-dlnorm(pars[2], -4.117028, 0.4751615, 1)
    parm3<-dlnorm(pars[3], -2.977965, 1.331832, 1)
    #parm2<-dlnorm(pars[2], log(mean.pull[2]), abs(log(sd.pull[2])), 1) 
    #parm3<-dlnorm(pars[3], log(mean.pull[3]), abs(log(sd.pull[3])), 1)
    
    #lprob <- dnorm(pars[5], mean.pull[5], sd.pull[5], 1) 
    #+ dt(pars[1], rm.burnin-1, ncp=log(mean.pull[1]), 1) 
    #+ dt(pars[2], rm.burnin-1, ncp=log(mean.pull[2]), 1) 
    #+ dt(pars[3], rm.burnin-1,ncp=log(mean.pull[3]), 1)

    lprob<-parm5+parm1+parm2+parm3
    return(lprob)
    }
  }

  ####################################################################################################
 
  burnin<-0
  out<-150000
  saveiteration<-1
  
  contact.ids<-list()
  
#print(7)

mcmc.result <- adaptive.mcmc(lprior = llprior, llikelihood=llikelihood, nburn=burnin, initial = initial.parameters, nbatch = out, blen = saveiteration,
 outfun= function() {contact.ids[[length(contact.ids)+1]]<<-current.contact.ids}, acceptfun= function() {current.contact.ids <<-proposed.contact.ids},agegrouplimits=age.group.limits, agegroupsizes=age.group.sizes, riskratios=risk.ratios.ce, season=sea, strain=ii, polymod=polymod)


save(contact.ids,file=ctname)
save(mcmc.result,file=rname)
}



