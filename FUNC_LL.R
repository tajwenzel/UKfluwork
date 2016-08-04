#Source file for fluEvidenceSynthesis for changing likelihood. Options for varying age group sizes, risk groups, (and susceptibility, and ascertainment through the starting parameters). Adapted from code by Baguelin 2013 and Edwin van Leevvan, by TajWenzel.

#pars<-initial.parameters;
#polymod<-set
#agegrouplimits<-age.group.limits
#agegroupsizes<-age.group.sizes
#riskratios<-risk.ratios.null


build.llikelihood<-function()
  {
llikelihood.f <- function(pars, agegrouplimits, agegroupsizes, riskratios,...)
  {
  proposed.contact.ids <<- current.contact.ids
    if (runif(1,0,1) < 0.1) {
      rs <<- round(runif(2,1,length(proposed.contact.ids)))
      proposed.contact.ids[rs[1]] <<- rs[2]
    }
    # Resample contact ids. T
    contacts<- contact.matrix(as.matrix(polymod[proposed.contact.ids,]),
                              age_sizes[,1], agegrouplimits)
    
    
    age.groups <- stratify_by_age(age_sizes[,1], agegrouplimits)
    
    # Population sizes in each age and risk group
    popv <- stratify_by_risk(age.groups,riskratios) #agegroups*risk groups
    
    
    epsilons <- c(pars[1], pars[1], pars[2], pars[2], pars[2],pars[3],pars[3]); #epsilon is         ascertainment
    
    initial.risk<-(rep(10^pars[9], length(riskratios)));
    initial.infected <- stratify_by_risk(initial.risk, riskratios) 
    odes <<- infectionODEs(popv, initial.infected,
                          vaccine_calendar,
                    contacts,c(pars[6], pars[6], pars[6],pars[7], pars[7], pars[7], pars[8]),
                          transmissibility = pars[5],
                          infection_delays=c(0.8,1.8), interval=7 ) #interval is in days
    
    
    # Ignore times row
    relevant.range<-(nrow(riskratios)*length(age.groups))+1
    converted.odes <- odes[,2:(relevant.range)];
    dateaxis<-odes[,1]
    
    #Convert age groups and risk groups
    #converted.odes[,1] <- rowSums(odes[,1])
    #converted.odes[,2] <- rowSums(odes[,c(2,3)])
    #converted.odes[,3] <- rowSums(odes[,c(4,5)])
    #converted.odes[,4] <- rowSums(odes[,c(6)])
    #converted.odes[,5] <- rowSums(odes[,c(7)])
    #converted.odes <- converted.odes[,1:5]
    
    matplot(dateaxis,converted.odes[,1:length(converted.odes)], type='l')
    legend('topleft',legend=1:length(converted.odes),col=1:length(converted.odes), pch=2)
    # For each week and each group sum log likelihood
    
    ll<-log_likelihood_cases(
      epsilons,pars[4], as.matrix(converted.odes),
      agegroupsizes, newili$ili, newili$total.monitored,
      newcs$positive, newcs$total.samples)
    return(ll)
  }
  return(llikelihood.f)
}

