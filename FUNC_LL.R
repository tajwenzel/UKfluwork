#Source file for fluEvidenceSynthesis for changing likelihood. Options for varying age group sizes, risk groups, (and susceptibility, and ascertainment through the starting parameters). Adapted from code by Baguelin 2013 and Edwin van Leevvan, by TajWenzel.
#riskratios<-risk.ratios.ce
#pars<-c(0.008093337, 0.01773221, 0.05030984, 3.170519e-03, 0.1670740, 0.3780008, 0.9224004, 0.5730906, -0.1608632)
#pars<-initial.parameters
#agegrouplimits<-age.group.limits
#########
#test.pars<- c(runif(1,0,0.01), runif(1,0.025,0.075), runif(1,0,0.00005), runif(1,0.01,0.25), runif(1,0.025,0.5), runif(1,0.5,1), runif(1,0.025,0.075),runif(1,0.25,0.75), runif(1,-0.4,-0.01))
#pars<-test.pars

#ll <- 0;
#for (j in 1:ncols(ili.cases) {
  #for (i in 1:nrow(ili.cases)) {
    #lp <- log_likelihood_cases(c(epsilon[j]), pars[4],
     #                          as.matrix(converted_odes[i,j]), c(age.groupLL[j]),
      #                         as.matrix(ili.cases[s.index,][i,j]), as.matrix(ili.total[s.index,][i,j]),
       #                        as.matrix(positive[s.index,][i,j]), as.matrix(total.samples[s.index,][i,j]))
    #print(lp)
    #print(converted.odes[i,j])
    #etc..
    #ll <- ll + lp
  #}
#}
###############################

data("age_sizes")

build.llikelihood<-function(season)
  {
llikelihood.f <- function(pars, agegrouplimits, agegroupsizes, riskratios,...)
  {
  proposed.contact.ids <<- current.contact.ids
    if (runif(1,0,1) < 1) {
      rs <- round(runif(2,1,length(proposed.contact.ids)))
      proposed.contact.ids[rs[1]] <<- rs[2]}
    # Resample contact ids.
    contacts<- contact_matrix(as.matrix(polymod[proposed.contact.ids,]),
                              age_sizes[,1], agegrouplimits)
    
    
    age.groups <- stratify_by_age(age_sizes[,1], agegrouplimits)
    
    # Population sizes in each age and risk group
    popv <- stratify_by_risk(age.groups,riskratios) #agegroups*risk groups
    #popv<-popv[1:length(riskratios)]
    #popv<-c(popv[12:22], popv[1:11],popv[23:33]) #high risk is first here
    
    epsilons <- c(pars[1], pars[1],pars[2],pars[2], pars[3]); #epsilon is ascertainment
    
    initial.risk<-(rep(10^pars[9], length(age.groups)));
    initial.infected <- stratify_by_risk(initial.risk, riskratios) 
    initial.infected <-c(initial.infected[12:22],initial.infected[1:11],initial.infected[23:33]) 

odes <- infectionODEs(popv, initial.infected,vaccine_calendar=vcalendar,contact_matrix=contacts,c(pars[6],pars[6],pars[6], pars[6], pars[6],pars[6],pars[7], pars[7], pars[7],pars[8], pars[8]),transmissibility = pars[5],infection_delays=c(0.8,1.8), interval=7) 
#interval is in days
    
    
    # Ignore times row
    #relevant.range<-(nrow(riskratios)*length(age.groups))+1
    #converted.odes <- odes[,2:(relevant.range)];
    converted.odes<-matrix(c(rep(0,52*5)),nrow=52,byrow=TRUE)
    dateaxis<-odes[,1]
    
    #Convert age groups and risk groups
    converted.odes[,1] <- rowSums(odes[,c(2,3,4,12,13,14)])
    converted.odes[,2] <- rowSums(odes[,c(5,6,15,16)])
    converted.odes[,3] <- rowSums(odes[,c(7,8,9,17,18,19)])
    converted.odes[,4] <- rowSums(odes[,c(10,20)])
    converted.odes[,5] <- rowSums(odes[,c(11,12,21,22)])
    converted.odes <- as.data.frame(converted.odes[,1:5])
    colnames(converted.odes)=c('V1','V2','V3','V4','V5')

    #matplot(dateaxis,converted.odes[,1:length(converted.odes)], type='l')
    #legend('topleft',legend=1:length(converted.odes),col=1:length(converted.odes), pch=2)
    # For each week and each group sum log likelihood
    
    #load in data from global
    #resort<-list(virological$pos.by.strain$H1,virological$pos.by.strain$H3,virological$pos.by.strain$B)
    #positive2<-as.matrix(resort[[ii]])
    #positive2<-positive2[-1,]; positive2<-rbind(positive, positive[987,]) #week 34 repeat
    
    colnames(ili.counts$ili)=c('V1','V2','V3','V4','V5')
    colnames(ili.counts$total.monitored)=c('V1','V2','V3','V4','V5')
    ili.cases<-as.matrix(ili.counts$ili)
    ili.total<-as.matrix(ili.counts$total.monitored)
    age.groupLL <- stratify_by_age(age_sizes[,1], c(5,15,45,65))
    
    
     #load('zero.compare')
    
    #if(sum(positive[s.index,1:5])<3) {ll<-0}
   #else {
    ll<-0
    ll<-log_likelihood_cases(
      epsilons,pars[4], as.matrix(converted.odes),
      age.groupLL, ili.cases[s.index,],ili_monitored=ili.total[s.index,],
     positive[s.index,], confirmed_samples=total.sampled[s.index,])
    return(ll)
  }
  return(llikelihood.f)
}

