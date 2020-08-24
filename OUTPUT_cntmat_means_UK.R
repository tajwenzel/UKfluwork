library(fluEvidenceSynthesis)
library(ggplot2)     # base plots are for Coursera professors
library(scales)      # pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   # a helper for arranging individual ggplot objects
library(ggthemes)    # has a clean theme for ggplot2
library(viridis)     # best. color. palette. evar.
library(knitr)       # kable : prettier data.frame output
library(data.table)  # MOST IMPORTANT PACKAGE SO THIS CODE ISN'T SLOW
library(abind)
library(xtable)
library(plyr)

##----start sampling prior and posterior contact matrix distribution

prior.posterior.cnt.matrix.avg<-function(country, strain, n.samp, version, ukstatus)
{
  post.means<-array()
  prior.means<-array()
  
  p1<-fread(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Inputdata', paste(allpolymod[i.country])))
  polymod<-cbind(p1,p1$V12) #add column for 75+ year olds
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
  
  if(ukstatus==T)
  {
  ct.samp<-list.files(pattern=glob2rx(paste0('UK','PostContactSample',version, strain.name[strain])))
  flu.samp<-list.files(pattern=glob2rx(paste0('UK', 'PostSample',version, strain.name[strain])))
  }else{   
  ct.samp<-list.files(pattern=glob2rx(paste0(fname[i.country],'PostContactSample',version, strain.name[strain])))
  flu.samp<-list.files(pattern=glob2rx(paste0(fname[i.country], 'PostSample',version, strain.name[strain])))}

  
if(strain==1){year.length<-length(ct.samp)}
if(strain==2){year.length<-length(ct.samp)}
if(strain==3){year.length<-length(ct.samp)}

for(i.season in 1:year.length)
{
  
  ####POSTERIOR
  load(ct.samp); 
  load(flu.samp)
  #if(strain==2) {load(ct.tablesH3[[i.season]]); load(flu.tablesH3[[i.season]])}
  #if(strain==3) {load(ct.tablesB[[i.season]]); load(flu.tablesB[[i.season]])}
  
  #mcmc.post<-mcmc.result$batch #rename mcmc output to simple data.table
  
  #iteration<-c(round(runif(n.samp,dim(mcmc.post)[1]*0.80,dim(mcmc.post)[1]))) 
  #take only last 20% of the runs to use for posterior distribution for calculations
  #rand.contact.samp[[i.season]]<-contact.ids[iteration]
  rand.contact.ids<-rand.contact.samp[[i.season]]
  
  #library(ggplot2)
  #colnames(mcmc.result$batch) <- c("eps1", "eps2", "eps3", "psi", "q",
   #                                "susc1", "susc2", "susc3", "I0")
  
  
  p.contacts<-NULL
  ct.function<-function(x) contact_matrix(as.matrix(polymod[x,]),age_sizes[,1], age.group.limits)
  poly.data<-lapply(rand.contact.ids,ct.function) #indexed list
  
  check<-array(unlist(poly.data), c(11,11,length(rand.contact.ids)))
  
  ifelse(i.season==1, post.means<-as.array(check), post.means<-abind(post.means, as.matrix(check), along=3))
  #####PRIOR 
  
  setwd(file.path("/Users","Natasha","Dropbox","UKfluworkGIT", 'Outputdata', paste(fname[i.country])));
  
  #load('GBPriorSample1H1N1')
  #load('GBsampled.prior')
  #load('prior.meansGBH1N1')
  
  proposed.contact.ids <- seq(1,nrow(polymod),1)
  n<-round(runif(n.samp,1,length(proposed.contact.ids)))
  
  for (yy in 1:length(n.samp))
  { 
    #actual id numbers of contacts
    #id.subtract <- round(n[yy],1,length(proposed.contact.ids)))
    id.subtract <- n[yy]
    id.replace <-round(runif(1,1,length(proposed.contact.ids)))
    proposed.contact.ids[id.subtract]<-id.replace
    
    sampled.ids<-lapply(proposed.contact.ids,ct.function)
  }
  
  #check2<-apply(array(unlist(sampled.ids), c(11,11,n.samp)), c(1,2), mean)
  check2<-array(unlist(sampled.ids), c(11,11,n.samp))
  
  ifelse(i.season==1, prior.means<-as.array(check2), prior.means<-abind(prior.means, as.matrix(check2), along=3))
} 
  
#These outputs are in lists so DO NOT save them as CSV or you will have to do the list/unlist stuff
save(prior.means,file=paste0('prior.means',fname[i.country], strain.name[strain]))
#save(prior.eigen,file=paste0('prior.eigen',fname[i.country], strain.name[strain]))
save(post.means,file=paste0('post.means',fname[i.country], strain.name[strain]))
}
