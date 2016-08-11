#rm(list = ls())
setwd("/Users/Natasha/Dropbox/UKfluworkGIT")
library(fluEvidenceSynthesis)
library(ggplot2)     # base plots are for Coursera professors
library(scales)      # pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   # a helper for arranging individual ggplot objects
library(ggthemes)    # has a clean theme for ggplot2
library(viridis)     # best. color. palette. evar.
library(knitr)       # kable : prettier data.frame output
library(data.table)


#install.packages('data.table', type='source', repos='https://Rdatatable.github.io/data.table', force=TRUE)

agegrouplimits<-age.group.limits

#Reconstruct contract matricies from posterior contact data 
temp = list.files(pattern="*table.csv")
allpolymod = lapply(temp, fread)
fname<-c('BEid.csv','DEid.csv', 'FIid.csv', 'GBid.csv', 'ITid.csv', 'LUid.csv', 'NLid.csv', 'PLid.csv')


#i.top<-1
#for(i.top in 1:(length(allpolymod)-1))
for(i.top in 1:length(fname))
  {
contact.ids<-fread(fname[i.top])

out<-10000
data("age_sizes")
polymod<-allpolymod[[i.top]]; #country levels BE DE FI GB IT LU NL PL
#contact.ids<-allcontact.ids[[i.top]];

p.contacts<-NULL
for (jj in 1:out)
{
  p.contacts[[jj]]<-contact.matrix(as.matrix(polymod[contact.ids[[jj]],]),
                             age_sizes[,1], agegrouplimits)
}

start.matrix<-contact.matrix(as.matrix(polymod[contact.ids[[jj]],]),
                             age_sizes[,1], agegrouplimits)

trial<-array(unlist(p.contacts), c(7,7,out))

set.mat<-matrix(
            data=c(1, 1, 1, 1, 1, 1, 1,
                   0, 1, 1, 1, 1, 1, 1,
                   0, 0, 1, 1, 1, 1, 1,
                   0, 0, 0, 1, 1, 1, 1,
                   0, 0, 0, 0, 1, 1, 1,
                   0, 0, 0, 0, 0, 1, 1,
                   0, 0, 0, 0, 0, 0, 1), nrow=7, ncol=7, byrow=TRUE)

  iter<-as.data.frame(which(set.mat==1, TRUE));
  
  #for (qq in 1:length(iter$row))
  #{
   #   hist(trial[as.numeric(iter[qq,1]),as.numeric(iter[qq,2]),],breaks=20,col=3)
  #}
  
  se.difference<-NULL
  se.sum<-NULL
  se.sum[[1]]<-matrix(rep(0,49),nrow=7,byrow=TRUE);
  for (ff in 2:out)
   {
  se.difference[[ff]]<-sqrt((abs(trial[,,ff]-trial[,,1]))^2)
  se.sum[[ff]]<-se.difference[[ff]]+se.sum[[ff-1]]
   }
  
  
  sd.set<-se.sum[[dim(trial)[3]]]/(dim(trial)[3])
  sd.set[lower.tri(sd.set,diag=FALSE)]<-NA
  #se.difference<-array(unlist(se.difference), c(7,7,out))

##################################################################
# Heatmaps for each country
##############################################################
  library(reshape2)
  #convert to molten data frame
  sd.4.ggplot <- melt(sd.set)
  
  #names
  sell<-c('<1','1-4','5-15','15-24','25-44','45-64','65+')
  gname<-c('BE Posterior', 'DE Posterior', 'FI Posterior', 'GB Posterior', 'IT Posterior', 'LU Posterior', 'NL Posterior', 'PL Posterior')
  heatname<-c('BE', 'DE', 'FI', 'GB', 'IT', 'LU', 'NL', 'PL')
  
  #paste('heatmap_',heatname[i.top], “.png”, sep=””)
  #i.top<-1
  #saving information
  heatpath <- file.path("/Users","Natasha","Dropbox","UKfluwork", paste('heatmap',heatname[i.top], '.png', sep=''))
  
  png(file=heatpath, width=1000,height=1000, bg='white')
heatmap<-ggplot(data = sd.4.ggplot, aes(x=Var1,y=Var2,fill=value))+labs(x="Contacts", y="Participants", title=gname[i.top])+ geom_tile(color="white", size=0.1)+coord_equal()+scale_fill_viridis(name="Standard Deviation", label=comma)+ theme(axis.text=element_text(size=7))+theme(plot.title=element_text(hjust=0)) 
  
print(heatmap)
    dev.off()

    #init.contacts<- contact.matrix(as.matrix(polymod[current.contact.ids,]),
    #age_sizes[,1], agegrouplimits)
###########################################################################
# Prior and Posterior kernal density plot    
#############################################################################
  #temp = list.files(pattern="*table.csv")
  #allpolymod = lapply(temp, read.csv)
  #polymod<-allpolymod[[4]]; #country levels BE DE FI GB IT LU NL PL

  proposed.contact.ids <- seq(1,nrow(polymod))
  
n<-round(runif(out,1,length(proposed.contact.ids)))
boot.ids<<-list()
for (yy in 1:length(n))
{ #actual id numbers of contacts
    id.subtract <<- round(runif(n[yy],1,length(proposed.contact.ids)))
    id.replace <<-round(runif(n[yy],1,length(proposed.contact.ids)))
    proposed.contact.ids[id.subtract]<-id.replace
    
    sampled.ids<- contact.matrix(as.matrix(polymod[proposed.contact.ids,]),
                              age_sizes[,1], agegrouplimits)
    
    boot.ids[[yy]]<-as.matrix(sampled.ids);
    }


reup<-array(unlist(boot.ids), c(7,7,out))

fwrite(as.data.frame(reup),paste('resample',gname[i.top],'.csv'))


#head(reup)
#iter<-as.data.frame(which(set.mat==1, TRUE));
#dev.off()

#graph order
graph.mat<-matrix(
  data=c(1, 2, 4, 7, 11, 16, 22,
         0, 3, 5, 8, 12, 17, 23,
         0, 0, 6, 9, 13, 18, 24,
         0, 0, 0,10, 14, 19, 25,
         0, 0, 0, 0, 15, 20, 26,
         0, 0, 0, 0, 0,  21, 27,
         0, 0, 0, 0, 0, 0,  28), nrow=7, ncol=7, byrow=TRUE)


#names
ppname<-c('BE', 'DE', 'FI', 'GB', 'IT', 'LU', 'NL', 'PL')

#saves
#paste('heatmap_',heatname[i.top], “.png”, sep=””)
pppath <- file.path("/Users","Natasha","Dropbox","UKfluwork", paste('prepost',ppname[i.top], '.png', sep=''))

dev.new()
png(filename=pppath,width=1500,height=1500)

layout(graph.mat)
par(mar=c(2,2,2,2))
for (qq in 1:length(iter$row))
 {
xmin<-min(trial[iter[qq,1],iter[qq,2],],reup[iter[qq,1],iter[qq,2],])
xmax<-max(trial[iter[qq,1],iter[qq,2],],reup[iter[qq,1],iter[qq,2],])

prior.col <- c('tomato') 
prior.trans <- adjustcolor(prior.col, alpha.f = 0.4) 
post.col<-c('lightskyblue')
post.trans<-adjustcolor(post.col,alpha.f=0.6)

  g<-density(reup[iter[qq,1],iter[qq,2],])
  d<-density(trial[iter[qq,1],iter[qq,2],], bw=g$bw)

  ymax<-max(d$y,g$y)  
  ymin<-min(d$y,g$y)
  plot(g,type='n',xlim=c(xmin, xmax),ylim=c(ymin,ymax),main=c(sell[iter[qq,1]],sell[iter[qq,2]]))
  polygon(g,col=prior.trans,border='gray')
  line(d)
  polygon(d, col=post.trans, border="gray")
 }#end graph loop for matrix
  dev.off()
}#end overall country loop