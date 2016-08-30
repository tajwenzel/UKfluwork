#Vstrategy with strain differences


library(reshape2)
library(tidyr)
library(dplyr)

#transform date list into first column of calendar list. Append calendar lists to each other. Note that dates here are untransformed in the saved files and they are not sorted by date. They need to be re-transformed and then sorted to avoid confusion (as excel cannot read the dates in that format.)

for(s in 1:3) 
{
     for (f in 1:length(vaccine.calendars))
      {
       if(f==1)
          {
        supple<-cbind(as.Date(vaccine.calendars[[f]][[s]]$dates,origin="1970-01-01"),
                      vaccine.calendars[[f]][[s]]$calendar)
        sup.eff<-matrix(c(rep(vaccine.calendars[[f]][[s]]$efficacy,
                              length(vaccine.calendars[[f]][[s]]$dates))),
                        nrow=length(vaccine.calendars[[f]][[s]]$dates),byrow=TRUE)
        vset<<-NULL
        vset[[f]]<-as.data.frame(cbind(supple,sup.eff))
        } 
       else
        {
        supple2<-cbind(as.Date(vaccine.calendars[[f]][[s]]$dates,origin="1970-01-01"),
                       vaccine.calendars[[f]][[s]]$calendar)
        sup.eff2<-matrix(c(rep(vaccine.calendars[[f]][[s]]$efficacy,length(vaccine.calendars[[f]][[s]]$dates))),nrow=length(vaccine.calendars[[f]][[s]]$dates),byrow=TRUE)
        colnames(supple2)<-NULL
        colnames(sup.eff2)<-NULL
        vset[[f]]<-as.data.frame(cbind(supple2,sup.eff2))
        }
     }
  
    unames<-c('H1','H3','B')
    save(vset,file=paste0('cov.eff',unames[s]))
}

load('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/cov.effH1'); cov.effH1<-vset;
load('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/cov.effH3'); cov.effH3<-vset;
load('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/cov.effB'); cov.effB<-vset;
#coverage data will have one additional column compared to efficacy as dates column was appended to coverage. We keep the coverage and efficacy datasets together so that we can sort by dates later on and have the two categories be aligned.
type[order(as.Date(type$V1))] #order with increasing date

as.Date(cov.effH1[[1]]$V1,origin="1970-01-01")
cov.effH1[order(as.Date(cov.effH1[[1]]$V1,origin="1970-01-01"))]

as.Date(type$V1,origin="1970-01-01")
type<-fread(input='cov.effB.csv')
type$V1<-as.Date(type$V1,origin="1970-01-01")
coverageB<-type[order(as.Date(type$V1))] #order with increasing date

#now cut and separate for coverage
vaccine.data<-NULL
vaccine.data$coverage<-coverageB[1:22]


as.Date(vaccine.calendars[[1]][[1]]$dates,origin="1970-01-01")
blue<-as.Date(vaccine.calendars[[11]][[2]]$dates,origin="1970-01-01")
set<-cbind(blue, vaccine.calendars[[1]][[2]]$calendar)
colnames(set)<-NULL

############################################################################
#Pull actual coverage from rates
###########################################################################

tx_start=as.Date(c("1995-09-01","1996-09-01","1997-09-01","1998-09-01","1999-09-01","2000-09-01",
                "2001-09-01","2002-09-01","2003-09-01","2004-09-01","2005-09-01","2006-09-01",
                "2007-09-01","2008-09-01","2009-09-01","2010-09-01","2011-09-01","2012-09-01",
                "2013-09-01"),origin="1970-01-01")

#load covrage rates and efficacy
load('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/cov.effH1'); cov.effH1<-vset;
load('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/cov.effH3'); cov.effH3<-vset;
load('/Users/Natasha/Dropbox/UKfluworkGIT/Inputdata/cov.effB'); cov.effB<-vset;

#getting actual coverage
for(s in 1:3)
{
  strain<-list(cov.effH1,cov.effH3,cov.effB)
  cov.eff<-strain[[s]]
for(kaw in 1:length(tx_start))
 {
filler<-c(rep(0,21),as.vector(cov.eff[[kaw]][1,23:43],mode='numeric'))
start=as.numeric(tx_start[kaw])
overj<-data.frame(diff(as.matrix(rbind(c(start,filler),cov.effB[[kaw]]))))
overl<-cumsum(abs(overj$V1*overj[,2:22]))
overl[overl>1]<-1
cov.eff[[kaw]][,2:22]<-overl
cov.eff[[kaw]]<-rbind(c(start,filler),cov.eff[[kaw]])

  unames<-c('H1','H3','B')
save(cov.eff,file=paste0('coverage',unames[s]))
  }
}

load('/Users/Natasha/Dropbox/UKfluworkGIT/coverageH1'); coverageH1<-cov.eff;
load('/Users/Natasha/Dropbox/UKfluworkGIT/coverageH3'); coverageH3<-cov.eff;
load('/Users/Natasha/Dropbox/UKfluworkGIT/coverageB'); coverageB<-cov.eff;
cov.eff<-NULL




first<-as.numeric((set[1]-tx_start[13]))*cov.effB[[13]][1,2:22]
add.cov<-first+abs(as.numeric((set[1]-set[2])))*cov.effB[[13]][2,2:22]
add.cov<-add.cov+abs(as.numeric((set[1]-set[2])))*cov.effB[[13]][2,2:22]

cov.effH3
cov.effB
starting_year<-1970
c(as.Date(paste0(starting_year, 
                 "-10-01")), as.Date(paste0(starting_year, "-11-01")), 
  as.Date(paste0(starting_year, "-12-01")), as.Date(paste0(starting_year + 1, "-01-01")), as.Date(paste0(starting_year + 1, "-02-01")))

vaccine.calendars$H3
vaccine.calendars$H1
vaccine.calendars$B