#Source file for fluEvidenceSynthesis with access to vaccine coverage options. Adapted from code by Baguelin 2013 and Edwin van Leevvan, by TajWenzel.
#riskratio=risk.ratios.ce
#scenario=1
#strain=2
#cov.eff<-coverageB
#season=4

vstrategy<-function(riskratio,scenario,cov.eff,season)
  {
  
   cov <- c(cov.eff[[season]][,2:22])
    #USE COV INSTEAD OF COVERAGE SO DON'T HAVE TO DO THE /100 at every step

##################################################################################################
#Efficacy
##################################################################################################
    #eff.pull<-cov.eff[[4]][1,23:43]
    eff.pull<-as.vector(cov.eff[[season]][1,23:43])
          #pull first row of efficacies from matrix, matrix starts at 23 as 21 risk*age groups in                  original data +1 column for dates. This encompasses 2 risk*age groups (14 columns)
          
    #efficacy <-eff.pull;
    dates <-as.Date(cov.eff[[season]]$V1,origin="1970-01-01");
     

          #set empty storage matrix on the dimension agegroups*time data=nrows, ncol=efficacy
          calendar <- matrix(rep(0),nrow=length(dates),ncol = length(eff.pull))


############SETTING UP VACCINE PROGRAM, DEFINE WHICH AGE GROUPS ARE VACCINATED 
# Set rate of vaccine uptake for different dates/age groups
          #additional analysis for amount of QALY's potentially lost due to mis-reporting from practices not responding
          ####0) STATUS QUO is vaccination of elderly (65+) and high risk proportion (ages 0-85)
          if (scenario==1)
          {
            t=1:length(dates);
            calendar[t,c(1)] <-cov[[1]] #0-1
            calendar[t,c(2)] <-cov[[2]] #1-4
            calendar[t,c(3)] <-cov[[3]] #5-11
            calendar[t,c(4)] <-cov[[3]] #12-14
            calendar[t,c(5)] <-cov[[3]] #15-16
            calendar[t,c(6)] <-cov[[4]] #17-24
            calendar[t,c(7)] <-cov[[5]] #25-44
            calendar[t,c(8)] <-cov[[6]] #45-64
            calendar[t,c(9)] <-cov[[7]] #65-74
            calendar[t,c(10)] <-cov[[7]] #75+
            calendar[t,c(11)] <-0 
            calendar[t,c(12)] <-0 
            calendar[t,c(13)] <-0 
            calendar[t,c(14)] <-0 
            calendar[t,c(15)] <-0 
            calendar[t,c(16)] <-0
            calendar[t,c(17)] <-0
            calendar[t,c(18)] <-0
            calendar[t,c(19)] <-cov[[14]]
            calendar[t,c(20)] <-cov[[14]] 
          }
            
####1) Preschool Vaccination only; ages 2-4     
          if(scenario==2)
          {
            t=1:length(dates);
            calendar[t,c(1)] <- 0 #<6 months old
            calendar[t,c(2)] <-coverage[,c("at.risk.under.65")]/100.0 #1-4 year olds
            calendar[t,c(3)] <-0 #5-11, inclusive
            calendar[t,c(4)] <-0 #11-16
            calendar[t,c(5)] <-0 #17-25
            calendar[t,c(6)] <-0 #adults 26-45 
            calendar[t,c(7)] <-0 #adults 46-65 
            calendar[t,c(8)] <-0 #seniors 65+
            calendar[t,c(9)] <- 0 #<6 months old
            calendar[t,c(10)] <-coverage[,c("Under.65")]/100.0 #1-4 year olds #1-4 year olds
            calendar[t,c(11)] <-0 #5-11, inclusive
            calendar[t,c(12)] <-0 #11-16
            calendar[t,c(13)] <-0 #17-25
            calendar[t,c(14)] <-0 #adults 26-45 
            calendar[t,c(15)] <-0 #adults 46-65 
            calendar[t,c(16)] <-0
          }
          
####2) Primary School Vaccination only; ages 5-11            
            if(scenario==3)    
              {
                t=1:length(dates);
                calendar[t,c(1)] <-0 #<6 months old
                calendar[t,c(2)] <-0 #1-4 year olds
                calendar[t,c(3)] <-coverage[,c("at.risk.under.65")]/100.0 #5-11
                calendar[t,c(4)] <-0 #11-16
                calendar[t,c(5)] <-0 #17-25
                calendar[t,c(6)] <-0 #adults 26-45 
                calendar[t,c(7)] <-0 #adults 46-65 
                calendar[t,c(8)] <-0 #seniors 65+
                calendar[t,c(9)] <- 0 #<6 months old
                calendar[t,c(10)] <-0 #1-4 year olds
                calendar[t,c(11)] <-coverage[,c("Under.65")]/100.0 #5-11
                calendar[t,c(12)] <-0 #11-16
                calendar[t,c(13)] <-0 #17-25
                calendar[t,c(14)] <-0 #adults 26-45 
                calendar[t,c(15)] <-0 #adults 46-65 
                calendar[t,c(16)] <-0
            }
####3) Secondary school vaccination only; ages 12-16        
              if(scenario==4)
              {
                t=1:length(dates);
                calendar[t,c(1)] <-0 #<6 months old
                calendar[t,c(2)] <-0 #1-4 year olds
                calendar[t,c(3)] <-0 #5-11, inclusive
                calendar[t,c(4)] <-coverage[,c("at.risk.under.65")]/100.0 #11-16
                calendar[t,c(5)] <-0 #17-25
                calendar[t,c(6)] <-0 #adults 26-45 
                calendar[t,c(7)] <-0 #adults 46-65 
                calendar[t,c(8)] <-0 #seniors 65+
                calendar[t,c(9)] <- 0 #<6 months old
                calendar[t,c(10)] <-0 #1-4 year olds
                calendar[t,c(11)] <-0 #5-11, inclusive
                calendar[t,c(12)] <-coverage[,c("Under.65")]/100.0 #11-16
                calendar[t,c(13)] <-0 #17-25
                calendar[t,c(14)] <-0 #adults 26-45 
                calendar[t,c(15)] <-0 #adults 46-65 
                calendar[t,c(16)] <-0
              }

####4) Preschool+Primary School Vaccination; ages 2-11          
            if(scenario==5)
               {
                  t=1:length(dates);
                  calendar[t,c(1)] <-0 #<6 months old
                  calendar[t,c(2)] <-coverage[,c("at.risk.under.65")]/100.0 #1-4 year olds
                  calendar[t,c(3)] <-coverage[,c("at.risk.under.65")]/100.0 #5-11, inclusive
                  calendar[t,c(4)] <-0 #11-16
                  calendar[t,c(5)] <-0 #17-25
                  calendar[t,c(6)] <-0 #adults 26-45 
                  calendar[t,c(7)] <-0 #adults 46-65 
                  calendar[t,c(8)] <-0 #seniors 65+
                  calendar[t,c(9)] <- 0 #<6 months old
                  calendar[t,c(10)] <-coverage[,c("Under.65")]/100.0 #1-4 year olds
                  calendar[t,c(11)] <-coverage[,c("Under.65")]/100.0 #5-11, inclusive
                  calendar[t,c(12)] <-0 #11-16
                  calendar[t,c(13)] <-0 #17-25
                  calendar[t,c(14)] <-0 #adults 26-45 
                  calendar[t,c(15)] <-0 #adults 46-65 
                  calendar[t,c(16)] <-0
            }
   
####5) Preschool+Primary+Secondary school vaccination; ages 2-16 
                if(scenario==6)
                  {
                    t=1:length(dates);
                    calendar[t,c(1)] <-0 #<6 months old
                    calendar[t,c(2)] <-coverage[,c('at.risk.under.65')]/100.0 #1-4 year olds
                    calendar[t,c(3)] <-coverage[,c('at.risk.under.65')]/100.0 #5-11, inclusive
                    calendar[t,c(4)] <-coverage[,c('at.risk.under.65')]/100.0 #12-16
                    calendar[t,c(5)] <-0 #17-25
                    calendar[t,c(6)] <-0 #adults 26-45 
                    calendar[t,c(7)] <-0 #adults 46-65 
                    calendar[t,c(8)] <-0 #seniors 65+
                    calendar[t,c(9)] <- 0 #<6 months old
                    calendar[t,c(10)] <-coverage[,c('Under.65')]/100.0 #1-4 year olds
                    calendar[t,c(11)] <-coverage[,c('Under.65')]/100.0 #5-11, inclusive
                    calendar[t,c(12)] <-coverage[,c('Under.65')]/100.0 #12-16
                    calendar[t,c(13)] <-0 #17-25
                    calendar[t,c(14)] <-0 #adults 26-45 
                    calendar[t,c(15)] <-0 #adults 46-65 
                    calendar[t,c(16)] <-0
                  }
 
####6) Preschool+Secondary school vaccination; ages 2-4, 12-16                   
                if(scenario==7)
                    {
                      t=1:length(dates);
                      calendar[t,c(1)] <-0 #<6 months old
                      calendar[t,c(2)] <-coverage[,c("at.risk.under.65")]/100.0 #1-4 year olds
                      calendar[t,c(3)] <-0 #5-11, inclusive
                      calendar[t,c(4)] <-coverage[,c("at.risk.under.65")]/100.0 #11-16
                      calendar[t,c(5)] <-0 #17-25
                      calendar[t,c(6)] <-0 #adults 26-45 
                      calendar[t,c(7)] <-0 #adults 46-65 
                      calendar[t,c(8)] <-0 #seniors 65+
                      calendar[t,c(9)] <- 0 #<6 months old
                      calendar[t,c(10)] <-coverage[,c("Under.65")]/100.0 #1-4 year olds
                      calendar[t,c(11)] <-0 #5-11, inclusive
                      calendar[t,c(12)] <-coverage[,c("Under.65")]/100.0 #11-16
                      calendar[t,c(13)] <-0 #17-25
                      calendar[t,c(14)] <-0 #adults 26-45 
                      calendar[t,c(15)] <-0 #adults 46-65 
                      calendar[t,c(16)] <-0
                }
          
####7) Primary+Secondary school vaccination; ages 4-16
                    if(scenario==8)
                      {
               t=1:length(dates);
               calendar[t,c(1)] <-0 #<6 months old
               calendar[t,c(2)] <- #1-4 year olds
               calendar[t,c(3)] <-coverage[,c("at.risk.under.65")]/100.0 #5-11, inclusive
               calendar[t,c(4)] <-coverage[,c("at.risk.under.65")]/100.0 #11-16
               calendar[t,c(5)] <-0 #17-25
               calendar[t,c(6)] <-0 #adults 26-45 
               calendar[t,c(7)] <-0 #adults 46-65 
               calendar[t,c(8)] <-0 #seniors 65+
               calendar[t,c(9)] <- 0 #<6 months old
               calendar[t,c(10)] <-0 #1-4 year olds
               calendar[t,c(11)] <-coverage[,c("Under.65")]/100.0 #5-11, inclusive
               calendar[t,c(12)] <-coverage[,c("Under.65")]/100.0 #11-16
               calendar[t,c(13)] <-0 #17-25
               calendar[t,c(14)] <-0 #adults 26-45 
               calendar[t,c(15)] <-0 #adults 46-65 
               calendar[t,c(16)] <-0 #seniors 65+
                      }
        #compile vaccine schedule
                   
vaccine3 <- as.vaccination.calendar(efficacy = eff.pull, dates = dates, 
                                    coverage = calendar,no_risk_groups=2,no_age_groups=10)

eff.out<-as.vector(unlist(vaccine3$efficacy),mode='numeric')
v.output<-list('efficacy'=eff.out,'dates'=vaccine3$dates,'calendar'=vaccine3$calendar)
return(v.output)
}