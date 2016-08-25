#Source file for fluEvidenceSynthesis with access to vaccine coverage options. Adapted from code by Baguelin 2013 and Edwin van Leevvan, by TajWenzel.
#riskratio=risk.ratios.null

#load data in
vstrategy<-function(riskratio,scenario)
  {
data<-data(coverage)
cov <- coverage[,c("Under.65","X65","at.risk.under.65","X65")]/100.0


  #efficacy
  efficacy <- rep(c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3), nrow(riskratio))
  dates <-coverage$Date

          #set empty storage matrix on the dimension agegroups*time data=nrows, ncol=efficacy
          calendar <- matrix(rep(0),nrow=length(dates),ncol = length(efficacy))


############SETTING UP VACCINE PROGRAM, DEFINE WHICH AGE GROUPS ARE VACCINATED 
# Set rate of vaccine uptake for different dates/age groups
          #additional analysis for amount of QALY's potentially lost due to mis-reporting from practices not responding
          #(1,5,12,17,26,46,65)
          ####0) STATUS QUO is vaccination of elderly (65+) and high risk proportion (ages 0-85)
          if (scenario==0)
          {
            t=1:length(dates);
            calendar[t,c(1)] <- 0 #<6 months old
            calendar[t,c(2)] <-coverage[,c("at.risk.under.65")]/100.0 #1-4 year olds
            calendar[t,c(3)] <-coverage[,c("at.risk.under.65")]/100.0 #5-11, inclusive
            calendar[t,c(4)] <-coverage[,c("at.risk.under.65")]/100.0 #11-16
            calendar[t,c(5)] <-coverage[,c("at.risk.under.65")]/100.0 #17-25
            calendar[t,c(6)] <-coverage[,c("at.risk.under.65")]/100.0 #adults 26-45 
            calendar[t,c(7)] <-coverage[,c("at.risk.under.65")]/100.0 #adults 46-65 
            calendar[t,c(8)] <-coverage[,c("X65")]/100.0 #seniors 65+
          }
            
####1) Preschool Vaccination only; ages 2-4     
          if(scenario==1)
          {
            t=1:length(dates);
            calendar[t,c(1)] <- 0 #<6 months old
            calendar[t,c(2)] <-coverage[,c("Under.65")]/100.0 #1-4 year olds
            calendar[t,c(3)] <-0 #5-11, inclusive
            calendar[t,c(4)] <-0 #11-16
            calendar[t,c(5)] <-0 #17-25
            calendar[t,c(6)] <-0 #adults 26-45 
            calendar[t,c(7)] <-0 #adults 46-65 
            calendar[t,c(8)] <-0 #seniors 65+
          }
          
####2) Primary School Vaccination only; ages 5-11            
            if(scenario==2)    
              {
                t=1:length(dates);
                calendar[t,c(1)] <-0 #<6 months old
                calendar[t,c(2)] <-0 #1-4 year olds
                calendar[t,c(3)] <-coverage[,c("Under.65")]/100.0 #5-11
                calendar[t,c(4)] <-0 #12-16
                calendar[t,c(5)] <-0 #17-25
                calendar[t,c(6)] <-0 #adults 26-45 
                calendar[t,c(7)] <-0 #adults 46-65 
                calendar[t,c(8)] <-0 #seniors 65+
            }
####3) Secondary school vaccination only; ages 12-16        
              if(scenario==3)
              {
                t=1:length(dates);
                calendar[t,c(1)] <-0 #<6 months old
                calendar[t,c(2)] <-0 #1-4 year olds
                calendar[t,c(3)] <-0 #5-11, inclusive
                calendar[t,c(4)] <-coverage[,c("Under.65")]/100.0 #11-16
                calendar[t,c(5)] <-0 #17-25
                calendar[t,c(6)] <-0 #adults 26-45 
                calendar[t,c(7)] <-0 #adults 46-65 
                calendar[t,c(8)] <-0 #seniors 65+
              }

####4) Preschool+Primary School Vaccination; ages 2-11          
            if(scenario==4)
               {
                  t=1:length(dates);
                  calendar[t,c(1)] <-0 #<6 months old
                  calendar[t,c(2)] <-coverage[,c("Under.65")]/100.0 #1-4 year olds
                  calendar[t,c(3)] <-coverage[,c("Under.65")]/100.0 #5-11, inclusive
                  calendar[t,c(4)] <-0 #11-16
                  calendar[t,c(5)] <-0 #17-25
                  calendar[t,c(6)] <-0 #adults 26-45 
                  calendar[t,c(7)] <-0 #adults 46-65 
                  calendar[t,c(8)] <-0 #seniors 65+
            }
   
####5) Preschool+Primary+Secondary school vaccination; ages 2-16 
                  if(scenario==5)
                  {
                    t=1:length(dates);
                    calendar[t,c(1)] <-0 #<6 months old
                    calendar[t,c(2)] <-coverage[,c("Under.65")]/100.0 #1-4 year olds
                    calendar[t,c(3)] <-coverage[,c("Under.65")]/100.0 #5-11, inclusive
                    calendar[t,c(4)] <-coverage[,c("Under.65")]/100.0 #12-16
                    calendar[t,c(5)] <-0 #17-25
                    calendar[t,c(6)] <-0 #adults 26-45 
                    calendar[t,c(7)] <-0 #adults 46-65 
                    calendar[t,c(8)] <-0 #seniors 65+
                  }
 
####6) Preschool+Secondary school vaccination; ages 2-4, 12-16                   
                if(scenario==6)
                    {
                      t=1:length(dates);
                      calendar[t,c(1)] <-0 #<6 months old
                      calendar[t,c(2)] <-coverage[,c("Under.65")]/100.0 #1-4 year olds
                      calendar[t,c(3)] <-0 #5-11, inclusive
                      calendar[t,c(4)] <-coverage[,c("Under.65")]/100.0 #11-16
                      calendar[t,c(5)] <-0 #17-25
                      calendar[t,c(6)] <-0 #adults 26-45 
                      calendar[t,c(7)] <-0 #adults 46-65 
                      calendar[t,c(8)] <-0 #seniors 65+
                }
          
####7) Primary+Secondary school vaccination; ages 4-16
             if(scenario==7)
             {
               t=1:length(dates);
               calendar[t,c(1)] <-0 #<6 months old
               calendar[t,c(2)] <- #1-4 year olds
               calendar[t,c(3)] <-coverage[,c("Under.65")]/100.0 #5-11, inclusive
               calendar[t,c(4)] <-coverage[,c("Under.65")]/100.0 #11-16
               calendar[t,c(5)] <-0 #17-25
               calendar[t,c(6)] <-0 #adults 26-45 
               calendar[t,c(7)] <-0 #adults 46-65 
               calendar[t,c(8)] <-0 #seniors 65+
             }
        #compile vaccine schedule
                   
vaccine3 <- as.vaccination.calendar(efficacy = efficacy, dates = dates, 
                                    coverage = calendar)

v.output<-list(vaccine3$efficacy,vaccine3$calendar)
return(v.output)
}