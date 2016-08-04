#Source file for fluEvidenceSynthesis with access to vaccine coverage options. Adapted from code by Baguelin 2013 and Edwin van Leevvan, by TajWenzel.
#riskratio=risk.ratios.null

#load data in
vstrategy<-function(riskratio)
  {
data<-data(coverage)
cov <- coverage[,c("at.risk.under.65","X65")]/100.0

  #efficacy
  efficacy <- rep(c(0.7, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3), nrow(riskratio))
  dates <-coverage$Date

          #set empty storage matrix on the dimension agegroups*time data=nrows, ncol=efficacy
          calendar <- matrix(rep(0),nrow=length(dates),ncol = length(efficacy))


############SETTING UP VACCINE PROGRAM, DEFINE WHICH AGE GROUPS ARE VACCINATED 
# Set rate of vaccine uptake for different dates/age groups

t=1:length(dates);
calendar[t,c(1)] <- 0 #<6 months old
calendar[t,c(2)] <-coverage[,c("Under.65")]/100.0 #2-4 year olds
calendar[t,c(3)] <-coverage[,c("Under.65")]/100.0 #5-10, inclusive
calendar[t,c(4)] <-coverage[,c("Under.65")]/100.0 #11-15
calendar[t,c(5)] <-0 #16-24
calendar[t,c(6,7)] <-0 #adults 25+, 45+ and seniors

        #compile vaccine schedule
vaccine3 <- as.vaccination.calendar(efficacy = efficacy, dates = dates, 
                                    coverage = calendar, no_risk_groups=NULL)

v.output<-list(vaccine3$efficacy,vaccine3$calendar)
return(v.output)
}