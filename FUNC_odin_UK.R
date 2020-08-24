infection_odin_uk <- function(population, initial_infected, vschedule, contact_matrix, susceptibility, transmissibility, interval) 
{
  # Extract the date used from the vaccine calendar
  begin_date <- as.Date(paste0(format(vschedule$dates[1], "%Y"),"-09-01"))
  #if(any(format(vcalendar1$dates[1],  "%Y")==c('1996','2000','2004','2008','2012'))) #for leap year
  #{t <- as.numeric(seq(begin_date, begin_date + 366, interval))}else{
  #t <- as.numeric(seq(begin_date, begin_date + 365, interval))}
  t <- as.numeric(seq(begin_date, begin_date + interval*52, interval))
  no_groups <- 22
  no_risk_groups <- no_groups/nrow(contact_matrix)
  no_age_groups <- no_groups/no_risk_groups
  
  #calendar <- vaccine_calendar$calendar[c(nrow(vaccine_calendar$calendar),1:nrow(vaccine_calendar$calendar)),]
  #dates <- as.numeric(c(t[1], vaccine_calendar$dates))
  
  
  calendar <- vschedule$calendar[c(nrow(vschedule$calendar),1:nrow(vschedule$calendar)),]
  dates <- as.numeric(c(t[1], vschedule$dates))
  
  # Contacts matrix only covers one set of age groups, here we "repeat" it to also cover 
  # risk groups
  new_cij <- matrix(rep(0,no_groups*no_groups), nrow = no_groups)
  for (k in 1:no_risk_groups) {
    for (l in 1:no_risk_groups) {
      lk <- (k - 1)*no_age_groups + 1
      ll <- (l - 1)*no_age_groups + 1
      new_cij[lk:(lk + no_age_groups - 1), ll:(ll + no_age_groups - 1)] <- contact_matrix
    }
  }
  # Set the parameter values
  mod <- gen_seeiir_ag_vacc(no_groups = 22, cij = new_cij, trans = transmissibility,
                            pop = population[1:no_groups],
                            I0 = initial_infected[1:no_groups],
                            susc = rep(susceptibility,no_risk_groups),
                            alpha=c(rep(vschedule$efficacy[1],9),
                                    rep(vschedule$efficacy[7],2),
                                    rep(vschedule$efficacy[8],9),
                                    rep(vschedule$efficacy[14],2)),
                            dates = dates,
                            calendar = calendar[,1:no_groups],
                            gamma1 = 2/0.8, gamma2 = 2/1.8)
  
  output <- mod$run(t, hmax = NULL, method = "euler", hini = 0.25, atol = 1)
  
  
  cumulative.cases <- output[, (ncol(output) - 1*no_groups + 1):(ncol(output))]
  if(colnames(cumulative.cases)[1]!='cumI[1]') stop
  #Returning the differences in cumulative infections from one week to the other
  incident <- data.frame(cumulative.cases[2:(nrow(cumulative.cases)), ] - cumulative.cases[1:(nrow(cumulative.cases) - 1), ])
  
  #Cleanup and add Time column
  #colnames(y) <- names(population)
  mutate(Time = as.Date(t[1:nrow(incident)], origin = "1970-01-01"), incident)
}

gen_seeiir_ag_vacc <- odin::odin({
  # Number of groups
  no_groups <- user()
  
  # INITIAL CONDITIONS
  # Population size by age/risk group
  pop[] <- user()
  # Initial infection by age/risk group
  I0[] <- user()
  
  # MODEL PARAMETERS
  # Susceptibility
  susc[] <- user()
  
  # Transmissibility
  trans <- user()
  
  # Latent periods
  gamma1 <- user()
  gamma2 <- user()
  
  # Vaccine related variables 
  dates[] <- user()
  calendar[,] <- user()
  
  # efficacy
  alpha[] <- user()
  
  # Contact matrix
  cij[,] <- user()
  
  # Force of infection
  lambda[] <- trans * susc[i] * (sum(sij[i,]))
  
  # Vaccination. The rate is a step function that changes at each date according
  # to the passed calendar
  vI[] <- interpolate(dates, calendar, "constant")
  # Vaccination is given as a fraction vaccination, here we scale it to 
  # a rate
  sumN[] <- if (vI[i]>0) (S[i]+E1[i]+E2[i]+I1[i]+I2[i]+R[i]) else 0
  v[] <- if (sumN[i]>0) vI[i]*pop[i]/sumN[i] else 0
  
  # Transmission matrix
  sij[,] <- cij[i,j] * (I1[j] + I2[j] + I1v[j] + I2v[j])
  
  # Newly infected
  newInf[] <- lambda[i] * S[i]
  newInfv[] <- lambda[i] * Sv[i]
  
  # THE DERIVATIVES OF THE SEEIIR MODEL
  # Derivatives of the not vaccinated group
  deriv(S[]) <- -newInf[i] - v[i] * S[i]
  deriv(E1[]) <- newInf[i] - gamma1 * E1[i] - v[i] * E1[i]
  deriv(E2[]) <- gamma1 * (E1[i] - E2[i]) - v[i] * E2[i]
  deriv(I1[]) <- gamma1 * E2[i]  - gamma2 * I1[i] - v[i] * I1[i]
  deriv(I2[]) <- gamma2 * (I1[i] - I2[i]) - v[i] * I2[i]
  deriv(R[]) <- gamma2 * I2[i] - v[i] * R[i]
  
  # Derivatives vaccination group
  deriv(Sv[]) <- -newInfv[i] + v[i] * (1-alpha[i]) * S[i]
  deriv(E1v[]) <- newInfv[i] - gamma1 * E1v[i] + v[i] * E1[i]
  deriv(E2v[]) <- gamma1 * (E1v[i] - E2v[i]) + v[i] * E2[i]
  deriv(I1v[]) <- gamma1 * E2v[i]  - gamma2 * I1v[i] + v[i] * I1[i]
  deriv(I2v[]) <- gamma2 * (I1v[i] - I2v[i]) + v[i] * I2[i]
  deriv(Rv[]) <- gamma2 * I2v[i] + v[i] * (R[i] + alpha[i] * S[i])
  
  # Tracking the cumulative amount of infections over time for output of incidence
  deriv(cumI[]) <- newInf[i] + newInfv[i]
  
  # Initial value of the variables
  initial(S[1:no_groups]) <- pop[i] - I0[i]
  initial(E1[1:no_groups]) <- 0
  initial(E2[1:no_groups]) <- 0
  initial(I1[1:no_groups]) <- I0[i]
  initial(I2[1:no_groups]) <- 0
  initial(R[1:no_groups]) <- 0
  
  initial(Sv[1:no_groups]) <- 0
  initial(E1v[1:no_groups]) <- 0
  initial(E2v[1:no_groups]) <- 0
  initial(I1v[1:no_groups]) <- 0
  initial(I2v[1:no_groups]) <- 0
  initial(Rv[1:no_groups]) <- 0
  initial(cumI[1:no_groups]) <- 0
  
  # Set dimension of all variables/parameters
  dim(dates) <- user()
  dim(calendar) <- user()
  
  dim(pop) <- no_groups
  dim(I0) <- no_groups
  dim(susc) <- no_groups
  dim(lambda) <- no_groups
  dim(v) <- no_groups
  dim(vI) <- no_groups
  dim(sumN) <- no_groups  
  dim(alpha) <- no_groups
  dim(cij) <- c(no_groups, no_groups)
  dim(sij) <- c(no_groups, no_groups)
  
  dim(S) <- no_groups
  dim(E1) <- no_groups
  dim(E2) <- no_groups
  dim(I1) <- no_groups
  dim(I2) <- no_groups
  dim(R) <- no_groups
  dim(Sv) <- no_groups
  dim(E1v) <- no_groups
  dim(E2v) <- no_groups
  dim(I1v) <- no_groups
  dim(I2v) <- no_groups
  dim(Rv) <- no_groups
  dim(cumI) <- no_groups
  dim(newInf) <- no_groups
  dim(newInfv) <- no_groups
}, verbose = F)
