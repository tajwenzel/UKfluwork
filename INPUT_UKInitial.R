#initial parameters for the United Kingdom
function(season) 
  { 
  if(season >2) {start.pars<-c(0.01188150, 0.01831852, 0.05434378,
                               1.049317e-05, 0.1657944,
                               0.3855279, 0.9269811, 0.5710709,
                               -0.1543508 ); 
  return(start.pars)} 
  else {
start.pars<-c(0.01188150, 0.01831852, 0.05434378,
              1.049317e-05, 0.1657944,
              0.3855279, 0.9269811, 0.5710709,
              -0.1543508 );

return(start.pars)
  }
}

  
  
#c(0.01188150, 0.01831852, 0.05434378,1.049317e-05, 0.1657944,0.3855279, 0.9269811, 0.5710709,-0.1543508 )
#1-3: Ascertainment probabilty for three age groups (ϵiϵi)
#4: Outside infection (ψψ)
#5: Transmissibility (qq)
#6-8: Susceptibility for three age groups (σiσi)
#9: Initial number of infections (II), 10^I, assume how many people infected with new strain on September 1st, 20XX
#log_likelihood_cases(c(0.16611236, 0.04915033, 0.042402057), 1.337120e-02,ili.cases[s.index,],ili_monitored=ili.total[s.index,],positive[s.index,1:5], confirmed_samples=total.sampled[s.index,1:5])