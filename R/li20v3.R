#' Covid spatPomp generator v3
#'
#' Generate a \sQuote{spatPomp} object for Covid in \code{U} cities in China.
#'
#' @param U A length-one numeric signifying the number of cities to be represented in the spatPomp object.
#' @param dt step size, in days, for the euler approximation.
#' @param sharedOneInterval estimated parameters that are equal for each unit and constant throughout time
#' @param sharedTwoInterval estimated parameters that are equal for each unit but change on day 15
#' @param version choose which set of parameters will be used for code testing
#' @param mob_modify_factor A multiplicative factor, which is used to control the magnitude of adjustment for the mobility data
#' @param days number of time points, in days
#' @param for_ibpf whether this object is generated for ibpf.
#' @import pomp
#' @import spatPomp
#' @return An object of class \sQuote{spatPomp} representing a \code{U}-dimensional spatially coupled Covid POMP model.
#'
#' @export
li20v3 <- function(U = 373, dt = 1/4,
		   sharedOneInterval = c("theta","tau","sigma_SE","E_0","Iu_0"),
                   sharedTwoInterval = c("alpha","Beta","mu","Z","D","Td"),
		   version = c("MLEperiod3","li20period1","li20period2","li20period3"),
		   mob_modify_factor = 20, days = 30, for_ibpf = T
) {
  if(version[1] == "MLEperiod3") {
    # mle from e3, loglik = -9200
    testPar <- c(
      alpha_be=0.1199543,
      Beta_be=0.8053392,
      alpha_af=0.4151188,
      Beta_af=0.3441618, 
      mu_be=0.9700681,
      Z_be=1.735166,
      D_be=4.635517,
      mu_af=0.7220044,
      Z_af=1.536596,
      D_af=2.714057,
      theta=1.752231,
      tau=0.3456199,
      sigma_SE=0.994764,
      Td_be=9,
      Td_af=6,
      E_0=934.2861,
      Iu_0=1435.945
    )
  } else if(version[1] == "li20period3") {
    testPar <- c(
      alpha_be=0.14,
      Beta_be=1.12,
      mu_be=0.55,
      Z_be=3.69,
      D_be=3.47,
      alpha_af=0.69,
      Beta_af=0.35,
      mu_af=0.43,
      Z_af=3.42,
      D_af=3.31,
      theta=1.36,
      tau=0.5,
      sigma_SE=0.5, 
      Td_be=9,
      Td_af=6,
      E_0=1000,
      Iu_0=1000
    )
  } else if(version[1] == "li20period2") {
    testPar <- c(
      alpha_be=0.14,
      Beta_be=1.12,
      mu_be=0.55,
      Z_be=3.69,
      D_be=3.47,
      alpha_af=0.65,
      Beta_af=0.52,
      mu_af=0.5,
      Z_af=3.6,
      D_af=3.14,
      theta=1.36,
      tau=0.5,
      sigma_SE=0.5, 
      Td_be=9,
      Td_af=6,
      E_0=1000,
      Iu_0=1000
    )
  } else if(version[1] == "li20period1") {
    testPar <- c(
      alpha_be=0.14,
      Beta_be=1.12,
      mu_be=0.55,
      Z_be=3.69,
      D_be=3.47,
      # period1 won't use _af, provide them to avoid an error
      alpha_af=0.69,
      Beta_af=0.35,
      mu_af=0.43,
      Z_af=3.42,
      D_af=3.31,
      theta=1.36,
      tau=0.5,
      sigma_SE=0.5, 
      Td_be=9,
      Td_af=6,
      E_0=1000,
      Iu_0=1000
    )
  }
  
 
  
  OneIntervalParNames <- c("theta","tau","sigma_SE","E_0","Iu_0")
  TwoIntervalParNames <- c("alpha","Beta","mu","Z","D","Td")
  sharedParNames <- c(
    if(length(sharedTwoInterval)>0) c(paste0(sharedTwoInterval,"_be"),paste0(sharedTwoInterval,"_af")) else NULL,
    sharedOneInterval)
  unitOneInterval <- setdiff(OneIntervalParNames,sharedOneInterval)
  unitTwoInterval <- setdiff(TwoIntervalParNames,sharedTwoInterval)
  unitParNames <- c(
    if(length(unitTwoInterval)>0) c(paste0(unitTwoInterval,"_be"),paste0(unitTwoInterval,"_af")) else NULL,
    unitOneInterval)
  if(length(unitParNames)==0) unitParNames <- NULL
  mobi <- metapoppkg::mobility+metapoppkg::v_by_g_day*mob_modify_factor
  popu <- metapoppkg::population[1:days+rep((0:372)*30,each=days),]
  incidence <- metapoppkg::incidence

  # choose the cities with largest number of cases
  index <- order(colSums(incidence)[-1], decreasing = T)[1:U]
  incidence <- incidence[1:days,c(1,index+1)]
  
  index_mob <- rep(index*14-13,each=14)+rep(0:13,U)
  mobi <- mobi[index_mob,index]
  popu <- popu[rep(index*days-days+1,each=days)+rep(0:(days-1),U),]

  to_C_array <- function(v)paste0("{",paste0(v,collapse=","),"}")
  mobi_rows <- apply(mobi,1,to_C_array)
  mobi_array <- to_C_array(mobi_rows)
  mobi_C <- Csnippet(paste0("const double mob[",U*14,"][",U,"] = ",mobi_array,"; "))

  covid_cases <- pomp::melt(incidence)[-c(1:days),]
  covid_covar <- as.data.frame(cbind(rep(incidence$Date[1:days],U),popu))
  colnames(covid_covar) <- c("date","city","pop")
  covid_covar$date <- as.numeric(covid_covar$date)
  covid_covar$pop <- as.integer(covid_covar$pop)
  colnames(covid_cases) <- c("city","date","cases")
  rownames(covid_covar) <- NULL
  
  if(for_ibpf == T) {
    set_fixed <- Csnippet(paste0("const int ", c(sharedOneInterval,sharedTwoInterval),
                                 "_expand = 1;\n", collapse=" "))
  } else {
    set_fixed <- Csnippet(paste0("const int ", c(sharedOneInterval,sharedTwoInterval),
                                 "_expand = 0;\n", collapse=" "))
  }
  set_expanded <- Csnippet(paste0("const int ", c(unitOneInterval, unitTwoInterval),
                                  "_expand = 1;\n", collapse=" "))
  if(!is.null(unitParNames) & !is.null(sharedParNames)){
    covid_globals <- Csnippet(
      paste(mobi_C, set_expanded, set_fixed, sep = "\n")
    )
  }
  if(is.null(unitParNames) & !is.null(sharedParNames)){
    covid_globals <- Csnippet(
      paste(mobi_C, set_fixed, sep = "\n")
    )
  }
  if(!is.null(unitParNames) & is.null(sharedParNames)){
    covid_globals <- Csnippet(
      paste(mobi_C, set_expanded, sep = "\n")
    )
  }

  covid_unit_statenames <- c('S','E','Iu','Ir','Ca','Cb','C')
  if(for_ibpf == T) {
    covid_paramnames <- c(
      if(!is.null(unitParNames)) paste0(rep(unitParNames,each=U),1:U) else NULL,
      if(!is.null(sharedParNames)) paste0(rep(sharedParNames,each=U),1:U) else NULL
    )
    covid_testPar <- c(
      if(!is.null(unitParNames)) rep(testPar[unitParNames],each=U) else NULL,
      if(!is.null(sharedParNames)) rep(testPar[sharedParNames],each=U) else NULL
    )
  } else {
    covid_paramnames <- c(
      if(!is.null(unitParNames)) paste0(rep(unitParNames,each=U),1:U) else NULL,
      if(!is.null(sharedParNames)) paste0(sharedParNames,"1") else NULL
    )
    covid_testPar <- c(
      if(!is.null(unitParNames)) rep(testPar[unitParNames],each=U) else NULL,
      if(!is.null(sharedParNames)) testPar[sharedParNames] else NULL
    )
  }
  names(covid_testPar) <- covid_paramnames

  covid_rprocess <- spatPomp_Csnippet(
    unit_statenames=c('S','E','Ir','Iu','Ca','Cb','C'),
    unit_paramnames=c("alpha_be","Beta_be","alpha_af","Beta_af","mu_be","Z_be",
      "D_be","mu_af","Z_af","D_af","theta","Td_be","Td_af","sigma_SE"),
    unit_covarnames=c('pop'),
    code='
      const int numTrans = 13; // number of transitions
      double rate[numTrans], dN[numTrans];
//      
//      |         |          |
//   TS |      TE |      TIu |
//      |         |          |
//      +    SE   +   EIu    +
//      S ------> E ------> Iu -----+ IuR
//      |         |  \       \ 
//   ST |      ET |   \       \ 
//      |         |    \       ---+ IuT
//      +         +     \        
//                       ---> Ir ----+ IrR
//                        \    
//                     EIr \    
//                          \    
//                           \
//                            \            CaCb       CbC
//                             -----+ Ca -------+ Cb ------+ C
//
      const int SE=0, ST=1;
      const int EIr=2, EIu=3, ET=4;
      const int IrR=5; 
      const int IuR=6, IuT=7;
      const int CaCb=8;
      const int CbC=9; 
      const int TS=10, TE=11, TIu=12; 
      int u,v,w;
      int day = (int) t+1;
      // note that day 1 corresponds to continuous time t in (0,1]

      const double* alpha = day<15 ? alpha_be : alpha_af;
      const double* Beta = day<15 ? Beta_be : Beta_af;
      const double* mu = day<15 ? mu_be : mu_af;
      const double* Z = day<15 ? Z_be : Z_af;
      const double* D = day<15 ? D_be : D_af;
      const double* Td = day<15 ? Td_be : Td_af;
      
      for (u = 0 ; u < U ; u++) {
        rate[ST]=rate[TS]=rate[ET]=rate[TE]=rate[IuT]=rate[TIu]=0;
	// immigration rates are total, for a Poisson entry model
        // emigration rates are per capita, for a multinomial exit model	
        if (day<15) {
          for (v=0; v < U ; v++) {
            rate[TS] +=mob[u*14-1+day][v]*S[v]/(pop[v]-Ir[v]);
            rate[ST] +=mob[v*14-1+day][u]/(pop[u]-Ir[u]);
            rate[TE] +=mob[u*14-1+day][v]*E[v]/(pop[v]-Ir[v]);
            rate[ET] +=mob[v*14-1+day][u]/(pop[u]-Ir[u]);
            rate[TIu]+=mob[u*14-1+day][v]*Iu[v]/(pop[v]-Ir[v]);
            rate[IuT]+=mob[v*14-1+day][u]/(pop[u]-Iu[u]);
          }
        } else {
          if (u==0) {
	    // u=0 is Wuhan: Wuhan travel multiplied by 0.02, elsewhere by 0.2
            for (v=0; v < U ; v++) {
              rate[TS] +=0.02*mob[u*14-1+14][v]*S[v]/(pop[v]-Ir[v]);
              rate[ST] +=0.02*mob[v*14-1+14][u]/(pop[u]-Ir[u]);
              rate[TE] +=0.02*mob[u*14-1+14][v]*E[v]/(pop[v]-Ir[v]);
              rate[ET] +=0.02*mob[v*14-1+14][u]/(pop[u]-Ir[u]);
              rate[TIu]+=0.02*mob[u*14-1+14][v]*Iu[v]/(pop[v]-Ir[v]);
              rate[IuT]+=0.02*mob[v*14-1+14][u]/(pop[u]-Iu[u]);
            }
          } else {
            rate[TS] +=0.02*mob[u*14-1+14][0]*S[0]/(pop[0]-Ir[0]);
            rate[ST] +=0.02*mob[0*14-1+14][u]/(pop[u]-Ir[u]);
            rate[TE] +=0.02*mob[u*14-1+14][0]*E[0]/(pop[0]-Ir[0]);
            rate[ET] +=0.02*mob[0*14-1+14][u]/(pop[u]-Ir[u]);
            rate[TIu]+=0.02*mob[u*14-1+14][0]*Iu[0]/(pop[0]-Ir[0]);
            rate[IuT]+=0.02*mob[0*14-1+14][u]/(pop[u]-Iu[u]);
          
            for (v=1; v < U ; v++) {
              rate[TS] +=0.2*mob[u*14-1+14][v]*S[v]/(pop[v]-Ir[v]);
              rate[ST] +=0.2*mob[v*14-1+14][u]/(pop[u]-Ir[u]);
              rate[TE] +=0.2*mob[u*14-1+14][v]*E[v]/(pop[v]-Ir[v]);
              rate[ET] +=0.2*mob[v*14-1+14][u]/(pop[u]-Ir[u]);
              rate[TIu]+=0.2*mob[u*14-1+14][v]*Iu[v]/(pop[v]-Ir[v]);
              rate[IuT]+=0.2*mob[v*14-1+14][u]/(pop[u]-Iu[u]);
            }
          }
        }

        rate[TS] = (theta[theta_expand*u])*rate[TS];
        rate[TE] = (theta[theta_expand*u])*rate[TE];
        rate[TIu] = (theta[theta_expand*u])*rate[TIu];
        rate[ST] = (theta[theta_expand*u])*rate[ST];
        rate[ET] = (theta[theta_expand*u])*rate[ET];
        rate[IuT] = (theta[theta_expand*u])*rate[IuT];
        rate[SE] = (Beta[Beta_expand*u]) * (Ir[u]+(mu[mu_expand*u])*Iu[u])/pop[u];
        rate[EIr] = (alpha[alpha_expand*u])/(Z[Z_expand*u]);
        rate[EIu] = (1-(alpha[alpha_expand*u]))/(Z[Z_expand*u]);
        rate[IrR] = 1/(D[D_expand*u]);
        rate[IuR] = 1/(D[D_expand*u]);
        rate[CaCb] = 2/(Td[Td_expand*u]);
        rate[CbC] = 2/(Td[Td_expand*u]);
        for (w=0; w < numTrans ; w++) {
          rate[w] = rate[w]>0 ? rate[w] : 0;
        }
        rate[SE] = rate[SE]*rgammawn(sigma_SE[sigma_SE_expand*u],dt)/dt;
        
        reulermultinom(2,S[u], &rate[SE],dt,&dN[SE]);
        reulermultinom(3,E[u], &rate[EIr],dt,&dN[EIr]);
        reulermultinom(1,Ir[u],&rate[IrR],dt,&dN[IrR]);
        reulermultinom(2,Iu[u],&rate[IuR],dt,&dN[IuR]);
        reulermultinom(1,Ca[u],&rate[CaCb],dt,&dN[CaCb]);
        reulermultinom(1,Cb[u],&rate[CbC],dt,&dN[CbC]);
	dN[TS] = rpois(rate[TS]*dt);
	dN[TE] = rpois(rate[TE]*dt);
	dN[TIu] = rpois(rate[TIu]*dt);

        // transitions between classes
        S[u]  += dN[TS]   - dN[ST]  - dN[SE];
        E[u]  += dN[SE]   + dN[TE]  - dN[ET] - dN[EIr] - dN[EIu];
        Ir[u] += dN[EIr]  - dN[IrR];
        Iu[u] += dN[EIu]  - dN[IuR] + dN[TIu] - dN[IuT];
        Ca[u] += dN[EIr]  - dN[CaCb];
        Cb[u] += dN[CaCb] - dN[CbC];
        C[u]  += dN[CbC];
      
     }
  ')

  covid_dmeasure <- spatPomp_Csnippet(
    unit_obsnames='cases',
    unit_paramnames='tau',
    unit_statenames='C',
    code="
      double m,v,A,B;
      int u;

      lik = 0;
      for (u = 0; u < U; u++) {
        m = C[u];
        v = m*tau[u*tau_expand]*m*tau[u*tau_expand]+m > 4 ? m*tau[u*tau_expand]*m*tau[u*tau_expand]+m : 4;

        if (cases[u] > m) {
	  A = pnorm(cases[u]-0.5,m,sqrt(v),0,1);
	  B = pnorm(cases[u]+0.5,m,sqrt(v),0,1);
	  lik += A + log(1-exp(B-A));
        } else if (cases[u] > 0) {
          lik += log(pnorm(cases[u]+0.5,m,sqrt(v),1,0)-
            pnorm(cases[u]-0.5,m,sqrt(v),1,0));
        } else {
          lik += log(pnorm(0.5,m,sqrt(v),1,0));
        }
        
      }
  
      if(!give_log) lik = exp(lik);
    ")
  
  covid_rmeasure <- spatPomp_Csnippet(
    unit_paramnames='tau',
    unit_statenames='C',
    code="
      double* cases=&cases1;
      double m,v;
      int u;

      for (u = 0; u < U; u++) {
        m = C[u];
        v = m*tau[u*tau_expand]*m*tau[u*tau_expand]+m > 4 ? m*tau[u*tau_expand]*m*tau[u*tau_expand]+m : 4;
        cases[u] = rnorm(m,sqrt(v));
        if (cases[u] > 0.0) {
          cases[u] = nearbyint(cases[u]);
        } else {
          cases[u] = 0.0;
        }
      } 
  ")
  
  covid_dunit_measure <- spatPomp_Csnippet(
    unit_paramnames='tau',
    code="
      double A,B;
      double m = C;
      double v = m*tau[u*tau_expand]*m*tau[u*tau_expand]+m > 4 ? m*tau[u*tau_expand]*m*tau[u*tau_expand]+m : 4;

      if (cases > m) {
	  A = pnorm(cases-0.5,m,sqrt(v),0,1);
	  B = pnorm(cases+0.5,m,sqrt(v),0,1);
	  lik = A + log(1-exp(B-A));
      } else if (cases > 0) {
	  A = pnorm(cases+0.5,m,sqrt(v),1,1);
	  B = pnorm(cases-0.5,m,sqrt(v),1,1);
	  lik = A + log(1-exp(B-A));
      } else {
        lik = pnorm(cases+0.5,m,sqrt(v),1,1);
      }
      if(!give_log) lik = exp(lik);
    "
  )
  
  covid_eunit_measure <- Csnippet("
    ey = C;
  ")
  
  covid_vunit_measure <- spatPomp_Csnippet(
    unit_paramnames='tau',
    code="
      double m;
      m = C;
      vc = m*tau[u*tau_expand]*m*tau[u*tau_expand]+m > 4 ? m*tau[u*tau_expand]*m*tau[u*tau_expand]+m : 4;
    ")

  covid_rinit <- spatPomp_Csnippet(
    unit_statenames=c('S','E','Ir','Iu','Ca','Cb','C'),
    unit_covarnames=c('pop'),
    unit_paramnames=c('E_0','Iu_0'),
    code="
      int u;
      for (u = 0; u < U; u++) {
        E[u] = Iu[u] = Ir[u] = Ca[u] = Cb[u] = C[u] = 0;
        S[u] = nearbyint(pop[u]);
      }
      E[0]=nearbyint(E_0[0]);
      Iu[0]=nearbyint(Iu_0[0]);
      for (u = 1; u < U; u++) {
        E[u]=nearbyint(3*mob[14*u][0]*E_0[u*E_0_expand]/(pop[0]));
        Iu[u]=nearbyint(3*mob[14*u][0]*Iu_0[u*Iu_0_expand]/(pop[0]));
      }
    ")
  
  if(for_ibpf == T) {
    covid_partrans <- parameter_trans(
      logit=c(
        paste0("alpha_be",seq(1,U)),
        paste0("alpha_af",seq(1,U)),
        paste0("mu_be",seq(1,U)),
        paste0("mu_af",seq(1,U))
      ),
      log=c(
        paste0("Beta_be",seq(1,U)),
        paste0("Beta_af",seq(1,U)),
        paste0("Z_be",seq(1,U)),
        paste0("Z_af",seq(1,U)),
        paste0("D_be",seq(1,U)),
        paste0("D_af",seq(1,U)),
        paste0("Td_be",seq(1,U)),
        paste0("Td_af",seq(1,U)),
        paste0("theta",seq(1,U)),
        paste0("tau",seq(1,U)),
        paste0("sigma_SE",seq(1,U)),
        paste0("E_0",seq(1,U)),
        paste0("Iu_0",seq(1,U))
      )
    )
  } else {
    covid_partrans <- parameter_trans(
      logit=c(
        paste0("alpha_be",seq(1,ifelse("alpha"%in%sharedTwoInterval,1,U))),
        paste0("alpha_af",seq(1,ifelse("alpha"%in%sharedTwoInterval,1,U))),
        paste0("mu_be",seq(1,ifelse("mu"%in%sharedTwoInterval,1,U))),
        paste0("mu_af",seq(1,ifelse("mu"%in%sharedTwoInterval,1,U)))
      ),
      log=c(
        paste0("Beta_be",seq(1,ifelse("Beta"%in%sharedTwoInterval,1,U))),
        paste0("Beta_af",seq(1,ifelse("Beta"%in%sharedTwoInterval,1,U))),
        paste0("Z_be",seq(1,ifelse("Z"%in%sharedTwoInterval,1,U))),
        paste0("Z_af",seq(1,ifelse("Z"%in%sharedTwoInterval,1,U))),
        paste0("D_be",seq(1,ifelse("D"%in%sharedTwoInterval,1,U))),
        paste0("D_af",seq(1,ifelse("D"%in%sharedTwoInterval,1,U))),
        paste0("Td_be",seq(1,ifelse("Td"%in%sharedTwoInterval,1,U))),
        paste0("Td_af",seq(1,ifelse("Td"%in%sharedTwoInterval,1,U))),
        paste0("theta",seq(1,ifelse("theta"%in%sharedOneInterval,1,U))),
        paste0("tau",seq(1,ifelse("tau"%in%sharedOneInterval,1,U))) ,
        paste0("sigma_SE",seq(1,ifelse("sigma_SE"%in%sharedOneInterval,1,U))),
        paste0("E_0",seq(1,ifelse("E_0"%in%sharedOneInterval,1,U))),
        paste0("Iu_0",seq(1,ifelse("Iu_0"%in%sharedOneInterval,1,U)))
      )
    )
  }
  

  po <- spatPomp(
    covid_cases,
    units = "city",
    times = "date",
    t0 = min(covid_cases$date),
    unit_statenames = covid_unit_statenames,
    covar = covid_covar,
    rprocess=euler(covid_rprocess, delta.t=dt),
    unit_accumvars = c("C"),
    paramnames=covid_paramnames,
    partrans=covid_partrans,
    globals=covid_globals,
    rinit=covid_rinit,
    dmeasure=covid_dmeasure,
    eunit_measure=covid_eunit_measure,
    vunit_measure=covid_vunit_measure,
    rmeasure=covid_rmeasure,
    dunit_measure=covid_dunit_measure
  )
  coef(po) <- covid_testPar
  
  # set up covariate table to use piecewise constant interpolation
  # noting that we want to look up the covariate for day n when t is in [n-1,n)
  # not currently convenient via spatPomp(), so we manipulate the resulting pomp
  po@t0 <- 0
  covar_tab <- po@covar
  covar_tab@order <- 0L
  covar_tab@times <- c(-1,covar_tab@times-1,max(covar_tab@times)+1) 
  covar_tab@table <- cbind(covar_tab@table[,1],covar_tab@table,covar_tab@table[,ncol(covar_tab@table)])
  po@covar <- covar_tab
  
  po
}
