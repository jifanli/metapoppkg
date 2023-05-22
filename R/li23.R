#' Build spatPomp object for li23
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
li23 <- function(U = 373, dt = 1/4,
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
      sigma_SE=0.78,
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
      sigma_SE=0.78,
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
      sigma_SE=0.78,
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
      int d,u,v,w;
      double factor_uv;
      int day = (int) t+1;  // day 1 corresponds to continuous time t in (0,1]

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
	factor_uv = 1;
        for (v=0; v < U ; v++) {
          if (day>=15) factor_uv = (u==0)||(v==0) ? 0.02 : 0.2 ;
	  d = (day>=15) ? 14 : day; // travel for days 15-30 is set equal to day 14
          rate[TS] +=factor_uv * mob[u*14-1+d][v]*S[v]/(pop[v]-Ir[v]);
          rate[ST] +=factor_uv * mob[v*14-1+d][u]/(pop[u]-Ir[u]);
          rate[TE] +=factor_uv * mob[u*14-1+d][v]*E[v]/(pop[v]-Ir[v]);
          rate[ET] +=factor_uv * mob[v*14-1+d][u]/(pop[u]-Ir[u]);
          rate[TIu]+=factor_uv * mob[u*14-1+d][v]*Iu[v]/(pop[v]-Ir[v]);
          rate[IuT]+=factor_uv * mob[v*14-1+d][u]/(pop[u]-Ir[u]);
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

  covid_skel <- spatPomp_Csnippet(
    unit_vfnames=c('S','E','Ir','Iu','Ca','Cb','C'),
    unit_paramnames=c("alpha_be","Beta_be","alpha_af","Beta_af","mu_be","Z_be",
      "D_be","mu_af","Z_af","D_af","theta","Td_be","Td_af","sigma_SE"),
    unit_covarnames=c('pop'),
    code='
      const double* S=&S1;
      const double* E=&E1;
      const double* Iu=&Iu1;
      const double* Ir=&Ir1;
      const double* Ca=&Ca1;
      const double* Cb=&Cb1;
      const double* C=&C1;
      const int numTrans = 13; // number of transitions
      double rate[numTrans];
      const int SE=0, ST=1;
      const int EIr=2, EIu=3, ET=4;
      const int IrR=5; 
      const int IuR=6, IuT=7;
      const int CaCb=8;
      const int CbC=9; 
      const int TS=10, TE=11, TIu=12; 
      int d,u,v,w;
      double factor_uv;
      int day = (int) t+1;  // day 1 corresponds to continuous time t in (0,1]

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
	factor_uv = 1;
        for (v=0; v < U ; v++) {
          if (day>=15) factor_uv = (u==0)||(v==0) ? 0.02 : 0.2 ;
	  d = (day>=15) ? 14 : day; // travel for days 15-30 is set equal to day 14
          rate[TS] +=factor_uv * mob[u*14-1+d][v]*S[v]/(pop[v]-Ir[v]);
          rate[ST] +=factor_uv * mob[v*14-1+d][u]/(pop[u]-Ir[u]);
          rate[TE] +=factor_uv * mob[u*14-1+d][v]*E[v]/(pop[v]-Ir[v]);
          rate[ET] +=factor_uv * mob[v*14-1+d][u]/(pop[u]-Ir[u]);
          rate[TIu]+=factor_uv * mob[u*14-1+d][v]*Iu[v]/(pop[v]-Ir[v]);
          rate[IuT]+=factor_uv * mob[v*14-1+d][u]/(pop[u]-Ir[u]);
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
       
        // rates are per capita, except for transitions from a reservoir
	// this is for consistency with the stochastic version of the model
        DS[u]  = rate[TS]   - (rate[ST] + rate[SE])* S[u];
        DE[u]  = rate[TE] + rate[SE]*S[u]   - (rate[ET] + rate[EIr] +rate[EIu])*E[u];
        DIr[u] = rate[EIr]*E[u]  - rate[IrR]*Ir[u];
        DIu[u] = rate[TIu] + rate[EIu]*E[u]  - (rate[IuR] + rate[IuT])*Iu[u];
        DCa[u] = rate[EIr]*E[u]  - rate[CaCb]*Ca[u];
        DCb[u] = rate[CaCb]*Ca[u] - rate[CbC]*Cb[u];
        DC[u]  = rate[CbC]*Cb[u];
      
     }
  ')

# We would like to remove tol at some point, and/or move
# to a negative binomial model. vtol has an interpretation 
# as a false positive rate, though number of false positives
# will depend on number of tests and hence cases.
  covid_dmeasure <- Csnippet("
    const double *C = &C1;
    const double *cases = &cases1;  // reported infected
    double m,v;
    double tol = 1e-300;
    double vtol = 1e-5;
    const double *tau = &tau1;
    int u;
    lik = 0;
    for (u = 0; u < U; u++) {
      m = C[u];
      v = m*tau[u*tau_expand]*m*tau[u*tau_expand]+m+vtol;
      if (C < 0) {lik += log(tol);} else {
        if (cases[u] > tol) {
          lik += log(pnorm(cases[u]+0.5,m,sqrt(v),1,0)-
            pnorm(cases[u]-0.5,m,sqrt(v)+tol,1,0) + tol);
        } else {
          lik += log(pnorm(cases[u]+0.5,m,sqrt(v)+tol,1,0)+tol);
        }
      }
    }
    if(!give_log) lik = (lik > log(tol)) ? exp(lik) : tol;
  ")

  covid_rmeasure <- spatPomp_Csnippet(
    unit_paramnames='tau',
    unit_statenames='C',
    code="
      double* cases=&cases1;
      double m,v;
      double vtol = 1e-5;
      int u;

      for (u = 0; u < U; u++) {
        m = C[u];
        v = m*tau[u*tau_expand]*m*tau[u*tau_expand]+m+vtol;
        cases[u] = rnorm(m,sqrt(v));
        if (cases[u] > 0.0) {
          cases[u] = nearbyint(cases[u]);
        } else {
          cases[u] = 0.0;
        }
      } 
    "
  )

  covid_dunit_measure <- Csnippet("
    double vtol = 1e-5;
    double m = C;
    const double *tau = &tau1;
    double v = m*tau[u*tau_expand]*m*tau[u*tau_expand]+m+vtol;
    double tol = 1e-300;
    if (C < 0) {lik = 0;} else {
      if (cases > tol) {
        lik = pnorm(cases+0.5,m,sqrt(v),1,0)-
          pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
      } else {
        lik = pnorm(cases+0.5,m,sqrt(v),1,0)+tol;
      }
    }
    if(give_log) lik = log(lik);
  ")
   
  covid_eunit_measure <- Csnippet("
    ey = C;
  ")

  covid_vunit_measure <- spatPomp_Csnippet(
    unit_paramnames='tau',
    code="
      double m;
      double vtol = 1e-5;
      m = C;
      vc = m*tau[u*tau_expand]*m*tau[u*tau_expand]+m+vtol;
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
    "
  )

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
    skeleton=vectorfield(covid_skel),
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
