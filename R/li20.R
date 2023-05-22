#' Covid spatPomp generator v1
#'
#' Generate a \sQuote{spatPomp} object for reproducing the results of Li et al.
#'
#' @param U A length-one numeric signifying the number of cities to be represented in the spatPomp object.
#' @param dt step size, in days, for the euler approximation.
#' @param sharedParNames estimated parameters that are equal for each unit.
#' @param unitParNames estimated parameters that are different for each unit.
#' @param for_ibpf whether this object is generated for ibpf.
#' @param overdispersion_dynamic whether to include overdispersion in dynamic model.
#' @param mob_modify_factor A multiplicative factor, which is used to control the magnitude of adjustment for the mobility data
#' @param version settings to reproduce results from Li et al (2020) or other selected parameter choices 
#' @param days number of time points, in days
#' @param measurement gives alternative sets of measurement model components
#' @import pomp
#' @import spatPomp
#' @return An object of class \sQuote{spatPomp} representing a \code{U}-dimensional spatially coupled Covid POMP model.
#'
#' @export

li20 <- function(U = 373, dt = 1, 
                 # tau is not used in the measurement model, but it is provided for the convenience of using the function par_li.
                 sharedParNames = c("alpha_be","Beta_be","alpha_af","Beta_af",
                                    "mu_be","Z_be","D_be","mu_af","Z_af","D_af",
                                    "theta","tau","Td_be","Td_af","sigma_SE",
                                    "E_0","Iu_0"),
                 unitParNames = NULL, for_ibpf = T, 
                 overdispersion_dynamic = F,
                 version = c("li20period3"), mob_modify_factor = 0, days = 30,
                 measurement = "simulation"
                 ) {
  if(version[1] == "MLEperiod3") {
    # mle from e3, loglik = -9200
    testPar <- c(
      alpha_be1=0.1199543,
      Beta_be1=0.8053392,
      alpha_af1=0.4151188,
      Beta_af1=0.3441618, 
      mu_be1=0.9700681,
      Z_be1=1.735166,
      D_be1=4.635517,
      mu_af1=0.7220044,
      Z_af1=1.536596,
      D_af1=2.714057,
      theta1=1.752231,
      tau1=0.3456199,
      sigma_SE1=0.994764,
      Td_be1=9,
      Td_af1=6,
      E_01=934.2861,
      Iu_01=1435.945
    )
  } else if(version[1] == "li20period3") {
    testPar <- c(
      alpha_be1=0.14,
      Beta_be1=1.12,
      mu_be1=0.55,
      Z_be1=3.69,
      D_be1=3.47,
      Td_be1=9,
      alpha_af1=0.69,
      Beta_af1=0.35,
      mu_af1=0.43,
      Z_af1=3.42,
      D_af1=3.31,
      Td_af1=6,
      theta1=1.36,
      tau1=0.5,
      sigma_SE1=0.5, 
      E_01=1000,
      Iu_01=1000
    )
  } else if(version[1] == "li20period2") {
    testPar <- c(
      alpha_be1=0.14,
      Beta_be1=1.12,
      mu_be1=0.55,
      Z_be1=3.69,
      D_be1=3.47,
      alpha_af1=0.65,
      Beta_af1=0.52,
      mu_af1=0.5,
      Z_af1=3.6,
      D_af1=3.14,
      theta1=1.36,
      tau1=0.5,
      sigma_SE1=0.5, 
      Td_be1=9,
      Td_af1=6,
      E_01=1000,
      Iu_01=1000
    )
  } else if(version[1] == "li20period1") {
    testPar = c(
      alpha_be1=0.14,
      Beta_be1=1.12,
      mu_be1=0.55,
      Z_be1=3.69,
      D_be1=3.47,
      # period1 won't use _af, provide them to avoid an error
      alpha_af1=0.69,
      Beta_af1=0.35,
      mu_af1=0.43,
      Z_af1=3.42,
      D_af1=3.31,
      theta1=1.36,
      tau1=0.5,
      sigma_SE1=0.5, 
      Td_be1=9,
      Td_af1=6,
      E_01=1000,
      Iu_01=1000
    )
  }
  testPar <- par_li(testPar,U=U,dir="to_li")
  
  mobi <- metapoppkg::mobility+metapoppkg::v_by_g_day*mob_modify_factor
  popu <- metapoppkg::population[1:days+rep((0:372)*30,each=days),]
  incidence <- metapoppkg::incidence # cases
  
  if(for_ibpf==T){
    unitParNames <- c("alpha_be","Beta_be","alpha_af","Beta_af","mu_be","Z_be",
                     "D_be","mu_af","Z_af","D_af","theta","tau","Td_be",
                     "Td_af","sigma_SE","E_0","Iu_0")
    sharedParNames <- NULL
  }
  
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
  covid_cases$city <- as.character(covid_cases$city)
  rownames(covid_covar) <- NULL
  
  set_expanded <- Csnippet(paste0("const int ", unitParNames,
                                 "_expand = 1;\n", collapse=" "))
  set_fixed <- Csnippet(paste0("const int ", sharedParNames,
                              "_expand = 0;\n", collapse=" "))
  dynamic_C <- Csnippet(paste0("const int overdispersion_dynamic = ",
                              as.numeric(overdispersion_dynamic),";\n"))
  
  if(!is.null(unitParNames) & !is.null(sharedParNames)){
    covid_globals <- Csnippet(
      paste(mobi_C, set_expanded, set_fixed, dynamic_C, sep = "\n")
    )
  }
  if(is.null(unitParNames) & !is.null(sharedParNames)){
    covid_globals <- Csnippet(
      paste(mobi_C, set_fixed, dynamic_C, sep = "\n")
    )
  }
  if(!is.null(unitParNames) & is.null(sharedParNames)){
    covid_globals <- Csnippet(
      paste(mobi_C, set_expanded, dynamic_C, sep = "\n")
    )
  }
  
  covid_unit_statenames <- c('S','E','Iu','Ir','Ca','Cb','C')
  if(!is.null(unitParNames) & !is.null(sharedParNames)){
    covid_paramnames <- c(paste0(rep(unitParNames,each=U),1:U),paste0(sharedParNames,"1"))
  }
  if(is.null(unitParNames) & !is.null(sharedParNames)){
    covid_paramnames <- c(paste0(sharedParNames,"1"))
  }
  if(!is.null(unitParNames) & is.null(sharedParNames)){
    covid_paramnames <- c(paste0(rep(unitParNames,each=U),1:U))
  }
  
  covid_testPar <- c(
    if(!is.null(unitParNames)) rep(testPar[paste0(unitParNames,"1")],each=U) else NULL,
    if(!is.null(sharedParNames)) testPar[paste0(sharedParNames,"1")] else NULL
  )
  names(covid_testPar) <- covid_paramnames

  covid_partrans <- parameter_trans(
    logit=covid_paramnames
  )
  
  covid_rprocess <- Csnippet('
    double trans[15];
    double *S = &S1;
    double *E = &E1;
    double *Ir = &Ir1;
    double *Iu = &Iu1;
    double *Ca = &Ca1; // Reporting delay states
    double *Cb = &Cb1;
    double *C = &C1;
    const double *pop = &pop1;
    const double *alpha_be = &alpha_be1;
    const double *alpha_af = &alpha_af1;
    const double *Beta_be = &Beta_be1;
    const double *Beta_af = &Beta_af1;
    const double *mu_af = &mu_af1;
    const double *Z_af = &Z_af1;
    const double *D_af = &D_af1;
    const double *mu_be = &mu_be1;
    const double *Z_be = &Z_be1;
    const double *D_be = &D_be1;
    const double *theta = &theta1;
    const double *Td_be = &Td_be1;
    const double *Td_af = &Td_af1;
    const double *sigma_SE = &sigma_SE1;

    double uc; // c means third, d means fourth
    double ud;
    double ug;
    double uh;
    double uk;
    double ul;
    double ska[U];
    double eka[U];
    double iska[U];
    double iaka[U];
    double caka[U];
    double cbka[U];
    double cka[U];
    double Tsa[U];
    double Tea[U];
    double Tisa[U];
    double Tiaa[U];
    double Tcaa[U];
    double Tcba[U];
    double skb[U];
    double ekb[U];
    double iskb[U];
    double iakb[U];
    double cakb[U];
    double cbkb[U];
    double ckb[U];
    double Tsb[U];
    double Teb[U];
    double Tisb[U];
    double Tiab[U];
    double Tcab[U];
    double Tcbb[U];
    double skc[U];
    double ekc[U];
    double iskc[U];
    double iakc[U];
    double cakc[U];
    double cbkc[U];
    double ckc[U];
    double Tsc[U];
    double Tec[U];
    double Tisc[U];
    double Tiac[U];
    double Tcac[U];
    double Tcbc[U];
    double skd[U];
    double ekd[U];
    double iskd[U];
    double iakd[U];
    double cakd[U];
    double cbkd[U];
    double ckd[U];
    int u,v,w;
    int day = (int) t+1;

    for (u = 0 ; u < U ; u++) {
      uc=0;
      ud=0;
      ug=0;
      uh=0;
      uk=0;
      ul=0;
      if(day<15){
        for (v=0; v < U ; v++) {
          uc+=mob[u*14-1+day][v]*S[v]/(pop[v]-Ir[v]);
          ud+=mob[v*14-1+day][u]*S[u]/(pop[u]-Ir[u]);
          ug+=mob[u*14-1+day][v]*E[v]/(pop[v]-Ir[v]);
          uh+=mob[v*14-1+day][u]*E[u]/(pop[u]-Ir[u]);
          uk+=mob[u*14-1+day][v]*Iu[v]/(pop[v]-Ir[v]);
          ul+=mob[v*14-1+day][u]*Iu[u]/(pop[u]-Iu[u]);
        }
      } else {
        if (u==0) {
	    // u=0 is Wuhan: Wuhan travel multiplied by 0.02, elsewhere by 0.2
          for (v=0; v < U ; v++) {
            uc +=0.02*mob[u*14-1+14][v]*S[v]/(pop[v]-Ir[v]);
            ud +=0.02*mob[v*14-1+14][u]*S[u]/(pop[u]-Ir[u]);
            ug +=0.02*mob[u*14-1+14][v]*E[v]/(pop[v]-Ir[v]);
            uh +=0.02*mob[v*14-1+14][u]*E[u]/(pop[u]-Ir[u]);
            uk +=0.02*mob[u*14-1+14][v]*Iu[v]/(pop[v]-Ir[v]);
            ul +=0.02*mob[v*14-1+14][u]*Iu[u]/(pop[u]-Iu[u]);
          }
        } else {
          uc +=0.02*mob[u*14-1+14][0]*S[0]/(pop[0]-Ir[0]);
          ud +=0.02*mob[0*14-1+14][u]*S[u]/(pop[u]-Ir[u]);
          ug +=0.02*mob[u*14-1+14][0]*E[0]/(pop[0]-Ir[0]);
          uh +=0.02*mob[0*14-1+14][u]*E[u]/(pop[u]-Ir[u]);
          uk +=0.02*mob[u*14-1+14][0]*Iu[0]/(pop[0]-Ir[0]);
          ul +=0.02*mob[0*14-1+14][u]*Iu[u]/(pop[u]-Iu[u]);
          
          for (v=1; v < U ; v++) {
            uc +=0.2*mob[u*14-1+14][v]*S[v]/(pop[v]-Ir[v]);
            ud +=0.2*mob[v*14-1+14][u]*S[u]/(pop[u]-Ir[u]);
            ug +=0.2*mob[u*14-1+14][v]*E[v]/(pop[v]-Ir[v]);
            uh +=0.2*mob[v*14-1+14][u]*E[u]/(pop[u]-Ir[u]);
            uk +=0.2*mob[u*14-1+14][v]*Iu[v]/(pop[v]-Ir[v]);
            ul +=0.2*mob[v*14-1+14][u]*Iu[u]/(pop[u]-Iu[u]);
          }
        }
      }
      
      trans[1] = (0.8+0.7*Beta_be[Beta_be_expand*u])*S[u]/pop[u]*Ir[u]*dt;
      trans[2] = (0.8+0.7*Beta_be[Beta_be_expand*u])*S[u]/pop[u]*(0.2+0.8*mu_be[mu_be_expand*u])*Iu[u]*dt;
      trans[3] = (1+0.75*theta[theta_expand*u])*uc*dt;
      trans[4] = ((1+0.75*theta[theta_expand*u])*ud*dt)<(S[u]*dt) ? (1+0.75*theta[theta_expand*u])*ud*dt : S[u]*dt;        
      trans[5] = (0.02+0.98*alpha_be[alpha_be_expand*u])*E[u]/(2+3*Z_be[Z_be_expand*u])*dt;
      trans[6] = (1-(0.02+0.98*alpha_be[alpha_be_expand*u]))*E[u]/(2+3*Z_be[Z_be_expand*u])*dt;
      trans[7] = (1+0.75*theta[theta_expand*u])*ug*dt;
      trans[8] = ((1+0.75*theta[theta_expand*u])*uh*dt)<(E[u]*dt) ? (1+0.75*theta[theta_expand*u])*uh*dt : E[u]*dt;
      trans[9] = Ir[u]/(2+3*D_be[D_be_expand*u])*dt;
      trans[10] = Iu[u]/(2+3*D_be[D_be_expand*u])*dt;
      trans[11] = (1+0.75*theta[theta_expand*u])*uk*dt;
      trans[12] = ((1+0.75*theta[theta_expand*u])*ul*dt)<(Iu[u]*dt) ? (1+0.75*theta[theta_expand*u])*ul*dt : Iu[u]*dt;
      trans[13] = Ca[u]*2/(5+5*Td_be[Td_be_expand*u])*dt;
      trans[14] = Cb[u]*2/(5+5*Td_be[Td_be_expand*u])*dt;
      for (w=1; w < 15 ; w++) {
        trans[w] = trans[w]>0 ? trans[w] : 0;
      }
      if (overdispersion_dynamic>0.5){
        trans[1] = rpois(trans[1]*rgammawn(5*sigma_SE[sigma_SE_expand*u],dt)/dt);
        trans[2] = rpois(trans[2]*rgammawn(5*sigma_SE[sigma_SE_expand*u],dt)/dt);
      // Rprintf("test if loop is run"); It has been tested, when overdispersion_dynamic = F, the loop will not be run
      } else {
        trans[1] = rpois(trans[1]);
        trans[2] = rpois(trans[2]);
      }
      for (w=3; w < 13 ; w++) {
        trans[w] = rpois(trans[w]);
      }

      // transitions between classes
      ska[u] = - trans[1] - trans[2] + trans[3] - trans[4];
      eka[u] = trans[1] + trans[2] - trans[5] - trans[6] + trans[7] - trans[8];
      iska[u] = trans[5] - trans[9];
      iaka[u] = trans[6] - trans[10] + trans[11] - trans[12];
      caka[u] = trans[5] - trans[13];
      cbka[u] = trans[13] - trans[14];
      cka[u] = trans[14];
      Tsa[u]=S[u]+ska[u]/2;
      Tea[u]=E[u]+eka[u]/2;
      Tisa[u]=Ir[u]+iska[u]/2;
      Tiaa[u]=Iu[u]+iaka[u]/2;
      Tcaa[u]=Ca[u]+caka[u]/2;
      Tcba[u]=Cb[u]+cbka[u]/2;
    }

    for (u = 0 ; u < U ; u++) {
      uc=0;
      ud=0;
      ug=0;
      uh=0;
      uk=0;
      ul=0;
      if(day<15){
        for (v=0; v < U ; v++) {
          uc+=mob[u*14-1+day][v]*Tsa[v]/(pop[v]-Tisa[v]);
          ud+=mob[v*14-1+day][u]*Tsa[u]/(pop[u]-Tisa[u]);
          ug+=mob[u*14-1+day][v]*Tea[v]/(pop[v]-Tisa[v]);
          uh+=mob[v*14-1+day][u]*Tea[u]/(pop[u]-Tisa[u]);
          uk+=mob[u*14-1+day][v]*Tiaa[v]/(pop[v]-Tisa[v]);
          ul+=mob[v*14-1+day][u]*Tiaa[u]/(pop[u]-Tiaa[u]);
        }
      } else {
        if (u==0) {
	    // u=0 is Wuhan: Wuhan travel multiplied by 0.02, elsewhere by 0.2
          for (v=0; v < U ; v++) {
            uc +=0.02*mob[u*14-1+14][v]*Tsa[v]/(pop[v]-Tisa[v]);
            ud +=0.02*mob[v*14-1+14][u]*Tsa[u]/(pop[u]-Tisa[u]);
            ug +=0.02*mob[u*14-1+14][v]*Tea[v]/(pop[v]-Tisa[v]);
            uh +=0.02*mob[v*14-1+14][u]*Tea[u]/(pop[u]-Tisa[u]);
            uk +=0.02*mob[u*14-1+14][v]*Tiaa[v]/(pop[v]-Tisa[v]);
            ul +=0.02*mob[v*14-1+14][u]*Tiaa[u]/(pop[u]-Tiaa[u]);
          }
        } else {
          uc +=0.02*mob[u*14-1+14][0]*Tsa[0]/(pop[0]-Tisa[0]);
          ud +=0.02*mob[0*14-1+14][u]*Tsa[u]/(pop[u]-Tisa[u]);
          ug +=0.02*mob[u*14-1+14][0]*Tea[0]/(pop[0]-Tisa[0]);
          uh +=0.02*mob[0*14-1+14][u]*Tea[u]/(pop[u]-Tisa[u]);
          uk +=0.02*mob[u*14-1+14][0]*Tiaa[0]/(pop[0]-Tisa[0]);
          ul +=0.02*mob[0*14-1+14][u]*Tiaa[u]/(pop[u]-Tiaa[u]);
          
          for (v=1; v < U ; v++) {
            uc +=0.2*mob[u*14-1+14][v]*Tsa[v]/(pop[v]-Tisa[v]);
            ud +=0.2*mob[v*14-1+14][u]*Tsa[u]/(pop[u]-Tisa[u]);
            ug +=0.2*mob[u*14-1+14][v]*Tea[v]/(pop[v]-Tisa[v]);
            uh +=0.2*mob[v*14-1+14][u]*Tea[u]/(pop[u]-Tisa[u]);
            uk +=0.2*mob[u*14-1+14][v]*Tiaa[v]/(pop[v]-Tisa[v]);
            ul +=0.2*mob[v*14-1+14][u]*Tiaa[u]/(pop[u]-Tiaa[u]);
          }
        }
      }
      trans[1] = (0.8+0.7*Beta_be[Beta_be_expand*u])*Tsa[u]/pop[u]*Tisa[u]*dt;
      trans[2] = (0.8+0.7*Beta_be[Beta_be_expand*u])*Tsa[u]/pop[u]*(0.2+0.8*mu_be[mu_be_expand*u])*Tiaa[u]*dt;
      trans[3] = (1+0.75*theta[theta_expand*u])*uc*dt;
      trans[4] = ((1+0.75*theta[theta_expand*u])*ud*dt)<(Tsa[u]*dt) ? (1+0.75*theta[theta_expand*u])*ud*dt : Tsa[u]*dt;
      trans[5] = (0.02+0.98*alpha_be[alpha_be_expand*u])*Tea[u]/(2+3*Z_be[Z_be_expand*u])*dt;
      trans[6] = (1-(0.02+0.98*alpha_be[alpha_be_expand*u]))*Tea[u]/(2+3*Z_be[Z_be_expand*u])*dt;
      trans[7] = (1+0.75*theta[theta_expand*u])*ug*dt;
      trans[8] = ((1+0.75*theta[theta_expand*u])*uh*dt)<(Tea[u]*dt) ? (1+0.75*theta[theta_expand*u])*uh*dt : Tea[u]*dt;
      trans[9] = Tisa[u]/(2+3*D_be[D_be_expand*u])*dt;
      trans[10] = Tiaa[u]/(2+3*D_be[D_be_expand*u])*dt;
      trans[11] = (1+0.75*theta[theta_expand*u])*uk*dt;
      trans[12] = ((1+0.75*theta[theta_expand*u])*ul*dt)<(Tiaa[u]*dt) ? (1+0.75*theta[theta_expand*u])*ul*dt : Tiaa[u]*dt;
      trans[13] = Tcaa[u]*2/(5+5*Td_be[Td_be_expand*u])*dt;
      trans[14] = Tcba[u]*2/(5+5*Td_be[Td_be_expand*u])*dt;
      for (w=1; w < 15 ; w++) {
        trans[w] = trans[w]>0 ? trans[w] : 0;
      }
      if (overdispersion_dynamic>0.5){
        trans[1] = rpois(trans[1]*rgammawn(5*sigma_SE[sigma_SE_expand*u],dt)/dt);
        trans[2] = rpois(trans[2]*rgammawn(5*sigma_SE[sigma_SE_expand*u],dt)/dt);
      } else {
        trans[1] = rpois(trans[1]);
        trans[2] = rpois(trans[2]);
      }
      for (w=3; w < 13 ; w++) {
        trans[w] = rpois(trans[w]);
      }
  
      // transitions between classes
      skb[u] = - trans[1] - trans[2] + trans[3] - trans[4];
      ekb[u] = trans[1] + trans[2] - trans[5] - trans[6] + trans[7] - trans[8];
      iskb[u] = trans[5] - trans[9];
      iakb[u] = trans[6] - trans[10] + trans[11] - trans[12];
      cakb[u] = trans[5] - trans[13];
      cbkb[u] = trans[13] - trans[14];
      ckb[u] = trans[14];
      Tsb[u]=S[u]+skb[u]/2;
      Teb[u]=E[u]+ekb[u]/2;
      Tisb[u]=Ir[u]+iskb[u]/2;
      Tiab[u]=Iu[u]+iakb[u]/2;
      Tcab[u]=Ca[u]+cakb[u]/2;
      Tcbb[u]=Cb[u]+cbkb[u]/2;
    }

    for (u = 0 ; u < U ; u++) {
      uc=0;
      ud=0;
      ug=0;
      uh=0;
      uk=0;
      ul=0;
      if(day<15){
        for (v=0; v < U ; v++) {
          uc+=mob[u*14-1+day][v]*Tsb[v]/(pop[v]-Tisb[v]);
          ud+=mob[v*14-1+day][u]*Tsb[u]/(pop[u]-Tisb[u]);
          ug+=mob[u*14-1+day][v]*Teb[v]/(pop[v]-Tisb[v]);
          uh+=mob[v*14-1+day][u]*Teb[u]/(pop[u]-Tisb[u]);
          uk+=mob[u*14-1+day][v]*Tiab[v]/(pop[v]-Tisb[v]);
          ul+=mob[v*14-1+day][u]*Tiab[u]/(pop[u]-Tiab[u]);
        }
      } else {
        if (u==0) {
	    // u=0 is Wuhan: Wuhan travel multiplied by 0.02, elsewhere by 0.2
          for (v=0; v < U ; v++) {
            uc +=0.02*mob[u*14-1+14][v]*Tsb[v]/(pop[v]-Tisb[v]);
            ud +=0.02*mob[v*14-1+14][u]*Tsb[u]/(pop[u]-Tisb[u]);
            ug +=0.02*mob[u*14-1+14][v]*Teb[v]/(pop[v]-Tisb[v]);
            uh +=0.02*mob[v*14-1+14][u]*Teb[u]/(pop[u]-Tisb[u]);
            uk +=0.02*mob[u*14-1+14][v]*Tiab[v]/(pop[v]-Tisb[v]);
            ul +=0.02*mob[v*14-1+14][u]*Tiab[u]/(pop[u]-Tiab[u]);
          }
        } else {
          uc +=0.02*mob[u*14-1+14][0]*Tsb[0]/(pop[0]-Tisb[0]);
          ud +=0.02*mob[0*14-1+14][u]*Tsb[u]/(pop[u]-Tisb[u]);
          ug +=0.02*mob[u*14-1+14][0]*Teb[0]/(pop[0]-Tisb[0]);
          uh +=0.02*mob[0*14-1+14][u]*Teb[u]/(pop[u]-Tisb[u]);
          uk +=0.02*mob[u*14-1+14][0]*Tiab[0]/(pop[0]-Tisb[0]);
          ul +=0.02*mob[0*14-1+14][u]*Tiab[u]/(pop[u]-Tiab[u]);
          
          for (v=1; v < U ; v++) {
            uc +=0.2*mob[u*14-1+14][v]*Tsb[v]/(pop[v]-Tisb[v]);
            ud +=0.2*mob[v*14-1+14][u]*Tsb[u]/(pop[u]-Tisb[u]);
            ug +=0.2*mob[u*14-1+14][v]*Teb[v]/(pop[v]-Tisb[v]);
            uh +=0.2*mob[v*14-1+14][u]*Teb[u]/(pop[u]-Tisb[u]);
            uk +=0.2*mob[u*14-1+14][v]*Tiab[v]/(pop[v]-Tisb[v]);
            ul +=0.2*mob[v*14-1+14][u]*Tiab[u]/(pop[u]-Tiab[u]);
          }
        }
      }

      trans[1] = (0.8+0.7*Beta_be[Beta_be_expand*u])*Tsb[u]/pop[u]*Tisb[u]*dt;
      trans[2] = (0.8+0.7*Beta_be[Beta_be_expand*u])*Tsb[u]/pop[u]*(0.2+0.8*mu_be[mu_be_expand*u])*Tiab[u]*dt;
      trans[3] = (1+0.75*theta[theta_expand*u])*uc*dt;
      trans[4] = ((1+0.75*theta[theta_expand*u])*ud*dt)<(Tsb[u]*dt) ? (1+0.75*theta[theta_expand*u])*ud*dt : Tsb[u]*dt;
      trans[5] = (0.02+0.98*alpha_be[alpha_be_expand*u])*Teb[u]/(2+3*Z_be[Z_be_expand*u])*dt;
      trans[6] = (1-(0.02+0.98*alpha_be[alpha_be_expand*u]))*Teb[u]/(2+3*Z_be[Z_be_expand*u])*dt;
      trans[7] = (1+0.75*theta[theta_expand*u])*ug*dt;
      trans[8] = ((1+0.75*theta[theta_expand*u])*uh*dt)<(Teb[u]*dt) ? (1+0.75*theta[theta_expand*u])*uh*dt : Teb[u]*dt;
      trans[9] = Tisb[u]/(2+3*D_be[D_be_expand*u])*dt;
      trans[10] = Tiab[u]/(2+3*D_be[D_be_expand*u])*dt;
      trans[11] = (1+0.75*theta[theta_expand*u])*uk*dt;
      trans[12] = ((1+0.75*theta[theta_expand*u])*ul*dt)<(Tiab[u]*dt) ? (1+0.75*theta[theta_expand*u])*ul*dt : Tiab[u]*dt;
      trans[13] = Tcab[u]*2/(5+5*Td_be[Td_be_expand*u])*dt;
      trans[14] = Tcbb[u]*2/(5+5*Td_be[Td_be_expand*u])*dt;
      for (w=1; w < 15 ; w++) {
        trans[w] = trans[w]>0 ? trans[w] : 0;
      }
      if (overdispersion_dynamic>0.5){
        trans[1] = rpois(trans[1]*rgammawn(5*sigma_SE[sigma_SE_expand*u],dt)/dt);
        trans[2] = rpois(trans[2]*rgammawn(5*sigma_SE[sigma_SE_expand*u],dt)/dt);
      } else {
        trans[1] = rpois(trans[1]);
        trans[2] = rpois(trans[2]);
      }
      for (w=3; w < 13 ; w++) {
        trans[w] = rpois(trans[w]);
      }

      // transitions between classes
      skc[u] = - trans[1] - trans[2] + trans[3] - trans[4];
      ekc[u] = trans[1] + trans[2] - trans[5] - trans[6] + trans[7] - trans[8];
      iskc[u] = trans[5] - trans[9];
      iakc[u] = trans[6] - trans[10] + trans[11] - trans[12];
      cakc[u] = trans[5] - trans[13];
      cbkc[u] = trans[13] - trans[14];
      ckc[u] = trans[14];
      Tsc[u]=S[u]+skc[u];
      Tec[u]=E[u]+ekc[u];
      Tisc[u]=Ir[u]+iskc[u];
      Tiac[u]=Iu[u]+iakc[u];
      Tcac[u]=Ca[u]+cakc[u];
      Tcbc[u]=Cb[u]+cbkc[u];
    }

    for (u = 0 ; u < U ; u++) {
      uc=0;
      ud=0;
      ug=0;
      uh=0;
      uk=0;
      ul=0;
      if(day<15){
        for (v=0; v < U ; v++) {
          uc+=mob[u*14-1+day][v]*Tsc[v]/(pop[v]-Tisc[v]);
          ud+=mob[v*14-1+day][u]*Tsc[u]/(pop[u]-Tisc[u]);
          ug+=mob[u*14-1+day][v]*Tec[v]/(pop[v]-Tisc[v]);
          uh+=mob[v*14-1+day][u]*Tec[u]/(pop[u]-Tisc[u]);
          uk+=mob[u*14-1+day][v]*Tiac[v]/(pop[v]-Tisc[v]);
          ul+=mob[v*14-1+day][u]*Tiac[u]/(pop[u]-Tiac[u]);
        }
      } else {
        if (u==0) {
	    // u=0 is Wuhan: Wuhan travel multiplied by 0.02, elsewhere by 0.2
          for (v=0; v < U ; v++) {
            uc +=0.02*mob[u*14-1+14][v]*Tsc[v]/(pop[v]-Tisc[v]);
            ud +=0.02*mob[v*14-1+14][u]*Tsc[u]/(pop[u]-Tisc[u]);
            ug +=0.02*mob[u*14-1+14][v]*Tec[v]/(pop[v]-Tisc[v]);
            uh +=0.02*mob[v*14-1+14][u]*Tec[u]/(pop[u]-Tisc[u]);
            uk +=0.02*mob[u*14-1+14][v]*Tiac[v]/(pop[v]-Tisc[v]);
            ul +=0.02*mob[v*14-1+14][u]*Tiac[u]/(pop[u]-Tiac[u]);
          }
        } else {
          uc +=0.02*mob[u*14-1+14][0]*Tsc[0]/(pop[0]-Tisc[0]);
          ud +=0.02*mob[0*14-1+14][u]*Tsc[u]/(pop[u]-Tisc[u]);
          ug +=0.02*mob[u*14-1+14][0]*Tec[0]/(pop[0]-Tisc[0]);
          uh +=0.02*mob[0*14-1+14][u]*Tec[u]/(pop[u]-Tisc[u]);
          uk +=0.02*mob[u*14-1+14][0]*Tiac[0]/(pop[0]-Tisc[0]);
          ul +=0.02*mob[0*14-1+14][u]*Tiac[u]/(pop[u]-Tiac[u]);
          
          for (v=1; v < U ; v++) {
            uc +=0.2*mob[u*14-1+14][v]*Tsc[v]/(pop[v]-Tisc[v]);
            ud +=0.2*mob[v*14-1+14][u]*Tsc[u]/(pop[u]-Tisc[u]);
            ug +=0.2*mob[u*14-1+14][v]*Tec[v]/(pop[v]-Tisc[v]);
            uh +=0.2*mob[v*14-1+14][u]*Tec[u]/(pop[u]-Tisc[u]);
            uk +=0.2*mob[u*14-1+14][v]*Tiac[v]/(pop[v]-Tisc[v]);
            ul +=0.2*mob[v*14-1+14][u]*Tiac[u]/(pop[u]-Tiac[u]);
          }
        }
      }
      trans[1] = (0.8+0.7*Beta_be[Beta_be_expand*u])*Tsc[u]/pop[u]*Tisc[u]*dt;
      trans[2] = (0.8+0.7*Beta_be[Beta_be_expand*u])*Tsc[u]/pop[u]*(0.2+0.8*mu_be[mu_be_expand*u])*Tiac[u]*dt;
      trans[3] = (1+0.75*theta[theta_expand*u])*uc*dt;
      trans[4] = ((1+0.75*theta[theta_expand*u])*ud*dt)<(Tsc[u]*dt) ? (1+0.75*theta[theta_expand*u])*ud*dt : Tsc[u]*dt;
      trans[5] = (0.02+0.98*alpha_be[alpha_be_expand*u])*Tec[u]/(2+3*Z_be[Z_be_expand*u])*dt;
      trans[6] = (1-(0.02+0.98*alpha_be[alpha_be_expand*u]))*Tec[u]/(2+3*Z_be[Z_be_expand*u])*dt;
      trans[7] = (1+0.75*theta[theta_expand*u])*ug*dt;
      trans[8] = ((1+0.75*theta[theta_expand*u])*uh*dt)<(Tec[u]*dt) ? (1+0.75*theta[theta_expand*u])*uh*dt : Tec[u]*dt;
      trans[9] = Tisc[u]/(2+3*D_be[D_be_expand*u])*dt;
      trans[10] = Tiac[u]/(2+3*D_be[D_be_expand*u])*dt;
      trans[11] = (1+0.75*theta[theta_expand*u])*uk*dt;
      trans[12] = ((1+0.75*theta[theta_expand*u])*ul*dt)<(Tiac[u]*dt) ? (1+0.75*theta[theta_expand*u])*ul*dt : Tiac[u]*dt;
      trans[13] = Tcac[u]*2/(5+5*Td_be[Td_be_expand*u])*dt;
      trans[14] = Tcbc[u]*2/(5+5*Td_be[Td_be_expand*u])*dt;
      for (w=1; w < 15 ; w++) {
        trans[w] = trans[w]>0 ? trans[w] : 0;
      }
      if (overdispersion_dynamic>0.5){
        trans[1] = rpois(trans[1]*rgammawn(5*sigma_SE[sigma_SE_expand*u],dt)/dt);
        trans[2] = rpois(trans[2]*rgammawn(5*sigma_SE[sigma_SE_expand*u],dt)/dt);
      } else {
        trans[1] = rpois(trans[1]);
        trans[2] = rpois(trans[2]);
      }
      for (w=3; w < 13 ; w++) {
        trans[w] = rpois(trans[w]);
      }

      // transitions between classes
      skd[u] = - trans[1] - trans[2] + trans[3] - trans[4];
      ekd[u] = trans[1] + trans[2] - trans[5] - trans[6] + trans[7] - trans[8];
      iskd[u] = trans[5] - trans[9];
      iakd[u] = trans[6] - trans[10] + trans[11] - trans[12];
      cakd[u] = trans[5] - trans[13];
      cbkd[u] = trans[13] - trans[14];
      ckd[u] = trans[14];
    }

    for (u = 0 ; u < U ; u++) {
      S[u]+=ska[u]/6+skb[u]/3+skc[u]/3+skd[u]/6;
      E[u]+=eka[u]/6+ekb[u]/3+ekc[u]/3+ekd[u]/6;
      Ir[u]+=iska[u]/6+iskb[u]/3+iskc[u]/3+iskd[u]/6;
      Iu[u]+=iaka[u]/6+iakb[u]/3+iakc[u]/3+iakd[u]/6;
      Ca[u]+=caka[u]/6+cakb[u]/3+cakc[u]/3+cakd[u]/6;
      Cb[u]+=cbka[u]/6+cbkb[u]/3+cbkc[u]/3+cbkd[u]/6;
      C[u]+=cka[u]/6+ckb[u]/3+ckc[u]/3+ckd[u]/6;
    }
    for (u = 0 ; u < U ; u++) {
      S[u] = S[u]>0 ? S[u] : 0;
      E[u] = E[u]>0 ? E[u] : 0;
      Ir[u] = Ir[u]>0 ? Ir[u] : 0;
      Iu[u] = Iu[u]>0 ? Iu[u] : 0;
      Ca[u] = Ca[u]>0 ? Ca[u] : 0;
      Cb[u] = Cb[u]>0 ? Cb[u] : 0;
      C[u] = C[u]>0 ? C[u] : 0;
    }
  ')

  covid_dmeasure <- Csnippet("
    const double *C = &C1;
    const double *cases = &cases1;  // reported infected
    double m,v,A,B;
    int u;

    lik = 0;
    for (u = 0; u < U; u++) {
      m = C[u];
      v = m*m/4 > 4 ? m*m/4 : 4;

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
  
  covid_dmeasure_simu <- Csnippet("
    const double *C = &C1;
    const double *cases = &cases1;  // reported infected
    double m,v,A,B;
    int u;

    lik = 0;
    for (u = 0; u < U; u++) {
      m = C[u];
      v = m*m/4;

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

  covid_rmeasure <- Csnippet("
    const double *C = &C1;
    double *cases = &cases1;
    double m,v;
    int u;

    for (u = 0; u < U; u++) {
      m = C[u];
      v = m*m/4 > 4 ? m*m/4 : 4;

      cases[u] = rnorm(m,sqrt(v));
      if (cases[u] > 0.0) {
        cases[u] = nearbyint(cases[u]);
      } else {
        cases[u] = 0.0;
      }
    }
  ")
  
  covid_rmeasure_simu <- Csnippet("
    const double *C = &C1;
    double *cases = &cases1;
    double m,v;
    int u;

    for (u = 0; u < U; u++) {
      m = C[u];
      v = m*m/4;

      cases[u] = rnorm(m,sqrt(v));
      if (cases[u] > 0.0) {
        cases[u] = nearbyint(cases[u]);
      } else {
        cases[u] = 0.0;
      }
    }
  ")
  
  covid_dunit_measure <- Csnippet('
      double A,B;
      double m = C;
      double v = m*m/4 > 4 ? m*m/4 : 4;

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
    ')
    
  covid_dunit_measure_simu <- Csnippet('
      double A,B;
      double m = C;
      double v = m*m/4;

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
    ')

  covid_eunit_measure <- Csnippet("
    ey = C;
  ")

  covid_vunit_measure <- Csnippet("
    double m;
    m = C;
    vc = m*m/4 > 4 ? m*m/4 : 4;
  ")
  
  covid_vunit_measure_simu <- Csnippet("
    double m;
    m = C;
    vc = m*m/4;
  ")
  
  covid_rinit <- Csnippet("
    double *S = &S1;
    double *E = &E1;
    double *Iu = &Iu1;
    double *Ir = &Ir1;
    double *Ca = &Ca1;
    double *Cb = &Cb1;
    double *C = &C1;
    const double *pop = &pop1;
    int u;
    const double *E_0 = &E_01;
    const double *Iu_0 = &Iu_01;

    for (u = 0; u < U; u++) {
      E[u] = 0;
      Iu[u] = 0;
      Ir[u] = 0;
      Ca[u] = 0;
      Cb[u] = 0;
      C[u] = 0;
      S[u] = nearbyint(pop[u]);
    }
      
    E[0]=2000*E_0[0];
    Iu[0]=2000*Iu_0[0];
    for (u = 1; u < U; u++) {
      E[u]=3*mob[14*u][0]*2000*E_0[u*E_0_expand]/(pop[0]);
      Iu[u]=3*mob[14*u][0]*2000*Iu_0[u*Iu_0_expand]/(pop[0]);
    }
  ")
  
  po <- spatPomp(covid_cases,
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
  
  if (measurement == "simulation"){
    po <- spatPomp(po, dmeasure=covid_dmeasure_simu,
                   vunit_measure=covid_vunit_measure_simu,
                   rmeasure=covid_rmeasure_simu,
                   dunit_measure=covid_dunit_measure_simu
    )
  }
  
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

