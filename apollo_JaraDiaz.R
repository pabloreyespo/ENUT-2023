
apollo_jaradiaz <- function(jaradiaz_settings, functionality="estimate"){
  # Rename input if necessary
  apollo_inputs <- tryCatch(get('apollo_inputs', envir=parent.frame(), inherits=FALSE),
                            error=function(e) list(silent=FALSE))

  # Copy variables from list to environment
  for(i in 1:length(jaradiaz_settings)) assign(names(jaradiaz_settings)[i], jaradiaz_settings[[i]])
  N <- nrow(apollo_inputs$database)

  if(functionality=="preprocess"){
    preproc_settings <- list(componentName = "jaraDiaz", gradient = TRUE)
    if(!is.null(jaradiaz_settings$componentName)){
      preproc_settings$componentName <- jaradiaz_settings$componentName
    } else if(!is.null(jaradiaz_settings$componentName2)){
      preproc_settings$componentName <- jaradiaz_settings$componentName2
    }
    return(preproc_settings)
  }

  # -------------- #
  #### VALIDATE ####
  # -------------- #
  if(functionality %in% c("validate")){ # (Engañito)
    return(invisible( rep(1, N) ))
  }

  # ------------- #
  #### ZERO LL ####
  # ------------- #
  if(functionality=="zero_LL"){
    ans <- rep(NA, N)
    return(ans)
  }

  # ------------------------------------ #
  #### ESTIMATE, CONDITIONALS AND RAW #### ----> Calculan la verosimilitud (La que tengo en likelihoods)
  # ------------------------------------ #
  if(functionality %in% c("estimate", "conditionals", "raw")){
    ec = Ec / (w * (ta-Tc))
    if ("Tw" %in% names(jaradiaz_settings)) {
      if ("alpha" %in% names(jaradiaz_settings)) {
        topt_work = (ta - Tc) *
           ( (beta + alpha*ec) + sqrt((beta + alpha*ec)^2 - (2*alpha+2*beta-1)*ec) )
        mu = (Tw - topt_work) /sigma
        ll = -0.5*mu^2 -log(sigma) #+ -0.5*log(2*pi)
        if(is.matrix(ll)) ll <- rowSums(ll)
        L <- exp(ll)
      } else {
        topt_work = ((ta-Tc) / (2*((PH + TH + theta_w)))) *
          (PH + theta_w + (TH + theta_w)*ec + sqrt((PH + theta_w + (TH + theta_w)*ec)^2 - 4*theta_w*ec*((PH + TH + theta_w))))
        mu = (Tw - topt_work) / sigma
        ll = -0.5*mu^2 -log(sigma) #+ -0.5*log(2*pi)
        if(is.matrix(ll)) ll <- rowSums(ll)
        L <- exp(ll)
        if (any(PH<0)) {L <- L*0.0001}
      }
    } else {
      if (any(is.na(corr))) corr = diag(nequations)
      Tc = obs_times[,ncol(obs_times)]
      Ec = obs_expenditures[,ncol(obs_expenditures)]

      if (nexpenditures > 0) { allocations = aux = as.matrix(cbind(obs_times[,1:ntimes, drop = F], obs_expenditures[,1:nexpenditures, drop = F]))
      } else { allocations = aux = as.matrix(obs_times[,1:ntimes, drop = F]) }
      aux[T] = 0
      opt_allocations    = aux
      conditional_mu     = aux
      conditional_sd     = rep(1, nequations)

      if ("alpha" %in% names(jaradiaz_settings)) {
        topt_work = (ta - Tc) *
           ( (beta + alpha*ec) + sqrt((beta + alpha*ec)^2 - (2*alpha+2*beta-1)*ec) )
        opt_allocations[,1] = topt_work

        if (ntimes >= 2) opt_allocations[,2:ntimes] = sapply(gammas, function(x) (x / (1-2*beta)) * (ta - topt_work - Tc))
        if (nexpenditures >= 1) opt_allocations[,(ntimes+1):nequations] = sapply(deltas, function(x) (x / (1-2*alpha)) * (w * topt_work - Ec))
        sd_errors = sweep(allocations - opt_allocations, MARGIN = 2, sigma, "/")
      } else {
        topt_work = ((ta-Tc) / (2*((PH + TH + theta_w)))) *
          (PH + theta_w + (TH + theta_w)*ec + sqrt((PH + theta_w + (TH + theta_w)*ec)^2 - 4*theta_w*ec*((PH + TH + theta_w))))
        opt_allocations[,1] = topt_work
        if (ntimes >= 2) opt_allocations[,2:ntimes] = sapply(thetas/TH, function(x) (x*(ta - topt_work - Tc)))
        if (nexpenditures >= 1) opt_allocations[,(ntimes+1):nequations] = sapply(phis/PH, function(x) (x * (w * topt_work - Ec)))
        sd_errors = sweep(allocations - opt_allocations, MARGIN = 2, sigma, "/")
      }

      if (!all(corr == diag(nequations))) {
        if (nequations > 1) for (j in 2:(nequations)) {
          if (j ==2) {
            conditional_mu[,2] = corr[2,1] * sd_errors[,1]
            conditional_sd[2]  = 1 - corr[2,1]^2
          } else if (j ==3) {
            conditional_mu[,3] = ((corr[2,3]-corr[1,3]*corr[1,2])* sd_errors[,2] + (corr[1,3]-corr[2,3]*corr[1,2])*sd_errors[,1]) / conditional_sd[2]
            conditional_sd[3]  = 1- (corr[2,3]^2 -2*corr[2,1]*corr[2,3]*corr[1,3] +corr[1,3]^2) / conditional_sd[2]
          } else {
            i = j-1
            conditional_mu[,j] = c(corr[j,i:1] %*% solve(corr[i:1,i:1]) %*% t(sd_errors[,i:1])) #(3|2,1)
            conditional_sd[j]  = c(corr[j,j] - corr[j,i:1] %*% solve(corr[i:1,i:1]) %*% corr[i:1,j]) }}
        conditional_errors = sweep(sd_errors - conditional_mu, MARGIN = 2, sqrt(conditional_sd), "/")
      } else conditional_errors = sd_errors

      lognorm = -0.5*(conditional_errors^2)
      ll = sweep(lognorm ,  MARGIN = 2, log(sigma * sqrt(conditional_sd)), "-")
      #invCorr <- solve(corr)
      #ll = -0.5*(log(det(corr)) + apply(sd_errors,1, function(x) t(x) %*% invCorr %*% x)) - sum(log(sigma))
      if(is.matrix(ll)) ll <- rowSums(ll)
      L <- exp(ll)
    }

    return(L)
  }

  # ------------ #
  #### OUTPUT #### ---> En base a los datos que entrego me hace un reporte, puedo hacerlo pero aun no es necesario
  # ------------ #
  if(functionality %in% c("output", "report")){
    ans <- apollo_jaradiaz(jaradiaz_settings, functionality="estimate")
    return(ans)
  }

  # ---------------- #
  #### PREDICTION ####   ---> Cómo predecir con el modelo
  # ---------------- #
  if(functionality=="prediction"){
    return(NA)
  }

  if(functionality=="gradient"){
    ec = Ec / (w * (ta-Tc))
    if ("Tw" %in% names(jaradiaz_settings)) {
      if ("alpha" %in% names(jaradiaz_settings)) {
        betaalphaec = (beta + alpha*ec)
        aux_sqrt = sqrt(betaalphaec^2 - (2*alpha+2*beta-1)*ec)
        topt_work = (ta - Tc) * (betaalphaec + aux_sqrt)
        mu = (Tw - topt_work)
        ll = -0.5*(mu/sigma)^2 -log(sigma) #+ -0.5*log(2*pi)
        if(is.matrix(ll)) ll <- rowSums(ll)
        L <- exp(ll)
        MuparcAlpha = -ec*(ta-topt_work-Tc)/ aux_sqrt
        MuparcBeta  = (topt_work - ec*(ta-Tc))/ aux_sqrt
        LLparcAlpha = (mu / sigma^2) * MuparcAlpha
        LLparcBeta  = (mu / sigma^2) * MuparcBeta
        LLparcSigma = (1/sigma) * ((mu/sigma)^2 - 1)

        G <- list()
        G[["alpha"]] =  LLparcAlpha * L
        G[["beta" ]] =  LLparcBeta  * L
        G[["sigma"]] =  LLparcSigma * L
      } else {
        x <- (PH + TH + theta_w)
        thetaphiec = PH + theta_w + (TH + theta_w)*ec
        aux_sqrt    = sqrt(thetaphiec^2 - 4*theta_w*ec*x)
        topt_work = (ta-Tc)*(thetaphiec + aux_sqrt) / (2*x)
        mu = (Tw - topt_work)
        ll = -0.5*(mu/sigma)^2 -log(sigma) #+ -0.5*log(2*pi)
        if(is.matrix(ll)) ll <- rowSums(ll)
        L <- exp(ll)

        MuparcPH   = (topt_work*(x-aux_sqrt) - theta_w*ec*(ta-Tc)) / (x *aux_sqrt)
        MuparcTHW  = (topt_work*(x*(1+ec)-aux_sqrt) - (PH + TH + 2*theta_w)*ec*(ta-Tc)) / (x *aux_sqrt)
        LLparcPH = (mu / sigma^2) * MuparcPH
        LLparcTHW  = (mu / sigma^2) * MuparcTHW
        LLparcSigma = (1/sigma) * ((mu/sigma)^2 - 1)

        G <- list()
        G[["PH"]]      =  LLparcPH * L
        G[["theta_w"]] =  LLparcTHW  * L
        G[["sigma"]]   =  LLparcSigma * L
      }
      output = list(like = L, grad = G)
    } else {
      output = list(like = NA, grad = NA)
    }
    return (output)
  }

  # End of function
  stop("Invalid value of argument 'functionality'")
}


#####################################################################################################
#####################################################################################################

apollo_jaradiaz_2pi <- function(jaradiaz_settings, functionality="estimate"){
  # Rename input if necessary
  apollo_inputs <- tryCatch(get('apollo_inputs', envir=parent.frame(), inherits=FALSE),
                            error=function(e) list(silent=FALSE))

  # Copy variables from list to environment
  for(i in 1:length(jaradiaz_settings)) assign(names(jaradiaz_settings)[i], jaradiaz_settings[[i]])
  N <- nrow(apollo_inputs$database)

  if(functionality=="preprocess"){
    preproc_settings <- list(componentName = "jaraDiaz", gradient = TRUE)
    if(!is.null(jaradiaz_settings$componentName)){
      preproc_settings$componentName <- jaradiaz_settings$componentName
    } else if(!is.null(jaradiaz_settings$componentName2)){
      preproc_settings$componentName <- jaradiaz_settings$componentName2
    }
    return(preproc_settings)
  }

  # -------------- #
  #### VALIDATE ####
  # -------------- #
  if(functionality %in% c("validate")){ # (Engañito)
    return(invisible( rep(1, N) ))
  }

  # ------------- #
  #### ZERO LL ####
  # ------------- #
  if(functionality=="zero_LL"){
    ans <- rep(NA, N)
    return(ans)
  }

  # ------------------------------------ #
  #### ESTIMATE, CONDITIONALS AND RAW #### ----> Calculan la verosimilitud (La que tengo en likelihoods)
  # ------------------------------------ #
  if(functionality %in% c("estimate", "conditionals", "raw")){
    ec = Ec / (w * (ta-Tc))
    if ("Tw" %in% names(jaradiaz_settings)) {
      if ("alpha" %in% names(jaradiaz_settings)) {
        topt_work = (ta - Tc) *
           ( (beta + alpha*ec) + sqrt((beta + alpha*ec)^2 - (2*alpha+2*beta-1)*ec) )
        mu = (Tw - topt_work) /sigma
        ll = -0.5*mu^2 -log(sigma) -0.5*log(2*base::pi)
        if(is.matrix(ll)) ll <- rowSums(ll)
        L <- exp(ll)
      } else {
        topt_work = ((ta-Tc) / (2*((PH + TH + theta_w)))) *
          (PH + theta_w + (TH + theta_w)*ec + sqrt((PH + theta_w + (TH + theta_w)*ec)^2 - 4*theta_w*ec*((PH + TH + theta_w))))
        mu = (Tw - topt_work) / sigma
        ll = -0.5*mu^2 -log(sigma) -0.5*log(2*base::pi)
        if(is.matrix(ll)) ll <- rowSums(ll)
        L <- exp(ll)
        if (any(PH<0)) {L <- L*0.0001}
      }
    } else {
      if (any(is.na(corr))) corr = diag(nequations)
      Tc = obs_times[,ncol(obs_times)]
      Ec = obs_expenditures[,ncol(obs_expenditures)]

      if (nexpenditures > 0) { allocations = aux = as.matrix(cbind(obs_times[,1:ntimes, drop = F], obs_expenditures[,1:nexpenditures, drop = F]))
      } else { allocations = aux = as.matrix(obs_times[,1:ntimes, drop = F]) }
      aux[T] = 0
      opt_allocations    = aux
      conditional_mu     = aux
      conditional_sd     = rep(1, nequations)

      if ("alpha" %in% names(jaradiaz_settings)) {
        topt_work = (ta - Tc) *
           ( (beta + alpha*ec) + sqrt((beta + alpha*ec)^2 - (2*alpha+2*beta-1)*ec) )
        opt_allocations[,1] = topt_work

        if (ntimes >= 2) opt_allocations[,2:ntimes] = sapply(gammas, function(x) (x / (1-2*beta)) * (ta - topt_work - Tc))
        if (nexpenditures >= 1) opt_allocations[,(ntimes+1):nequations] = sapply(deltas, function(x) (x / (1-2*alpha)) * (w * topt_work - Ec))
        sd_errors = sweep(allocations - opt_allocations, MARGIN = 2, sigma, "/")
      } else {
        topt_work = ((ta-Tc) / (2*((PH + TH + theta_w)))) *
          (PH + theta_w + (TH + theta_w)*ec + sqrt((PH + theta_w + (TH + theta_w)*ec)^2 - 4*theta_w*ec*((PH + TH + theta_w))))
        opt_allocations[,1] = topt_work
        if (ntimes >= 2) opt_allocations[,2:ntimes] = sapply(thetas/TH, function(x) (x*(ta - topt_work - Tc)))
        if (nexpenditures >= 1) opt_allocations[,(ntimes+1):nequations] = sapply(phis/PH, function(x) (x * (w * topt_work - Ec)))
        sd_errors = sweep(allocations - opt_allocations, MARGIN = 2, sigma, "/")
      }

      if (!all(corr == diag(nequations))) {
        if (nequations > 1) for (j in 2:(nequations)) {
          if (j ==2) {
            conditional_mu[,2] = corr[2,1] * sd_errors[,1]
            conditional_sd[2]  = 1 - corr[2,1]^2
          } else if (j ==3) {
            conditional_mu[,3] = ((corr[2,3]-corr[1,3]*corr[1,2])* sd_errors[,2] + (corr[1,3]-corr[2,3]*corr[1,2])*sd_errors[,1]) / conditional_sd[2]
            conditional_sd[3]  = 1- (corr[2,3]^2 -2*corr[2,1]*corr[2,3]*corr[1,3] +corr[1,3]^2) / conditional_sd[2]
          } else {
            i = j-1
            conditional_mu[,j] = c(corr[j,i:1] %*% solve(corr[i:1,i:1]) %*% t(sd_errors[,i:1])) #(3|2,1)
            conditional_sd[j]  = c(corr[j,j] - corr[j,i:1] %*% solve(corr[i:1,i:1]) %*% corr[i:1,j]) }}
        conditional_errors = sweep(sd_errors - conditional_mu, MARGIN = 2, sqrt(conditional_sd), "/")
      } else conditional_errors = sd_errors

      lognorm = -0.5*(conditional_errors^2)
      ll = sweep(lognorm ,  MARGIN = 2, log(sigma * sqrt(conditional_sd)), "-") - 0.5*log(2*base::pi)
      #invCorr <- solve(corr)
      #ll = -0.5*(log(det(corr)) + apply(sd_errors,1, function(x) t(x) %*% invCorr %*% x)) - sum(log(sigma))
      if(is.matrix(ll)) ll <- rowSums(ll)
      L <- exp(ll)
    }

    return(L)
  }

  # ------------ #
  #### OUTPUT #### ---> En base a los datos que entrego me hace un reporte, puedo hacerlo pero aun no es necesario
  # ------------ #
  if(functionality %in% c("output", "report")){
    ans <- apollo_jaradiaz(jaradiaz_settings, functionality="estimate")
    return(ans)
  }

  # ---------------- #
  #### PREDICTION ####   ---> Cómo predecir con el modelo
  # ---------------- #
  if(functionality=="prediction"){
    return(NA)
  }

  if(functionality=="gradient"){
    ec = Ec / (w * (ta-Tc))
    if ("Tw" %in% names(jaradiaz_settings)) {
      if ("alpha" %in% names(jaradiaz_settings)) {
        betaalphaec = (beta + alpha*ec)
        aux_sqrt = sqrt(betaalphaec^2 - (2*alpha+2*beta-1)*ec)
        topt_work = (ta - Tc) * (betaalphaec + aux_sqrt)
        mu = (Tw - topt_work)
        ll = -0.5*(mu/sigma)^2 -log(sigma) #+ -0.5*log(2*pi)
        if(is.matrix(ll)) ll <- rowSums(ll)
        L <- exp(ll)
        MuparcAlpha = -ec*(ta-topt_work-Tc)/ aux_sqrt
        MuparcBeta  = (topt_work - ec*(ta-Tc))/ aux_sqrt
        LLparcAlpha = (mu / sigma^2) * MuparcAlpha
        LLparcBeta  = (mu / sigma^2) * MuparcBeta
        LLparcSigma = (1/sigma) * ((mu/sigma)^2 - 1)

        G <- list()
        G[["alpha"]] =  LLparcAlpha * L
        G[["beta" ]] =  LLparcBeta  * L
        G[["sigma"]] =  LLparcSigma * L
      } else {
        x <- (PH + TH + theta_w)
        thetaphiec = PH + theta_w + (TH + theta_w)*ec
        aux_sqrt    = sqrt(thetaphiec^2 - 4*theta_w*ec*x)
        topt_work = (ta-Tc)*(thetaphiec + aux_sqrt) / (2*x)
        mu = (Tw - topt_work)
        ll = -0.5*(mu/sigma)^2 -log(sigma) #+ -0.5*log(2*pi)
        if(is.matrix(ll)) ll <- rowSums(ll)
        L <- exp(ll)

        MuparcPH   = (topt_work*(x-aux_sqrt) - theta_w*ec*(ta-Tc)) / (x *aux_sqrt)
        MuparcTHW  = (topt_work*(x*(1+ec)-aux_sqrt) - (PH + TH + 2*theta_w)*ec*(ta-Tc)) / (x *aux_sqrt)
        LLparcPH = (mu / sigma^2) * MuparcPH
        LLparcTHW  = (mu / sigma^2) * MuparcTHW
        LLparcSigma = (1/sigma) * ((mu/sigma)^2 - 1)

        G <- list()
        G[["PH"]]      =  LLparcPH * L
        G[["theta_w"]] =  LLparcTHW  * L
        G[["sigma"]]   =  LLparcSigma * L
      }
      output = list(like = L, grad = G)
    } else {
      output = list(like = NA, grad = NA)
    }
    return (output)
  }

  # End of function
  stop("Invalid value of argument 'functionality'")
}