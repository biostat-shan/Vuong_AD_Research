
source('cov_Donanemab Phase 3.R')

library(mvtnorm)
library(parallel)
library(nlme)

set.seed(20262026)

Month <- c(3,6,9,12,15,18)
Month.followup <- Month[Month != 0]
followupN <- length(Month.followup)

Narm   <- c(50,100,150,200,250,300,350,400)
NsimSD <- 1000  # number of simulations per (shape, N)

## -------------------- helper: stable inverse --------------------
safe_solve <- function(S, ridge = 1e-8){
  S <- as.matrix(S)
  K <- ncol(S)
  solve(S + diag(ridge, K))
}

## -------------------- your original helper functions --------------------
MatchPointA <- function(Month,Tv,Pv) {
  Monthv <- c(Month)
  if (sum(Pv > 0) > sum(Pv < 0)) {
    Pv <- c(Pv); Tv <- c(Tv)
  } else {
    Pv <- c(-Pv); Tv <- c(-Tv)
  }
  indexCeiling <- length(Monthv)
  Index <- sum(Pv < Tv[indexCeiling])
  if(Index > indexCeiling - 0.5){
    Monthv[indexCeiling+1] <- Month[indexCeiling]*2
    Pv[indexCeiling+1] <- Pv[indexCeiling] + (Pv[indexCeiling]-Pv[indexCeiling-1]) * Month[indexCeiling] / (Month[indexCeiling]-Month[indexCeiling-1])
    slope <- (Pv[Index+1] - Pv[Index]) / (Monthv[Index+1] - Monthv[Index])
    cross.x <- (Tv[indexCeiling]-Pv[Index])/slope + Monthv[Index]
  } else if (Index == 0) {
    slope <- (Pv[length(Monthv)] - Pv[1]) / (Monthv[length(Monthv)] - Monthv[1])
    cross.x <- Tv[indexCeiling]/slope
  } else {
    slope <- (Pv[Index+1] - Pv[Index]) / (Monthv[Index+1] - Monthv[Index])
    cross.x <- (Tv[indexCeiling]-Pv[Index])/slope + Monthv[Index]
  }
  cross.x
}

adjust_mean_actA <- function(mean_actA, mean_pbo) {
  if (sum(mean_actA >0) > sum(mean_actA<0) & sum(mean_actA<0)>0.5 ) {
    NegIndex <- max(which(mean_actA<0))
    mean_actA[c(1:(NegIndex+1))] <- seq(min(mean_pbo[1]/2, mean_actA[NegIndex+1]/2), mean_actA[NegIndex+1], length.out=(NegIndex+1))
  }
  if (sum(mean_actA >0) < sum(mean_actA<0) & sum(mean_actA>0)>0.5 ) {
    PosiIndex <- max(which(mean_actA>0))
    mean_actA[c(1:(PosiIndex+1))] <- seq(min(mean_pbo[1]/2, mean_actA[PosiIndex+1]/2), mean_actA[PosiIndex+1], length.out=(PosiIndex+1))
  }
  if (sum(mean_actA >0) == sum(mean_actA<0) & tail(mean_pbo,1)>0 ) {
    PosiIndex <- min(which(mean_actA>0))
    mean_actA[1:PosiIndex] <- seq(mean_pbo[1]/2, min(mean_actA[PosiIndex],mean_pbo[PosiIndex]), length.out=PosiIndex)
  }
  if (sum(mean_actA >0) == sum(mean_actA<0) & tail(mean_pbo,1)<0 ) {
    NegIndex <- min(which(mean_actA<0))
    mean_actA[1:NegIndex] <- seq(mean_pbo[1]/2, max(mean_actA[NegIndex],mean_pbo[NegIndex]), length.out=NegIndex)
  }
  mean_actA
}

## -------------------- objects from your sourced file --------------------
DRUG <- "iADRS"

seReal   <- get(paste0("cov_", DRUG))        # your "standard error" object
mean_pbo <- get(paste0(DRUG, ".Placebo"))    # placebo mean vector

sePlace <- seReal
seTreat <- seReal

Mean.Placebo <- mean_pbo
MeanTrue.placeReal <- Mean.Placebo

## keep your scaling (as you used)
covPlaceReal <- sePlace * 600
covTreatReal <- seTreat * 600

## -------------------- build shapes (your original style) --------------------
TRTshapeA <- data.frame()
TRTshapeB <- data.frame()
TRTshapeC <- data.frame()
TRTshapeD <- data.frame()

## Shape D: delay (x-axis shift via spline)
Dpara <- c(0.92, 0.90, 0.88, 0.86, 0.84, 0.82, 0.80)

define_mu_treatment_solver <- function(t, mu_placebo, t_c) {
  mu_placebo_func <- splinefun(t, mu_placebo, method = "natural")
  mu_treatment <- function(s) mu_placebo_func(s * t_c)
  return(mu_treatment)
}

for(Cdelay in Dpara){
  mean_actD_fun <- define_mu_treatment_solver(Month, Mean.Placebo, t_c = Cdelay)
  mean_actD <- mean_actD_fun(Month)
  TRTshapeD <- rbind(TRTshapeD, mean_actD)
}
TRTshapeD$shape <- Dpara
TRTshapeD <- TRTshapeD[, c(7, 1:length(Month))]  # shape + 6 visits

## Shape C: y-direction reduction (kept)
rate <- 1 - seq(0.08, 0.2, by = 0.02)
for(cvalueij in 1:7){
  mean_actC <- Mean.Placebo * rate[cvalueij]
  TRTshapeC <- rbind(TRTshapeC, mean_actC)
}
TRTshapeC$shape <- rate
TRTshapeC <- TRTshapeC[, c(7, 1:6)]

define_mu_treatment_shift <- function(t, mu_placebo, td, clamp = FALSE) {
  fp <- splinefun(t, mu_placebo, method = "natural")  # F_p(·)
  
  if(!clamp){
    
    ft <- function(s) fp(s + td)
  } else {
    
    tmin <- min(t); tmax <- max(t)
    ft <- function(s) fp(pmin(pmax(s + td, tmin), tmax))
  }
  return(ft)
}

TDpara <- c(-0.4, -0.8, -1.2, -1.6, -2, -2.4)  

TRTshape_shift <- data.frame()
for(td in TDpara){
  mean_act_fun <- define_mu_treatment_shift(Month, Mean.Placebo, td = td, clamp = FALSE)
  mean_act     <- mean_act_fun(Month)
  TRTshape_shift <- rbind(TRTshape_shift, mean_act)
}
TRTshape_shift$shape <- TDpara
TRTshape_shift <- TRTshape_shift[, c(ncol(TRTshape_shift), 1:length(Month))]  # shape + 6 visits


mix_chisq_q_cached <- local({
  cache <- list()
  function(alpha, w, nsim = 50000L, seed = 42L) {
    key <- paste(round(w, 4), collapse = "_")
    if (!is.null(cache[[key]])) return(cache[[key]])
    
    if (exists(".Random.seed", envir = .GlobalEnv))
      old_seed <- .Random.seed
    else
      old_seed <- NULL
    
    set.seed(seed)
    x <- w[1] * rchisq(nsim, df=1) + w[2] * rchisq(nsim, df=1)
    
    if (!is.null(old_seed))
      .Random.seed <<- old_seed
    
    val <- unname(quantile(x, probs = 1 - alpha, type = 8))
    cache[[key]] <<- val
    val
  }
})



shapeTYPE <- TRTshapeD

## -------------------- output container (final dataframe) --------------------
vuong_results <- data.frame()
overlap_results <- data.frame()

dec_levels  <- c("CDF","CDT")
type_levels <- c("Choice","Recommend")

SHAPEpara = 3
i = 3
## ============================================================
## MAIN LOOPS (same structure as yours)
## ============================================================
for(SHAPEpara in 1:nrow(shapeTYPE)){
  
  ## keep your indexing style: shapeTYPE has shape in col1, visits in col2..col7
  mean_act <- as.numeric(shapeTYPE[SHAPEpara, 2:(followupN+1)])
  
  for(i in 1:length(Narm)){
    
    Narm0 <- Narm[i]
    Narm1 <- Narm[i]
    N_current <- Narm0
    
    parallel.Vuong <- function(fsimu){
      
      ## ----------------------------------------------------------
      ## 1) Generate simulated data
      ## IMPORTANT FIX: treatment should use covTreatReal (you used covPlaceReal before)
      ## ----------------------------------------------------------
      dataX <- rmvnorm(n = Narm0, mean = MeanTrue.placeReal, sigma = covPlaceReal)
      dataY <- rmvnorm(n = Narm1, mean = mean_act,            sigma = covTreatReal)
      
      Mean.Placebo <- colMeans(dataX)
      Mean.Treat   <- colMeans(dataY)
      
      covPsimu <- cov(dataX)
      covTsimu <- cov(dataY)
      
      np <- Narm0
      nt <- Narm1
      K  <- length(Month)
      
      ## stable inverses (avoid singular covariance crash)
      Sigma_p_inv <- safe_solve(covPsimu, ridge = 1e-8)
      Sigma_t_inv <- safe_solve(covTsimu, ridge = 1e-8)
      
      Y_bar_p <- matrix(Mean.Placebo, ncol = 1)
      Y_bar_t <- matrix(Mean.Treat,   ncol = 1)
      
      ## ==========================================================
      ## 2) PP family (keep your grid + optim logic)
      ## ==========================================================
      #logLRmin <- function(lambda){
      
      #partA <- np*Sigma_p_inv + nt*(lambda^2)*Sigma_t_inv
      #partB <- np*Sigma_p_inv %*% Y_bar_p + nt*lambda*Sigma_t_inv %*% Y_bar_t
      
      #Fv <- solve(partA, partB)
      
      #r1 <- Y_bar_p - Fv
      #r2 <- Y_bar_t - lambda*Fv
      
      #q1 <- np/2 * as.numeric(t(r1) %*% Sigma_p_inv %*% r1)
      #q2 <- nt/2 * as.numeric(t(r2) %*% Sigma_t_inv %*% r2)
      #(q1 + q2)
      #}
      
      #lambda_seq_pp <- seq(0, 3, by = 0.01)
      #logLR_seq_pp  <- sapply(lambda_seq_pp, logLRmin)
      #lambda_hat_pp <- lambda_seq_pp[ which.min(logLR_seq_pp) ]
      #lambda_hat_pp <- optim(lambda_hat_pp, logLRmin, method="L-BFGS-B", lower=0, upper=3)$par
      
      #partA <- np*Sigma_p_inv + nt*(lambda_hat_pp^2)*Sigma_t_inv
      #partB <- np*Sigma_p_inv %*% Y_bar_p + nt*lambda_hat_pp*Sigma_t_inv %*% Y_bar_t
      #Fp_hat_pp <- drop(solve(partA, partB))
      
      #logLR_PP <- function(Fp, lambda){
      #Fp <- matrix(Fp, ncol=1)
      #r1 <- Y_bar_p - Fp
      #r2 <- Y_bar_t - lambda*Fp
      #q1 <- np/2 * as.numeric(t(r1)%*%Sigma_p_inv%*%r1)
      #q2 <- nt/2 * as.numeric(t(r2)%*%Sigma_t_inv%*%r2)
      #-(q1 + q2)
      #}
      
      ## ==========================================================
      ## 3) CDF family (your sign-flip style, but with fallback instead of stop)
      ## ==========================================================
      generate_Gtd <- function(K, td) {
        G <- matrix(0, nrow = K, ncol = K)
        for (k in 1:K) {
          k_td <- k * td
          if (k_td >= 1 && k_td <= K) {
            t_bar <- floor(k_td)
            t_star <- k_td - t_bar
            if (t_bar + 1 <= K) G[k, t_bar + 1] <- t_star
            if (t_bar >= 1)     G[k, t_bar]     <- 1 - t_star
          } else if (k_td > K) {
            t_bar <- K-1
            t_star <- k_td - t_bar
            G[k, t_bar]   <- 1 - t_star
            G[k, t_bar+1] <- t_star
          } else {
            t_bar <- 1
            floor_k_td <- floor(k_td)
            t_star <- k_td - floor_k_td
            G[k, t_bar]     <- 1 - (floor_k_td + t_star - 1)
            G[k, t_bar + 1] <- (floor_k_td + t_star) - 1
          }
        }
        G
      }
      
      generate_Gbar_td <- function(K, td) {
        Gbar <- matrix(0, nrow = K, ncol = K)
        if(abs(td-1) < 1e-8) td <- 0.99999
        for(k in 1:K){
          k_td <- k * td
          if(k_td >= 1 && k_td <= K){
            t_bar <- floor(k_td)
            if(t_bar + 1 <= K) Gbar[k, t_bar + 1] <-  k
            if(t_bar >= 1)     Gbar[k, t_bar]     <- -k
          } else if (k_td > K){
            t_bar <- K-1
            Gbar[k, t_bar + 1] <-  k
            Gbar[k, t_bar]     <- -k
          } else {
            t_bar <- 1
            Gbar[k, t_bar + 1] <-  k
            Gbar[k, t_bar]     <- -k
          }
        }
        Gbar
      }
      
      max_iter <- 300
      tol <- 1e-4
      diff <- 1
      iter <- 0
      
      Fp_hat <- Y_bar_p
      td_hat_current <- 2.998
      
      posschange1 <- sort(1 - c((1:5)/6, (1:4)/5, 1/4, 3/4))
      posschange  <- c(posschange1, 0.999999, 1.0001, 1.566, 2.998)
      
      score_at <- function(td, Fp){
        Gtd  <- generate_Gtd(K, td)
        Gbar <- generate_Gbar_td(K, td)
        as.numeric(nt * t(Fp) %*% t(Gbar) %*% Sigma_t_inv %*% (Y_bar_t - Gtd %*% Fp))
      }
      
      while(diff > tol && iter < max_iter){
        iter <- iter + 1
        
        sc <- vapply(posschange, function(td) score_at(td, Fp_hat), numeric(1))
        flip_index <- which(diff(sign(sc)) != 0)[1]
        
        if(is.na(flip_index)){
          ## Fallback: profile search on a grid (avoid stopping)
          td_grid <- seq(min(posschange), max(posschange), by = 0.002)
          ll_grid <- vapply(td_grid, function(td){
            Gtd <- generate_Gtd(K, td)
            A <- (np/nt)*Sigma_p_inv + t(Gtd)%*%Sigma_t_inv%*%Gtd
            b <- (np/nt)*Sigma_p_inv%*%Y_bar_p + t(Gtd)%*%Sigma_t_inv%*%Y_bar_t
            Fp_tmp <- solve(A, b)
            r1 <- Y_bar_p - Fp_tmp
            r2 <- Y_bar_t - Gtd%*%Fp_tmp
            q1 <- (np/2)*as.numeric(t(r1)%*%Sigma_p_inv%*%r1)
            q2 <- (nt/2)*as.numeric(t(r2)%*%Sigma_t_inv%*%r2)
            -(q1+q2)
          }, numeric(1))
          td_hat_new <- td_grid[ which.max(ll_grid) ]
        } else {
          lo <- posschange[flip_index]
          hi <- posschange[flip_index+1]
          
          Dtd72 <- function(td_hat_est){
            Gtd_est  <- generate_Gtd(K, td_hat_est)
            Gbar_est <- generate_Gbar_td(K, td_hat_est)
            v <- t(Fp_hat) %*% t(Gbar_est) %*% Sigma_t_inv %*% (Y_bar_t - Gtd_est %*% Fp_hat)
            as.numeric(v^4)
          }
          Find72 <- optimise(Dtd72, lower=lo, upper=hi)
          td_hat_new <- as.numeric(Find72$minimum)
        }
        
        Gtd_current <- generate_Gtd(K, td_hat_new)
        A <- (np/nt)*Sigma_p_inv + t(Gtd_current)%*%Sigma_t_inv%*%Gtd_current
        b <- (np/nt)*Sigma_p_inv%*%Y_bar_p + t(Gtd_current)%*%Sigma_t_inv%*%Y_bar_t
        Fp_hat_new <- solve(A, b)
        
        diff <- abs(td_hat_new - td_hat_current)
        td_hat_current <- td_hat_new
        Fp_hat <- Fp_hat_new
      }
      
      Fp_hat_cdf <- drop(Fp_hat)
      lambda_hat_cdf <- td_hat_current
      
      logLR_CDF <- function(Fp, td){
        Fp <- matrix(Fp, ncol=1)
        Gtd <- generate_Gtd(K, td)
        r1 <- Y_bar_p - Fp
        r2 <- Y_bar_t - Gtd%*%Fp
        q1 <- (np/2)*as.numeric(t(r1)%*%Sigma_p_inv%*%r1)
        q2 <- (nt/2)*as.numeric(t(r2)%*%Sigma_t_inv%*%r2)
        -(q1+q2)
      }
      
      ## ==========================================================
      ## 4) CDT family (WRITE INSIDE LOOP, eq(57)+eq(59), keep Gbar)
      ## ==========================================================
      ## ---- G_td for CDT: k_td = k + td ----
      generate_Gtd_CDT_local <- function(K, td) {
        G <- matrix(0, nrow = K, ncol = K)
        for (k in 1:K) {
          k_td <- k + td
          if (k_td >= 1 && k_td <= K) {
            t_bar <- floor(k_td)
            t_star <- k_td - t_bar
            if (t_bar + 1 <= K) G[k, t_bar + 1] <- t_star
            if (t_bar >= 1)     G[k, t_bar]     <- 1 - t_star
          } else if (k_td > K) {
            t_bar <- K-1
            t_star <- k_td - t_bar
            G[k, t_bar]   <- 1 - t_star
            G[k, t_bar+1] <- t_star
          } else {
            t_bar <- 1
            floor_k_td <- floor(k_td)
            t_star <- k_td - floor_k_td
            G[k, t_bar]     <- 1 - (floor_k_td + t_star - 1)
            G[k, t_bar + 1] <- (floor_k_td + t_star) - 1
          }
        }
        G
      }
      
      ## ---- Gbar_td for CDT: derivative wrt td (piecewise constant) ----
      generate_Gbar_td_CDT_local <- function(K, td) {
        Gbar <- matrix(0, nrow = K, ncol = K)
        if(abs(td) < 1e-10) td <- 1e-6  # avoid exact breakpoint
        for(k in 1:K){
          k_td <- k + td
          if(k_td >= 1 && k_td <= K){
            t_bar <- floor(k_td)
            if(t_bar + 1 <= K) Gbar[k, t_bar + 1] <-  1
            if(t_bar >= 1)     Gbar[k, t_bar]     <- -1
          } else if (k_td > K){
            t_bar <- K-1
            Gbar[k, t_bar + 1] <-  1
            Gbar[k, t_bar]     <- -1
          } else {
            t_bar <- 1
            Gbar[k, t_bar + 1] <-  1
            Gbar[k, t_bar]     <- -1
          }
        }
        Gbar
      }
      
      ## ---- eq(59): Fp_hat given td ----
      Fp_hat_given_td_CDT <- function(td){
        Gtd <- generate_Gtd_CDT_local(K, td)
        A <- (np/nt)*Sigma_p_inv + t(Gtd) %*% Sigma_t_inv %*% Gtd
        b <- (np/nt)*Sigma_p_inv %*% Y_bar_p + t(Gtd) %*% Sigma_t_inv %*% Y_bar_t
        solve(A, b)
      }
      
      ## ---- profile logLik(td): used only as fallback when eq(57) fails ----
      logLik_td_CDT <- function(td){
        Gtd <- generate_Gtd_CDT_local(K, td)
        Fp  <- Fp_hat_given_td_CDT(td)
        r1 <- Y_bar_p - Fp
        r2 <- Y_bar_t - Gtd %*% Fp
        q1 <- (np/2) * as.numeric(t(r1) %*% Sigma_p_inv %*% r1)
        q2 <- (nt/2) * as.numeric(t(r2) %*% Sigma_t_inv %*% r2)
        -(q1 + q2)
      }
      
      ## ---- CDT iteration using eq(57) and eq(59) ----
      max_iter_cdt <- 300
      tol_cdt <- 1e-4
      diff_cdt <- 1
      iter_cdt <- 0
      
      Fp_hat_cdt_mat <- Y_bar_p
      td_hat_current_cdt <- -1e-3
      
      m_grid <- 0:-10
      tiny_den <- 1e-12
      
      while(diff_cdt > tol_cdt && iter_cdt < max_iter_cdt){
        iter_cdt <- iter_cdt + 1
        
        ## Step 1: update td using eq(57) on each integer interval m
        td_hat_new_cdt <- NA_real_
        
        for(m in m_grid){
          Gm    <- generate_Gtd_CDT_local(K, td = m)
          GbarM <- generate_Gbar_td_CDT_local(K, td = m)
          
          num <- as.numeric(t(Fp_hat_cdt_mat) %*% t(GbarM) %*% Sigma_t_inv %*% (Y_bar_t - Gm %*% Fp_hat_cdt_mat))
          den <- as.numeric(t(Fp_hat_cdt_mat) %*% t(GbarM) %*% Sigma_t_inv %*% (GbarM %*% Fp_hat_cdt_mat))
          
          ## If denominator is too small, skip this m
          if(!is.finite(num) || !is.finite(den) || abs(den) < tiny_den) next
          
          td_star <- num / den   # eq(57)
          if(is.finite(td_star) && td_star >= 0 && td_star < 1){
            td_hat_new_cdt <- m + td_star
            break
          }
        }
        
        ## Fallback: if eq(57) cannot produce a valid td in any interval
        if(is.na(td_hat_new_cdt)){
          td_grid <- seq(min(m_grid), 0, by = 0.05)
          ll_grid <- vapply(td_grid, logLik_td_CDT, numeric(1))
          td0 <- td_grid[ which.max(ll_grid) ]
          
          ## refine within [floor(td0), floor(td0)+1)
          eps <- 1e-6
          m0 <- floor(td0)
          lo <- max(min(m_grid), m0 + eps)
          hi <- min(0, m0 + 1 - eps)
          if(lo < hi){
            opt <- optimize(function(td) -logLik_td_CDT(td), lower = lo, upper = hi)
            td_hat_new_cdt <- opt$minimum
          } else {
            td_hat_new_cdt <- td0
          }
        }
        
        ## Step 2: update Fp using eq(59)
        Gtd_current_cdt <- generate_Gtd_CDT_local(K, td = td_hat_new_cdt)
        A <- (np/nt) * Sigma_p_inv + t(Gtd_current_cdt) %*% Sigma_t_inv %*% Gtd_current_cdt
        b <- (np/nt) * Sigma_p_inv %*% Y_bar_p + t(Gtd_current_cdt) %*% Sigma_t_inv %*% Y_bar_t
        Fp_hat_new_cdt <- solve(A, b)  # eq(59)
        
        ## Step 3: check convergence on td
        diff_cdt <- abs(td_hat_new_cdt - td_hat_current_cdt)
        
        td_hat_current_cdt <- td_hat_new_cdt
        Fp_hat_cdt_mat <- Fp_hat_new_cdt
      }
      
      lambda_hat_cdt <- td_hat_current_cdt
      Fp_hat_cdt <- drop(Fp_hat_cdt_mat)
      
      logLR_CDT <- function(Fp, td){
        Fp <- matrix(Fp, ncol=1)
        Gtd <- generate_Gtd_CDT_local(K, td)
        r1 <- Y_bar_p - Fp
        r2 <- Y_bar_t - Gtd%*%Fp
        q1 <- (np/2)*as.numeric(t(r1)%*%Sigma_p_inv%*%r1)
        q2 <- (nt/2)*as.numeric(t(r2)%*%Sigma_t_inv%*%r2)
        -(q1+q2)
      }
      
      ## ==========================================================
      ## 5) Vuong statistics (keep your formula)
      ## ==========================================================
      Lambda_p_fun <- function(Fp1, Fp2, idx) {
        Y <- dataX[idx, ]
        d1 <- Y - Fp1
        d2 <- Y - Fp2
        as.numeric((-0.5)*t(d1)%*%Sigma_p_inv%*%d1 - (-0.5)*t(d2)%*%Sigma_p_inv%*%d2)
      }
      
      n_all <- np + nt
      
      ## CDF vs CDT
      lambda_t_fun_CDF_CDT <- function(Fp1, lambda1, Fp2, lambda2, idx) {
        Y <- dataY[idx, ]
        Gcdf <- generate_Gtd(K, lambda1)
        Gcdt <- generate_Gtd_CDT_local(K, lambda2)
        d1 <- Y - as.numeric(Gcdf %*% matrix(Fp1, ncol=1))
        d2 <- Y - as.numeric(Gcdt %*% matrix(Fp2, ncol=1))
        as.numeric((-0.5)*t(d1)%*%Sigma_t_inv%*%d1 - (-0.5)*t(d2)%*%Sigma_t_inv%*%d2)
      }
      
      Lambda_p <- vapply(1:np, function(ii) Lambda_p_fun(Fp_hat_cdf, Fp_hat_cdt, ii), numeric(1))
      Lambda_t <- vapply(1:nt, function(ii) lambda_t_fun_CDF_CDT(Fp_hat_cdf, lambda_hat_cdf, Fp_hat_cdt, lambda_hat_cdt, ii), numeric(1))
      LR_CDF_CDT <- (sum(Lambda_p) + sum(Lambda_t)) / n_all
      omega2 <- (sum((Lambda_p - LR_CDF_CDT)^2) + sum((Lambda_t - LR_CDF_CDT)^2)) / n_all
      omega <- sqrt(omega2)
      T_cdf_cdt <- ifelse(omega <= 1e-14, 0, sqrt(n_all) * LR_CDF_CDT / omega)
      nw2 <- n_all * omega2
      
      Fp_given_G <- function(G){
        M <- (np/nt) * Sigma_p_inv + t(G) %*% Sigma_t_inv %*% G
        b <- (np/nt) * Sigma_p_inv %*% Y_bar_p + t(G) %*% Sigma_t_inv %*% Y_bar_t
        solve(M, b)
      }
      
      dFp_dlam <- function(G, Gdot, Fp){
        M <- np * Sigma_p_inv + nt * t(G) %*% Sigma_t_inv %*% G 
        b <- 2 * nt * t(G) %*% Sigma_t_inv %*% Gdot %*% Fp
        -solve(M, b)
      }
      
      score_subjects <- function(dX, dY, Fp, dFp, G, Gdot){
        Ft_prime <- Gdot %*% Fp + G %*% dFp
        mu_t <- G %*% Fp
        Rp <- sweep(dX, 2, drop(Fp), "-")
        Rt <- sweep(dY, 2, mu_t, "-")
        sp       <- as.numeric(Rp %*% (Sigma_p_inv %*% dFp))
        st       <- as.numeric(Rt %*% (Sigma_t_inv %*% Ft_prime))
        list(sp = sp, st = st, Ft_prime = Ft_prime)
      }
      
      G0_cdf <- generate_Gtd(K, 1)
      Gb0_cdf <- generate_Gbar_td(K, 1)
      Fp0_cdf <- Fp_given_G(G0_cdf)
      dFp0_cdf <- dFp_dlam(G0_cdf, Gb0_cdf, Fp0_cdf)
      sc_cdf <- score_subjects(dataX, dataY, Fp0_cdf, dFp0_cdf, G0_cdf, Gb0_cdf)
      
      G0_cdt <- generate_Gtd_CDT_local(K, 0)
      Gb0_cdt <- generate_Gbar_td_CDT_local(K, 0)
      Fp0_cdt <- Fp_given_G(G0_cdt)
      dFp0_cdt  <- dFp_dlam(G0_cdt, Gb0_cdt, Fp0_cdt)
      sc_cdt    <- score_subjects(dataX, dataY, Fp0_cdt, dFp0_cdt, G0_cdt, Gb0_cdt)
      
      A11 <- -(np/n_all) * as.numeric(t(dFp0_cdf) %*% Sigma_p_inv %*% dFp0_cdf) +
        -(nt/n_all) * as.numeric(t(sc_cdf$Ft_prime) %*% Sigma_t_inv %*% sc_cdf$Ft_prime)
      
      A22 <- +(np/n_all) * as.numeric(t(dFp0_cdt) %*% Sigma_p_inv %*% dFp0_cdt) +
        +(nt/n_all) * as.numeric(t(sc_cdt$Ft_prime) %*% Sigma_t_inv %*% sc_cdt$Ft_prime)
      
      A_hat <- diag(c(A11, A22))
      
      Gmat  <- rbind(cbind(sc_cdf$sp, -sc_cdt$sp),
                     cbind(sc_cdf$st, -sc_cdt$st))
      B_hat <- crossprod(Gmat) / n_all
      
      ev <- Re(eigen(A_hat %*% B_hat, only.values = TRUE)$values)
      wa <- ev[1]; wb <- ev[2]
      if (!is.finite(wa) || !is.finite(wb)) { wa <- 0.5; wb <- 0.5 }
      
      crit_mix   <- mix_chisq_q_cached(0.05, c(wa, wb), nsim = 50000L, seed = 42L)
      reject_H00 <- (nw2 > crit_mix)
      
      
      decision <- data.frame(selected = "", decision_type = "")
      if(T_cdf_cdt >= 1.96) {
        decision$selected <- "CDF"
        decision$decision_type <- "Choice"
      }else if(T_cdf_cdt <= -1.96){
        decision$selected <- "CDT"
        decision$decision_type <- "Choice"
      }else if(T_cdf_cdt > 0 && T_cdf_cdt < 1.96 ){
        decision$selected <- "CDF"
        decision$decision_type <- "Recommend"
      }else if(T_cdf_cdt > -1.96 && T_cdf_cdt < 0) {
        decision$selected <- "CDT"
        decision$decision_type <- "Recommend"
      }
      
      c(decision$selected, decision$decision_type, T_cdf_cdt, reject_H00, LR_CDF_CDT, nw2, crit_mix)
    }
    
    ## run simulations (keep your structure)
    NsimSD=100
    results <- mclapply(1:NsimSD, function(x) parallel.Vuong(fsimu = x), mc.cores = 1)
    
    mat <- do.call(rbind, results)
    mat <- as.data.frame(mat)
    mat_T <- mat[mat$V4 == TRUE, ]
    T_results <- data.frame(Decision = mat_T$V1, Type = mat_T$V2)
    #T_results <- data.frame(Decision = mat[,1], Type = mat[,2], stringsAsFactors = FALSE)
    
    ## force levels so table includes zeros
    T_results$Decision <- factor(T_results$Decision, levels = dec_levels)
    T_results$Type     <- factor(T_results$Type,     levels = type_levels)
    
    tab <- table(T_results$Decision, T_results$Type) / NsimSD
    tab
    
    out <- as.data.frame(t(c(tab)))  
    names(out) <- c(outer(dec_levels, type_levels, paste, sep = "_"))
    out
    
    out$N <- N_current
    out$shape <- shapeTYPE$shape[SHAPEpara]
    
    out <- out[, c("shape", "N", "CDF_Choice", "CDF_Recommend",
                   "CDT_Choice", "CDT_Recommend")]
    
    mat_H0 <- mat[mat$V4 == FALSE,]
    H0_result <- mat_H0[,4:5]
    
    H0_result$Decision <- ifelse(H0_result$V5 >= 0, "CDF_Recommend", "CDT_Recommend")
    H0_result$Decision <- factor(H0_result$Decision, levels = c("CDF_Recommend", "CDT_Recommend"))
    tab_H0 <- table(H0_result$Decision)/NsimSD
    out_H0 <- as.data.frame(t(c(tab_H0)))
    
    out_H0$N <- N_current 
    out_H0$shape <- shapeTYPE$shape[SHAPEpara]
    
    out_H0 <- out_H0[, c("shape", "N", "CDF_Recommend", "CDT_Recommend")]
    overlap_results <- rbind(overlap_results, out_H0)
    
    ## append to final dataframe
    vuong_results <- rbind(vuong_results, out)
    
    cat("Done shape =", shapeTYPE$shape[SHAPEpara], " N =", N_current, "\n")
  }
}

## final output
vuong_results
overlap_results













