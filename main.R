library(vars)
## DGP
# AR(1) model parameters
sigma_ar1 <- diag(1,2);
phi_ar1 <- matrix(c(-1,0.96,-1.5,1.4), nrow = 2, ncol=2, byrow = T);

# AR(2) model parameters
sigma_ar2 <- matrix(c(1,-0.08,-0.08,1), nrow = 2, ncol=2, byrow = T);
phi_ar2 <- matrix(c(0.5,-0.3,0.2,0.65,-0.5,0.3,0,0.4), nrow = 4, ncol=2, byrow = T);

# sample size = 40, # of simulation = 5000
sample_size <- 40;
nsim <- 5000;
max_lag = 6;

# matrices to store selected model numbers
ar1_model_select <- matrix(NaN,ncol=4,nrow = nsim)
ar2_model_select <- matrix(NaN,ncol=4,nrow = nsim)
ar1_delta <- matrix(NaN,ncol=6,nrow = nsim )
ar2_delta <- matrix(NaN,ncol=6,nrow = nsim )

ar1_aic_select <- matrix(NaN,ncol=nsim,nrow = 6)
ar1_aicc_select <- matrix(NaN,ncol=nsim,nrow = 6)
ar1_aicbd_select <- matrix(NaN,ncol=nsim,nrow = 6)
ar1_sic_select <- matrix(NaN,ncol=nsim,nrow = 6)


ar2_aic_select <- matrix(NaN,ncol=nsim,nrow = 6)
ar2_aicc_select <- matrix(NaN,ncol=nsim,nrow = 6)
ar2_aicbd_select <- matrix(NaN,ncol=nsim,nrow = 6)
ar2_sic_select <- matrix(NaN,ncol=nsim,nrow = 6)

# Monte Carlo
for (jj in 1:nsim) {
  y_ar1 <- dgp(sample_size,sigma_ar1,phi_ar1);
  y_ar2 <- dgp(sample_size,sigma_ar2,phi_ar2);

  for (ii in 1:max_lag) {
    result <- VAR(y_ar1, p = ii, type = "none")
    residual <- resid(result)
    sigma_hat_ar1<- t(residual) %*% residual / sample_size
    ar1_delta[jj,ii] = (-2)*logLik(result)
    ar1_aic_select[ii,jj] <- AIC(sigma_hat_ar1,ii,2,sample_size)
    ar1_aicc_select[ii,jj] <- AICc(sigma_hat_ar1,ii,2,sample_size)
    ar1_aicbd_select[ii,jj] <- AICc.BD(sigma_hat_ar1,ii,2,sample_size)
    ar1_sic_select[ii,jj] <- SIC(sigma_hat_ar1,ii,2,sample_size)
  }
  ar1_model_select[jj,1] <- which.min(ar1_aic_select[,jj])
  ar1_model_select[jj,2] <- which.min(ar1_aicc_select[,jj])
  ar1_model_select[jj,3] <- which.min(ar1_aicbd_select[,jj])
  ar1_model_select[jj,4] <- which.min(ar1_sic_select[,jj])

  for (ii in 1:max_lag) {
    result <- VAR(y_ar2, p = ii, type = "none")
    residual <- resid(result)
    sigma_hat_ar2<- t(residual) %*% residual / sample_size
    ar2_delta[jj,ii] = (-2)*logLik(result)
    ar2_aic_select[ii,jj] <- AIC(sigma_hat_ar2,ii,2,sample_size)
    ar2_aicc_select[ii,jj] <- AICc(sigma_hat_ar2,ii,2,sample_size)
    ar2_aicbd_select[ii,jj] <- AICc.BD(sigma_hat_ar2,ii,2,sample_size)
    ar2_sic_select[ii,jj] <- SIC(sigma_hat_ar2,ii,2,sample_size)
  }
  ar2_model_select[jj,1] <- which.min(ar2_aic_select[,jj])
  ar2_model_select[jj,2] <- which.min(ar2_aicc_select[,jj])
  ar2_model_select[jj,3] <- which.min(ar2_aicbd_select[,jj])
  ar2_model_select[jj,4] <- which.min(ar2_sic_select[,jj])
}

# Product tables and figures
figure1_aic <- rowSums(ar1_aic_select)/nsim
figure1_aicc <- rowSums(ar1_aicc_select)/nsim
figure1_aicbd <- rowSums(ar1_aicbd_select)/nsim
figure1_sic <- rowSums(ar1_sic_select)/nsim
figure1_delta <- colSums(ar1_delta)/nsim

figure2_aic <- rowSums(ar2_aic_select)/nsim
figure2_aicc <- rowSums(ar2_aicc_select)/nsim
figure2_aicbd <- rowSums(ar2_aicbd_select)/nsim
figure2_sic <- rowSums(ar2_sic_select)/nsim
figure2_delta <- colSums(ar2_delta)/nsim

table1 <- data.frame(ar1_model_select)
table1_counts <- matrix(NaN,4,6)
for (ii in 1:max_lag) {
  table1_counts[1,ii] <- sum(with(table1, X1==ii))
  table1_counts[2,ii] <- sum(with(table1, X2==ii))
  table1_counts[3,ii] <- sum(with(table1, X3==ii))
  table1_counts[4,ii] <- sum(with(table1, X4==ii))
}

table2 <- data.frame(ar2_model_select)
table2_counts <- matrix(NaN,4,6)
for (ii in 1:max_lag) {
  table2_counts[1,ii] <- sum(with(table2, X1==ii))
  table2_counts[2,ii] <- sum(with(table2, X2==ii))
  table2_counts[3,ii] <- sum(with(table2, X3==ii))
  table2_counts[4,ii] <- sum(with(table2, X4==ii))
}
