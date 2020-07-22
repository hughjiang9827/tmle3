devtools::load_all("Moss")
devtools::load_all("tmle3")
library(sl3)
library(ggplot2)

reshape_long_data = function(long_data, t_max) {
      n <- length(long_data) / t_max
      # TODO: assume long_data is a list
      rs <- list()
      for (i in 1:t_max) {
        current <- long_data[seq(1 + (i - 1) * n, i * n)]
        rs <- c(rs, list(current))
      }
      rs <- do.call(cbind, rs)
      return(rs)
}

vet_data <- read.csv("https://raw.githubusercontent.com/tlverse/deming2019-workshop/master/data/veteran.csv")
vet_data$trt <- vet_data$trt - 1
# TODO: check
vet_data$time <- ceiling(vet_data$time / 20) # make fewer times for testing

T_tilde <- vet_data$time
Delta <- vet_data$status
A <- vet_data$trt
W <- vet_data[, c(3, 6:9)]
t_max <- max(T_tilde)

# TODO: lrnr
sl_lib_decent <- c("SL.xgboost")
sl_fit <- initial_sl_fit(T_tilde, Delta, A, W, t_max,
                              sl_treatment = sl_lib_decent, 
                              sl_censoring = sl_lib_decent, 
                              sl_failure = sl_lib_decent)

k_grid <- 1:t_max
sl_fit$density_failure_1$hazard_to_survival()
sl_fit$density_failure_0$hazard_to_survival()
sl_fit$density_failure_1$t <- k_grid
sl_fit$density_failure_0$t <- k_grid

test_that("sl_1 results should not be NA", {
  expect_true(all(!sapply(sl_fit$density_failure_1$survival, is.na)))
})
test_that("sl_0 results should not be NA", {
  expect_true(all(!sapply(sl_fit$density_failure_0$survival, is.na)))
})

################################################################################
# tlverse
tmax <- max(vet_data$time)
all_times <- lapply(seq_len(tmax), function(t_current){
  vet_data_time <- copy(vet_data)
  # TODO: check
  vet_data_time$N <- ifelse(t_current == vet_data_time$time & vet_data_time$status == 1, 1, 0)
  vet_data_time$A_c <- ifelse(t_current == vet_data_time$time & vet_data_time$status == 0, 1, 0)
  vet_data_time$t <- t_current

  return(vet_data_time)
})
vet_data_long <- rbindlist(all_times)

node_list <- list(W = c("celltype", "karno", "diagtime", "age", "prior"), A = "trt", T_tilde = "time", Delta = "status", 
  t = "t", N = "N", A_c = "A_c")

# TODO: lrnr
lrnr_xgb <- make_learner(Lrnr_xgboost)
learner_list <- list(A = lrnr_xgb, N = lrnr_xgb, A_c = lrnr_xgb)
# lrnr_glm <- make_learner(Lrnr_glm)
# learner_list <- list(A = lrnr_glm, N = lrnr_glm, A_c = lrnr_glm)
# lrnr_mean <- make_learner(Lrnr_mean)
# learner_list <- list(A = lrnr_mean, N = lrnr_mean, A_c = lrnr_mean)

# TODO: check
var_types <- list(T_tilde = Variable_Type$new("continuous"))
survival_spec <- tmle_survival(treatment_level = 1, control_level = 0, variable_types = var_types)
survival_task <- survival_spec$make_tmle_task(vet_data_long, node_list)

likelihood <- survival_spec$make_initial_likelihood(survival_task, learner_list)

# Update Process
initial_likelihood <- likelihood
# TODO: check
up <- tmle3_Update_survival$new(maxit = 1e1, clipping = 1e-2)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater = up)
# targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)
tmle_task <- survival_task
tmle_params <- survival_spec$make_params(survival_task, targeted_likelihood)

ps <- tmle_params[[1]]
cf_task <- ps$cf_likelihood$enumerate_cf_tasks(tmle_task)[[1]]
pN1 <- ps$observed_likelihood$get_likelihoods(cf_task, "N")
pS_N1 <- ps$hazards_to_survival(pN1, tmax)
r_pN1 <- reshape_long_data(pN1, tmax)
r_pS_N1 <- reshape_long_data(pS_N1, tmax)
sum(sl_fit$density_failure_1$survival - r_pS_N1)

psi0_moss <- colMeans(sl_fit$density_failure_1$survival)
psi0_tl <- ps$get_psi(pS_N1, tmax)
dt <- data.table(psi0_moss,psi0_tl)
dt[,t:=.I]
long <- melt(dt,id="t")
ggplot(long,aes(x=t, y=value, color=variable))+geom_line()+theme_bw()

################################################################################
# moss hazard submodel
moss_hazard_l2 <- MOSS_hazard$new(
  A = A,
  T_tilde = T_tilde,
  Delta = Delta,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  A_intervene = 1,
  k_grid = k_grid
)
psi1_moss <- moss_hazard_l2$iterate_onestep(
  method = "l2", epsilon = 1e-2, max_num_interation = 1e1, verbose = FALSE
)

# tlverse update process
tmle_fit_manual <- fit_tmle3(
    tmle_task, targeted_likelihood, tmle_params,
    targeted_likelihood$updater
)
rs <- tmle_fit_manual$estimates[[1]]
psi1_tl <- rs$psi
sum(psi1_moss - psi1_tl)

dt <- data.table(psi1_moss,psi1_tl)
dt[,t:=.I]
long <- melt(dt,id="t")
ggplot(long,aes(x=t, y=value, color=variable))+geom_line()+theme_bw()
