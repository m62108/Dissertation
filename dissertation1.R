library(dplyr)
library(nnet)
library(mice)
library(MASS)
library(caret)
library(tidyr)
library(reshape2)
library(VGAM)
library(randomForest)
library(pROC)
library(e1071)

data <- read.csv("mimiccohort.csv")
dat.mimic <- data %>%
  mutate(AKI = case_when(
    creatinine_max_day2_3 >= (3 * creatinine_min) |
      creatinine_max_day2_3 > 4 |
      (age < 18 & baseline_egfr < 35) ~ "Stage 3",
    creatinine_max_day2_3 >= (2 * creatinine_min) ~ "Stage 2",
    creatinine_max_day2_3 >= (1.5 * creatinine_min) |
      creatinine_max_day2_3 - creatinine_min >= 0.3 ~ "Stage 1",
    is.na(creatinine_max_day2_3) | is.na(creatinine_min) ~ NA_character_,
    TRUE ~ "No AKI")) %>%
  filter(!is.na(AKI))

set.seed(123)

#data preprocessing
dat.mimic <- dat.mimic %>% dplyr::select(-deathtime)
dat.mimic$age <- replace(dat.mimic$age, dat.mimic$age > 89, 90)

#mice
dat.mimic$gender<- as.factor(dat.mimic$gender)
dat.mimic$ethnicity_grouped <- as.factor(dat.mimic$ethnicity_grouped)
dry_run_imp <- mice(dat.mimic, m = 5, maxit = 0)
dry_run_imp$predictorMatrix[, c("subject_id", "hadm_id", "icustay_id", "hadm_id", "hospstay_seq")] <- 0
imputed_data <- mice(dat.mimic, m = 10, method = "pmm", dry_run_imp$predictorMatrix, maxit = 5)
dat.mimic_com <- complete(imputed_data, 1)
dat.mimic_com <- dat.mimic_com %>% dplyr::select(-subject_id, -hadm_id, -icustay_id, -intime, -outtime, -icu_los_days, -icustay_seq, -dob, 
                                                 -admittime, -dischtime, -hospital_expire_flag, -hospstay_seq, -admission_type)
dat.mimic_com$AKI_num <- as.numeric(factor(dat.mimic_com$AKI, levels = c("No AKI", "Stage 1", "Stage 2", "Stage 3")))

#divide training set
m <- 5
n <- nrow(dat.mimic_com)
folds <- sample(rep(1:m, length.out = n))
train_data_list <- list()
test_data_list <- list()
for (i in 1:m) {
  test_indices <- which(folds == i, arr.ind = TRUE)
  test_data <- dat.mimic_com[test_indices, ]
  train_data <- dat.mimic_com[-test_indices, ]
  cat("Fold", i, "Training Set Rows:", nrow(train_data), "\n")
  cat("Fold", i, "Test Set Rows:", nrow(test_data), "\n")
  train_data_list[[i]] <- train_data
  test_data_list[[i]] <- test_data
}

#select variables
summary(glm(AKI_num ~ age, data = dat.mimic_com))
initial_model <- multinom(AKI_num ~ age + factor(gender) + factor(ethnicity_grouped) + bicarbonate_min + bicarbonate_max + creatinine_min + creatinine_max + chloride_min + chloride_max + hemoglobin_min +
                            hemoglobin_max + platelet_min + platelet_max + potassium_min + potassium_max + ptt_min + ptt_max + inr_min + inr_max + pt_min + pt_max + bun_min + bun_max + wbc_min + wbc_max + baseline_egfr +
                            heartrate_min + heartrate_max + sysbp_min + sysbp_max + diasbp_min + diasbp_max + meanbp_min + meanbp_max + resprate_min + resprate_max + tempc_min + tempc_max + spo2_min + spo2_max + glucose_min + glucose_max, data = dat.mimic_com)
mlr_model <- stepAIC(initial_model, direction = "backward")
summary(mlr_model)

#calculate minimum sample size
dat.mimic1 <- dat.mimic
dat.mimic1$AKI_num <- as.numeric(factor(dat.mimic1$AKI, levels = c("No AKI", "Stage 1", "Stage 2", "Stage 3")))
p1 <- sum(dat.mimic1$AKI_num == 1)/nrow(dat.mimic1)
p2 <- sum(dat.mimic1$AKI_num == 2)/nrow(dat.mimic1)
p3 <- sum(dat.mimic1$AKI_num == 3)/nrow(dat.mimic1)
p4 <- sum(dat.mimic1$AKI_num == 4)/nrow(dat.mimic1)

p2_1 <- (sum(dat.mimic1$AKI_num == 1)+sum(dat.mimic1$AKI_num == 2))/nrow(dat.mimic1)
p3_1 <- (sum(dat.mimic1$AKI_num == 1)+sum(dat.mimic1$AKI_num == 3))/nrow(dat.mimic1)
p4_1 <- (sum(dat.mimic1$AKI_num == 1)+sum(dat.mimic1$AKI_num == 4))/nrow(dat.mimic1)
p3_2 <- (sum(dat.mimic1$AKI_num == 3)+sum(dat.mimic1$AKI_num == 2))/nrow(dat.mimic1)
p4_2 <- (sum(dat.mimic1$AKI_num == 4)+sum(dat.mimic1$AKI_num == 2))/nrow(dat.mimic1)
p4_3 <- (sum(dat.mimic1$AKI_num == 4)+sum(dat.mimic1$AKI_num == 3))/nrow(dat.mimic1)

f2_1 <- sum(dat.mimic1$AKI_num == 2)/(sum(dat.mimic1$AKI_num == 1)+sum(dat.mimic1$AKI_num == 2))
f3_1 <- sum(dat.mimic1$AKI_num == 3)/(sum(dat.mimic1$AKI_num == 1)+sum(dat.mimic1$AKI_num == 3))
f4_1 <- sum(dat.mimic1$AKI_num == 4)/(sum(dat.mimic1$AKI_num == 1)+sum(dat.mimic1$AKI_num == 4))
f3_2 <- sum(dat.mimic1$AKI_num == 3)/(sum(dat.mimic1$AKI_num == 2)+sum(dat.mimic1$AKI_num == 3))
f4_2 <- sum(dat.mimic1$AKI_num == 4)/(sum(dat.mimic1$AKI_num == 2)+sum(dat.mimic1$AKI_num == 4))
f4_3 <- sum(dat.mimic1$AKI_num == 4)/(sum(dat.mimic1$AKI_num == 3)+sum(dat.mimic1$AKI_num == 4))

maxR <- 1- (p2_1^p2_1*p3_1^p3_1*p4_1^p4_1*p3_2^p3_2*p4_2^p4_2*p4_3^p4_3)^2
R <- 0.15* maxR

R2_1 <- 0.15 * (1- (f2_1^f2_1* (1-f2_1)^(1-f2_1))^2)
R3_1 <- 0.15 * (1- (f3_1^f3_1* (1-f3_1)^(1-f3_1))^2)
R4_1 <- 0.15 * (1- (f4_1^f2_1* (1-f4_1)^(1-f4_1))^2)
R3_2 <- 0.15 * (1- (f3_2^f3_2* (1-f3_2)^(1-f3_2))^2)
R4_2 <- 0.15 * (1- (f4_2^f4_2* (1-f4_2)^(1-f4_2))^2)
R4_3 <- 0.15 * (1- (f4_3^f4_3* (1-f4_3)^(1-f4_3))^2)

#criteria 1
n2_1 <- (23/((0.9-1)*log(1-R2_1/0.9)))/p2_1
n3_1 <- (23/((0.9-1)*log(1-R3_1/0.9)))/p3_1
n4_1 <- (23/((0.9-1)*log(1-R4_1/0.9)))/p4_1
n3_2 <- (23/((0.9-1)*log(1-R3_2/0.9)))/p3_2
n4_2 <- (23/((0.8-1)*log(1-R4_2/0.9)))/p4_2
n4_3 <- (23/((0.6-1)*log(1-R4_3/0.9)))/p4_3

#criteria 2
sample2 <- 3*23/(((R/(R+0.05*maxR))-1)*log(1-R-0.05*maxR))

#criteria 3
x1_2 <- qchisq(p = 0.05/4, df = 1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
n1 <- (x1_2*p1*(1-p1))/0.05^2
n2 <- (x1_2*p2*(1-p2))/0.05^2
n3 <- (x1_2*p3*(1-p3))/0.05^2
n4 <- (x1_2*p4*(1-p4))/0.05^2

#mlr
mlr_models <- list()
mlr_pred <- list()
mlr_log <- list()
for (i in 1:5) {
  train_new <- train_data_list[[i]]
  test_new <- test_data_list[[i]]
  mlr_model <- multinom(AKI_num ~ age + gender + bicarbonate_min + 
                             creatinine_min + creatinine_max + hemoglobin_min + hemoglobin_max + 
                             platelet_max + potassium_min + potassium_max + ptt_min + 
                             ptt_max + inr_min + bun_min + bun_max + baseline_egfr + sysbp_min + 
                             sysbp_max + diasbp_max + resprate_min + resprate_max + tempc_min + 
                         glucose_min, data = train_new)
  mlr_models[[i]] <- mlr_model
  mlr_predictions <- predict(mlr_model, test_new, type = "prob")
  mlr_ref_predictions <- mlr_predictions[, 1]
  log <- log(mlr_predictions / mlr_ref_predictions)
  log <- log[, -1]
  mlr_pred[[i]] <- mlr_predictions
  mlr_log[[i]] <- log
}

#rf
rf_models <- list()
rf_pred <- list()
rf_log <- list()
for (i in 1:5) {
  train_new <- train_data_list[[i]]
  test_new <- test_data_list[[i]]
  train_new$AKI_num <- as.factor(train_new$AKI_num)
  test_new$AKI_num <- as.factor(test_new$AKI_num)
  rf_model <- randomForest(AKI_num ~ age + gender + bicarbonate_min + 
                             creatinine_min + creatinine_max + hemoglobin_min + hemoglobin_max + 
                             platelet_max + potassium_min + potassium_max + ptt_min + 
                             ptt_max + inr_min + bun_min + bun_max + baseline_egfr + sysbp_min + 
                             sysbp_max + diasbp_max + resprate_min + resprate_max + tempc_min + 
                             glucose_min, data = train_new)
  rf_models[[i]] <- rf_model
  rf_predictions <- predict(rf_model, test_new, type = "prob")
  rf_ref_predictions <- rf_predictions[, 1]
  log <- log(rf_predictions / rf_ref_predictions)
  log <- log[, -1]
  rf_pred[[i]] <- rf_predictions
  rf_log[[i]] <- log
}

#ovo dataset
train_pc_list <- list()
test_pc_list <- list()
for (i in 1:5) {
  train_new <- train_data_list[[i]]
  test_new <- test_data_list[[i]]
  train_new_pc <- train_new %>%
    mutate(Y1.2 = case_when(AKI_num == 1 ~ 1, AKI_num == 2 ~ 0, TRUE ~ NA_real_),
           Y1.3 = case_when(AKI_num == 1 ~ 1, AKI_num == 3 ~ 0, TRUE ~ NA_real_),
           Y1.4 = case_when(AKI_num == 1 ~ 1, AKI_num == 4 ~ 0, TRUE ~ NA_real_),
           Y2.3 = case_when(AKI_num == 2 ~ 1, AKI_num == 3 ~ 0, TRUE ~ NA_real_),
           Y2.4 = case_when(AKI_num == 2 ~ 1, AKI_num == 4 ~ 0, TRUE ~ NA_real_),
           Y3.4 = case_when(AKI_num == 3 ~ 1, AKI_num == 4 ~ 0, TRUE ~ NA_real_)) %>%
    mutate(across(starts_with("Y"), ~ factor(.))) 
  test_new_pc <- test_new %>%
    mutate(Y1.2 = case_when(AKI_num == 1 ~ 1, AKI_num == 2 ~ 0, TRUE ~ NA_real_),
           Y1.3 = case_when(AKI_num == 1 ~ 1, AKI_num == 3 ~ 0, TRUE ~ NA_real_),
           Y1.4 = case_when(AKI_num == 1 ~ 1, AKI_num == 4 ~ 0, TRUE ~ NA_real_),
           Y2.3 = case_when(AKI_num == 2 ~ 1, AKI_num == 3 ~ 0, TRUE ~ NA_real_),
           Y2.4 = case_when(AKI_num == 2 ~ 1, AKI_num == 4 ~ 0, TRUE ~ NA_real_),
           Y3.4 = case_when(AKI_num == 3 ~ 1, AKI_num == 4 ~ 0, TRUE ~ NA_real_)) %>%
    mutate(across(starts_with("Y"), ~ factor(.))) 
  train_pc_list[[i]] <- train_new_pc
  test_pc_list[[i]] <- test_new_pc
}

#lr-ovo
fit.OvO.list <- list()
pairwise.p.list <- list()
lr_pc_predictions.list <- list()
lr_pc_log.list <- list()
for (i in 1:5) {
  train_new_pc <- train_pc_list[[i]]
  test_new_pc <- test_pc_list[[i]]
  test_features <- test_new_pc %>%
    dplyr::select(age, gender, bicarbonate_min, creatinine_min, 
                  creatinine_max, hemoglobin_min, hemoglobin_max, platelet_max, 
                  potassium_min, potassium_max, ptt_min, ptt_max, inr_min, bun_min, 
                  bun_max, baseline_egfr, sysbp_min, sysbp_max, diasbp_max, 
                  resprate_min, resprate_max, tempc_min, glucose_min)
  fit.OvO.list[[i]] <- vector("list", 6)
  fit.OvO.list[[i]][[1]] <- glm(Y1.2 ~ age + gender + bicarbonate_min + 
                                  creatinine_min + creatinine_max + hemoglobin_min + hemoglobin_max + 
                                  platelet_max + potassium_min + potassium_max + ptt_min + 
                                  ptt_max + inr_min + bun_min + bun_max + baseline_egfr + sysbp_min + 
                                  sysbp_max + diasbp_max + resprate_min + resprate_max + tempc_min + 
                                  glucose_min, family = binomial(link = "logit"), 
                                data = subset(train_new_pc, AKI_num %in% c(1,2)) %>%
                                  dplyr::select(-c(Y1.3, Y1.4, Y2.3, Y2.4, Y3.4)))
  fit.OvO.list[[i]][[2]] <- glm(Y1.3 ~ age + gender + bicarbonate_min + 
                                  creatinine_min + creatinine_max + hemoglobin_min + hemoglobin_max + 
                                  platelet_max + potassium_min + potassium_max + ptt_min + 
                                  ptt_max + inr_min + bun_min + bun_max + baseline_egfr + sysbp_min + 
                                  sysbp_max + diasbp_max + resprate_min + resprate_max + tempc_min + 
                                  glucose_min, family = binomial(link = "logit"), 
                                data = subset(train_new_pc, AKI_num %in% c(1,3)) %>%
                                  dplyr::select(-c(Y1.2, Y1.4, Y2.3, Y2.4, Y3.4)))
  fit.OvO.list[[i]][[3]] <- glm(Y1.4 ~ age + gender + bicarbonate_min + 
                                  creatinine_min + creatinine_max + hemoglobin_min + hemoglobin_max + 
                                  platelet_max + potassium_min + potassium_max + ptt_min + 
                                  ptt_max + inr_min + bun_min + bun_max + baseline_egfr + sysbp_min + 
                                  sysbp_max + diasbp_max + resprate_min + resprate_max + tempc_min + 
                                  glucose_min, family = binomial(link = "logit"), 
                                data = subset(train_new_pc, AKI_num %in% c(1,4)) %>%
                                  dplyr::select(-c(Y1.2, Y1.3, Y2.3, Y2.4, Y3.4)))
  fit.OvO.list[[i]][[4]] <- glm(Y2.3 ~ age + gender + bicarbonate_min + 
                                  creatinine_min + creatinine_max + hemoglobin_min + hemoglobin_max + 
                                  platelet_max + potassium_min + potassium_max + ptt_min + 
                                  ptt_max + inr_min + bun_min + bun_max + baseline_egfr + sysbp_min + 
                                  sysbp_max + diasbp_max + resprate_min + resprate_max + tempc_min + 
                                  glucose_min, family = binomial(link = "logit"), 
                                data = subset(train_new_pc, AKI_num %in% c(2,3)) %>%
                                  dplyr::select(-c(Y1.2, Y1.3, Y1.4, Y2.4, Y3.4)))
  fit.OvO.list[[i]][[5]] <- glm(Y2.4 ~ age + gender + bicarbonate_min + 
                                  creatinine_min + creatinine_max + hemoglobin_min + hemoglobin_max + 
                                  platelet_max + potassium_min + potassium_max + ptt_min + 
                                  ptt_max + inr_min + bun_min + bun_max + baseline_egfr + sysbp_min + 
                                  sysbp_max + diasbp_max + resprate_min + resprate_max + tempc_min + 
                                  glucose_min, family = binomial(link = "logit"), 
                                data = subset(train_new_pc, AKI_num %in% c(2,4)) %>%
                                  dplyr::select(-c(Y1.2, Y1.3, Y1.4, Y2.3, Y3.4)))
  fit.OvO.list[[i]][[6]] <- glm(Y3.4 ~ age + gender + bicarbonate_min + 
                                  creatinine_min + creatinine_max + hemoglobin_min + hemoglobin_max + 
                                  platelet_max + potassium_min + potassium_max + ptt_min + 
                                  ptt_max + inr_min + bun_min + bun_max + baseline_egfr + sysbp_min + 
                                  sysbp_max + diasbp_max + resprate_min + resprate_max + tempc_min + 
                                  glucose_min, family = binomial(link = "logit"), 
                                data = subset(train_new_pc, AKI_num %in% c(3,4)) %>%
                                  dplyr::select(-c(Y1.2, Y1.3, Y1.4, Y2.3, Y2.4)))
  pairwise.p <- lapply(fit.OvO.list[[i]], function(x){
    predict(x, newdata = test_features, type = "response")
  })
  pairwise.p <- do.call(cbind, pairwise.p, 1)
  pairwise.p.list[[i]] <- pairwise.p
}
pairwise.couple <- function(probs.in, K){
  probs <- array(probs.in, dim = c(1, K*(K-1)/2))
  Q <- matrix(0,K,K)
  Q[lower.tri(Q)] <- 1 - probs
  Qt <- t(Q)
  Q[upper.tri(Q)] <- 1 - Qt[upper.tri(Qt)]
  diag(Q) <- rowSums(Q)
  Q <- Q / (K-1)
  p <- rbeta(K,1,1)
  p <- p/sum(p)
  for(i in 1:1000) p <- Q %*% p
  return(t(p))
}
for (i in 1:5) {
  lr_pc_predictions <- t(apply(pairwise.p.list[[i]], 1, function(x) pairwise.couple(x, k)))
  colnames(lr_pc_predictions) <- as.character(1:ncol(lr_pc_predictions))
  lr_pc_predictions.list[[i]] <- lr_pc_predictions
  lr_pc_ref_predictions <- lr_pc_predictions[, 1]
  lr_pc_log <- log(lr_pc_predictions / lr_pc_ref_predictions)
  lr_pc_log <- lr_pc_log[, -1]
  lr_pc_log.list[[i]] <- lr_pc_log
}

#SVM-ovo
fit.OvO.list1 <- list()
svm_pairwise.p.list <- list()
svm_pc_predictions.list <- list()
svm_pc_log.list <- list()
for (i in 1:5) {
    train_new_pc <- train_pc_list[[i]]
    test_new_pc <- test_pc_list[[i]]
    test_features <- test_new_pc %>%
      dplyr::select(age, gender, bicarbonate_min, creatinine_min, creatinine_max, 
                    hemoglobin_min, hemoglobin_max, platelet_max, potassium_min, 
                    potassium_max, ptt_min, ptt_max, inr_min, bun_min, bun_max, 
                    baseline_egfr, sysbp_min, sysbp_max, diasbp_max, resprate_min, 
                    resprate_max, tempc_min, glucose_min)
    fit.OvO.list1[[i]] <- vector("list", 6)
    fit.OvO.list1[[i]][[1]] <- svm(Y1.2 ~ age + gender + bicarbonate_min + 
                                     creatinine_min + creatinine_max + hemoglobin_min + hemoglobin_max + 
                                     platelet_max + potassium_min + potassium_max + ptt_min + 
                                     ptt_max + inr_min + bun_min + bun_max + baseline_egfr + sysbp_min + 
                                     sysbp_max + diasbp_max + resprate_min + resprate_max + tempc_min + 
                                     glucose_min, 
                                   data = subset(train_new_pc, AKI_num %in% c(1,2)) %>%
                                     dplyr::select(-c(Y1.3, Y1.4, Y2.3, Y2.4, Y3.4)), probability = TRUE)
    fit.OvO.list1[[i]][[2]] <- svm(Y1.3 ~ age + gender + bicarbonate_min + 
                                     creatinine_min + creatinine_max + hemoglobin_min + hemoglobin_max + 
                                     platelet_max + potassium_min + potassium_max + ptt_min + 
                                     ptt_max + inr_min + bun_min + bun_max + baseline_egfr + sysbp_min + 
                                     sysbp_max + diasbp_max + resprate_min + resprate_max + tempc_min + 
                                     glucose_min, 
                                   data = subset(train_new_pc, AKI_num %in% c(1,3)) %>%
                                     dplyr::select(-c(Y1.2, Y1.4, Y2.3, Y2.4, Y3.4)), probability = TRUE)
    fit.OvO.list1[[i]][[3]] <- svm(Y1.4 ~ age + gender + bicarbonate_min + 
                                     creatinine_min + creatinine_max + hemoglobin_min + hemoglobin_max + 
                                     platelet_max + potassium_min + potassium_max + ptt_min + 
                                     ptt_max + inr_min + bun_min + bun_max + baseline_egfr + sysbp_min + 
                                     sysbp_max + diasbp_max + resprate_min + resprate_max + tempc_min + 
                                     glucose_min, 
                                   data = subset(train_new_pc, AKI_num %in% c(1,4)) %>%
                                     dplyr::select(-c(Y1.2, Y1.3, Y2.3, Y2.4, Y3.4)), probability = TRUE)
    fit.OvO.list1[[i]][[4]] <- svm(Y2.3 ~ age + gender + bicarbonate_min + 
                                     creatinine_min + creatinine_max + hemoglobin_min + hemoglobin_max + 
                                     platelet_max + potassium_min + potassium_max + ptt_min + 
                                     ptt_max + inr_min + bun_min + bun_max + baseline_egfr + sysbp_min + 
                                     sysbp_max + diasbp_max + resprate_min + resprate_max + tempc_min + 
                                     glucose_min, 
                                   data = subset(train_new_pc, AKI_num %in% c(2,3)) %>%
                                     dplyr::select(-c(Y1.2, Y1.3, Y1.4, Y2.4, Y3.4)), probability = TRUE)
    fit.OvO.list1[[i]][[5]] <- svm(Y2.4 ~ age + gender + bicarbonate_min + 
                                     creatinine_min + creatinine_max + hemoglobin_min + hemoglobin_max + 
                                     platelet_max + potassium_min + potassium_max + ptt_min + 
                                     ptt_max + inr_min + bun_min + bun_max + baseline_egfr + sysbp_min + 
                                     sysbp_max + diasbp_max + resprate_min + resprate_max + tempc_min + 
                                     glucose_min, 
                                   data = subset(train_new_pc, AKI_num %in% c(2,4)) %>%
                                     dplyr::select(-c(Y1.2, Y1.3, Y1.4, Y2.3, Y3.4)), probability = TRUE)
    fit.OvO.list1[[i]][[6]] <- svm(Y3.4 ~ age + gender + bicarbonate_min + 
                                     creatinine_min + creatinine_max + hemoglobin_min + hemoglobin_max + 
                                     platelet_max + potassium_min + potassium_max + ptt_min + 
                                     ptt_max + inr_min + bun_min + bun_max + baseline_egfr + sysbp_min + 
                                     sysbp_max + diasbp_max + resprate_min + resprate_max + tempc_min + 
                                     glucose_min, 
                                   data = subset(train_new_pc, AKI_num %in% c(3,4)) %>%
                                     dplyr::select(-c(Y1.2, Y1.3, Y1.4, Y2.3, Y2.4)), probability = TRUE)
    svm_pairwise.p <- lapply(fit.OvO.list1[[i]], function(x){
      pred_prob <- predict(x, newdata = test_features, probability = TRUE)
      attr(pred_prob, "probabilities")[, 1]
    })
    svm_pairwise.p <- do.call(cbind, svm_pairwise.p, 1)
    svm_pc_predictions <- t(apply(svm_pairwise.p, 1, function(x) pairwise.couple(x, k)))
    colnames(svm_pc_predictions) <- as.character(1:ncol(svm_pc_predictions))
    svm_pc_ref_predictions <- svm_pc_predictions[, 1]
    svm_pc_log <- log(svm_pc_predictions / svm_pc_ref_predictions)
    svm_pc_log <- svm_pc_log[, -1]
    
    svm_pairwise.p.list[[i]] <- svm_pairwise.p
    svm_pc_predictions.list[[i]] <- svm_pc_predictions
    svm_pc_log.list[[i]] <- svm_pc_log
}

#sequential dataset
train_seq_list <- list()
test_seq_list <- list()
for (i in 1:5) {
  train_new <- train_data_list[[i]]
  test_new <- test_data_list[[i]]
  train_new <- train_new %>%
    mutate(AKI_num = as.numeric(as.character(AKI_num)))
  test_new <- test_new %>%
    mutate(AKI_num = as.numeric(as.character(AKI_num)))
  
  train_new_seq<- train_new %>% mutate(Y1 = case_when(AKI_num == 1 ~ 1, AKI_num > 1 ~ 0, TRUE ~ NA_real_),
                                       Y2 = case_when(AKI_num == 2 ~ 1, AKI_num > 2 ~ 0, TRUE ~ NA_real_),
                                       Y3 = case_when(AKI_num == 3 ~ 1, AKI_num > 3 ~ 0, TRUE ~ NA_real_)) %>%
    mutate(across(starts_with("Y"), ~ factor(.))) 
  test_new_seq<- test_new %>% mutate(Y1 = case_when(AKI_num == 1 ~ 1, AKI_num > 1 ~ 0, TRUE ~ NA_real_),
                                     Y2 = case_when(AKI_num == 2 ~ 1, AKI_num > 2 ~ 0, TRUE ~ NA_real_),
                                     Y3 = case_when(AKI_num == 3 ~ 1, AKI_num > 3 ~ 0, TRUE ~ NA_real_)) %>%
    mutate(across(starts_with("Y"), ~ factor(.))) 
  
  train_seq_list[[i]] <- train_new_seq
  test_seq_list[[i]] <- test_new_seq
}

#sequential LR
fit.OvO.list2 <- vector("list", 5)
lr_seq.p.list <- vector("list", 5)
lr_seq_predictions1.list <- vector("list", 5)
lr_seq_log.list <- vector("list", 5)
for (i in 1:5) {
  train_new_seq <- train_seq_list[[i]]
  test_new_seq <- test_seq_list[[i]]
  test_features <- test_new_seq %>%
    dplyr::select(age, gender, creatinine_min, creatinine_max, 
                  hemoglobin_min, platelet_min, potassium_min, potassium_max, 
                  ptt_min, ptt_max, inr_min, bun_min, bun_max, wbc_max, 
                  baseline_egfr, sysbp_min, sysbp_max, diasbp_max, meanbp_min, 
                  resprate_min, resprate_max, tempc_min, glucose_min)
  fit.OvO.list2[[i]] <- vector("list", 3)
  fit.OvO.list2[[i]][[1]] <- glm(Y1 ~ age + gender + creatinine_min + 
                                   creatinine_max + hemoglobin_min + platelet_min + potassium_min + 
                                   potassium_max + ptt_min + ptt_max + inr_min + bun_min + bun_max + 
                                   wbc_max + baseline_egfr + sysbp_min + sysbp_max + diasbp_max + 
                                   meanbp_min + resprate_min + resprate_max + tempc_min + glucose_min,
                                 family = binomial(link = "logit"), 
                                 data = train_new_seq %>% 
                                   dplyr::select(-c(Y2, Y3)))
  
  fit.OvO.list2[[i]][[2]] <- glm(Y2 ~ age + gender + creatinine_min + 
                                   creatinine_max + hemoglobin_min + platelet_min + potassium_min + 
                                   potassium_max + ptt_min + ptt_max + inr_min + bun_min + bun_max + 
                                   wbc_max + baseline_egfr + sysbp_min + sysbp_max + diasbp_max + 
                                   meanbp_min + resprate_min + resprate_max + tempc_min + glucose_min, 
                                 family = binomial(link = "logit"), 
                                 data = subset(train_new_seq, AKI_num %in% c(2,3,4)) %>% 
                                   dplyr::select(-c(Y1, Y3)))
  
  fit.OvO.list2[[i]][[3]] <- glm(Y3 ~ age + gender + creatinine_min + 
                                   creatinine_max + hemoglobin_min + platelet_min + potassium_min + 
                                   potassium_max + ptt_min + ptt_max + inr_min + bun_min + bun_max + 
                                   wbc_max + baseline_egfr + sysbp_min + sysbp_max + diasbp_max + 
                                   meanbp_min + resprate_min + resprate_max + tempc_min + glucose_min, 
                                 family = binomial(link = "logit"), 
                                 data = subset(train_new_seq, AKI_num %in% c(3,4)) %>% 
                                   dplyr::select(-c(Y1, Y2)))
  lr_seq.p <- lapply(fit.OvO.list2[[i]], function(x){
    predict(x, newdata = test_features, type = "response")
  })
  lr_seq.p <- do.call(cbind, lr_seq.p, 1)
  
  lr_seq_predictions1 <- matrix(0, nrow = nrow(test_new_seq), ncol = 4)
  lr_seq_predictions1[, 1] <- lr_seq.p[, 1]
  lr_seq_predictions1[, 2] <- (1 - lr_seq.p[, 1]) * lr_seq.p[, 2]
  lr_seq_predictions1[, 3] <- (1 - lr_seq.p[, 1]) * (1 - lr_seq.p[, 2]) * lr_seq.p[, 3]
  lr_seq_predictions1[, 4] <- (1 - lr_seq.p[, 1]) * (1 - lr_seq.p[, 2]) * (1 - lr_seq.p[, 3])
  
  colnames(lr_seq_predictions1) <- as.character(1:ncol(lr_seq_predictions1))
  lr_seq_ref_predictions1 <- lr_seq_predictions1[, 1]
  lr_seq_log <- log(lr_seq_predictions1 / lr_seq_ref_predictions1)
  lr_seq_log <- lr_seq_log[, -1]
  
  lr_seq.p.list[[i]] <- lr_seq.p
  lr_seq_predictions1.list[[i]] <- lr_seq_predictions1
  lr_seq_log.list[[i]] <- lr_seq_log
}

#sequential SVM
fit.OvO.list3 <- vector("list", 5)
svm_seq.p.list <- vector("list", 5)
svm_seq_predictions.list <- vector("list", 5)
svm_seq_log.list <- vector("list", 5)
for (i in 1:5) {
  train_new_seq <- train_seq_list[[i]]
  test_new_seq <- test_seq_list[[i]]
  test_features <- test_new_seq %>%
    dplyr::select(age, gender, creatinine_min, creatinine_max, 
                  hemoglobin_min, platelet_min, potassium_min, potassium_max, 
                  ptt_min, ptt_max, inr_min, bun_min, bun_max, wbc_max, 
                  baseline_egfr, sysbp_min, sysbp_max, diasbp_max, meanbp_min, 
                  resprate_min, resprate_max, tempc_min, glucose_min)
  fit.OvO.list3[[i]] <- vector("list", 3)
  fit.OvO.list3[[i]][[1]] <- svm(Y1 ~ age + gender + creatinine_min + 
                                   creatinine_max + hemoglobin_min + platelet_min + potassium_min + 
                                   potassium_max + ptt_min + ptt_max + inr_min + bun_min + bun_max + 
                                   wbc_max + baseline_egfr + sysbp_min + sysbp_max + diasbp_max + 
                                   meanbp_min + resprate_min + resprate_max + tempc_min + glucose_min,
                                 data = train_new_seq %>% 
                                   dplyr::select(-c(Y2, Y3)), probability = TRUE)
  fit.OvO.list3[[i]][[2]] <- svm(Y2 ~ age + gender + creatinine_min + 
                                   creatinine_max + hemoglobin_min + platelet_min + potassium_min + 
                                   potassium_max + ptt_min + ptt_max + inr_min + bun_min + bun_max + 
                                   wbc_max + baseline_egfr + sysbp_min + sysbp_max + diasbp_max + 
                                   meanbp_min + resprate_min + resprate_max + tempc_min + glucose_min, 
                                 data = subset(train_new_seq, AKI_num %in% c(2,3,4)) %>% 
                                   dplyr::select(-c(Y1, Y3)), probability = TRUE)
  fit.OvO.list3[[i]][[3]] <- svm(Y3 ~ age + gender + creatinine_min + 
                                   creatinine_max + hemoglobin_min + platelet_min + potassium_min + 
                                   potassium_max + ptt_min + ptt_max + inr_min + bun_min + bun_max + 
                                   wbc_max + baseline_egfr + sysbp_min + sysbp_max + diasbp_max + 
                                   meanbp_min + resprate_min + resprate_max + tempc_min + glucose_min, 
                                 data = subset(train_new_seq, AKI_num %in% c(3,4)) %>% 
                                   dplyr::select(-c(Y1, Y2)), probability = TRUE)
  svm_seq.p <- lapply(fit.OvO.list3[[i]], function(model){
    pred_prob <- predict(model, newdata = test_features, probability = TRUE)
    attr(pred_prob, "probabilities")[, 1]
  })
  svm_seq.p <- do.call(cbind, svm_seq.p, 1)
  svm_seq_predictions <- matrix(0, nrow = nrow(test_new_seq), ncol = 4)
  svm_seq_predictions[, 1] <- svm_seq.p[, 1]
  svm_seq_predictions[, 2] <- (1 - svm_seq.p[, 1]) * svm_seq.p[, 2]
  svm_seq_predictions[, 3] <- (1 - svm_seq.p[, 1]) * (1 - svm_seq.p[, 2]) * svm_seq.p[, 3]
  svm_seq_predictions[, 4] <- (1 - svm_seq.p[, 1]) * (1 - svm_seq.p[, 2]) * (1 - svm_seq.p[, 3])
  
  colnames(svm_seq_predictions) <- as.character(1:ncol(svm_seq_predictions))
  svm_seq_ref_predictions <- svm_seq_predictions[, r]
  svm_seq_log <- log(svm_seq_predictions / svm_seq_ref_predictions)
  svm_seq_log <- svm_seq_log[, -1]
  
  svm_seq.p.list[[i]] <- svm_seq.p
  svm_seq_predictions.list[[i]] <- svm_seq_predictions
  svm_seq_log.list[[i]] <- svm_seq_log
}

#linear prediction matrix (Zij)
mlr_probs_list <- vector("list", 5)
mlr_lps_list <- vector("list", 5)
rf_probs_list <- vector("list", 5)
rf_lps_list <- vector("list", 5)
lr_pc_probs_list <- vector("list", 5)
lr_pc_lps_list <- vector("list", 5)
svm_pc_probs_list <- vector("list", 5)
svm_pc_lps_list <- vector("list", 5)
lr_seq_probs_list <- vector("list", 5)
lr_seq_lps_list <- vector("list", 5)
svm_seq_probs_list <- vector("list", 5)
svm_seq_lps_list <- vector("list", 5)

for (i in 1:5) {
  mlr_probs_list[[i]] <- split(mlr_pred[[i]][, order(colnames(mlr_pred[[i]]))], col(mlr_pred[[i]][, order(colnames(mlr_pred[[i]]))]))
  mlr_lps_list[[i]] <- split(mlr_log[[i]], col(mlr_log[[i]]))
  
  rf_probs_list[[i]] <- split(rf_pred[[i]][, order(colnames(rf_pred[[i]]))], col(rf_pred[[i]][, order(colnames(rf_pred[[i]]))]))
  rf_lps_list[[i]] <- split(rf_log[[i]], col(rf_log[[i]]))

  lr_pc_probs_list[[i]] <- split(lr_pc_predictions.list[[i]][, order(colnames(lr_pc_predictions.list[[i]]))], col(lr_pc_predictions.list[[i]][, order(colnames(lr_pc_predictions.list[[i]]))]))
  lr_pc_lps_list[[i]] <- split(lr_pc_log.list[[i]], col(lr_pc_log.list[[i]]))
  
  svm_pc_probs_list[[i]] <- split(svm_pc_predictions.list[[i]][, order(colnames(svm_pc_predictions.list[[i]]))], col(svm_pc_predictions.list[[i]][, order(colnames(svm_pc_predictions.list[[i]]))]))
  svm_pc_lps_list[[i]] <- split(svm_pc_log.list[[i]], col(svm_pc_log.list[[i]]))
  
  lr_seq_probs_list[[i]] <- split(lr_seq_predictions1.list[[i]][, order(colnames(lr_seq_predictions1.list[[i]]))], col(lr_seq_predictions1.list[[i]][, order(colnames(lr_seq_predictions1.list[[i]]))]))
  lr_seq_lps_list[[i]] <- split(lr_seq_log.list[[i]], col(lr_seq_log.list[[i]]))
  
  svm_seq_probs_list[[i]] <- split(svm_seq_predictions.list[[i]][, order(colnames(svm_seq_predictions.list[[i]]))], col(svm_seq_predictions.list[[i]][, order(colnames(svm_seq_predictions.list[[i]]))]))
  svm_seq_lps_list[[i]] <- split(svm_seq_log.list[[i]], col(svm_seq_log.list[[i]]))
}

for (j in 1:5) {
  for (i in seq_along(mlr_lps_list[[j]])) {
    assign(paste("mlr_lp", j, "_", i, sep = ""), mlr_lps_list[[j]][[i]])
  }
  for (i in seq_along(rf_lps_list[[j]])) {
    assign(paste("rf_lp", j, "_", i, sep = ""), rf_lps_list[[j]][[i]])
  }
  for (i in seq_along(lr_pc_lps_list[[j]])) {
    assign(paste("lr_pc_lp", j, "_", i, sep = ""), lr_pc_lps_list[[j]][[i]])
  }
  for (i in seq_along(svm_pc_lps_list[[j]])) {
    assign(paste("svm_pc_lp", j, "_", i, sep = ""), svm_pc_lps_list[[j]][[i]])
  }
  for (i in seq_along(lr_seq_lps_list[[j]])) {
    assign(paste("lr_seq_lp", j, "_", i, sep = ""), lr_seq_lps_list[[j]][[i]])
  }
  for (i in seq_along(svm_seq_lps_list[[j]])) {
    assign(paste("svm_seq_lp", j, "_", i, sep = ""), svm_seq_lps_list[[j]][[i]])
  }
}

#OE ratio
calculate_OE <- function(dat.valid.pred, prediction_matrix){
  pred_df <- as.data.frame(prediction_matrix)
  colnames(pred_df) <- paste0("p", 1:ncol(pred_df))
  dat.valid.pred <- cbind(dat.valid.pred, pred_df)
  K <- as.numeric(max(dat.valid.pred$AKI_num))
  output.object <- rep(NA, K)
  names(output.object) <- paste("p", 1:K, sep = "")
  for (i in 1:K){
    observed <- sum(dat.valid.pred$AKI_num == i)/length(dat.valid.pred$AKI_num)
    expected <- mean(dat.valid.pred[,paste("p", i, sep = "")])
    output.object[i] <- observed/expected
  }
  return(output.object)
}
mlr_ratio <- list()
rf_ratio <- list()
lr_pc_ratio <- list()
svm_pc_ratio <- list()
lr_seq_ratio <- list()
svm_seq_ratio <- list()
for (i in 1:5) {
  test_new <- test_data_list[[i]]
  mlr_predictions <- mlr_pred[[i]]
  mlr_ratio[[i]] <- calculate_OE(test_new, mlr_predictions)
  rf_predictions <- rf_pred[[i]]
  rf_ratio[[i]] <- calculate_OE(test_new, rf_predictions)
  lr_pc_predictions <- lr_pc_predictions.list[[i]]
  lr_pc_ratio[[i]] <- calculate_OE(test_new, lr_pc_predictions)
  svm_pc_predictions <- svm_pc_predictions.list[[i]]
  svm_pc_ratio[[i]] <- calculate_OE(test_new, svm_pc_predictions)
  lr_seq_predictions1 <- lr_seq_predictions1.list[[i]]
  lr_seq_ratio[[i]] <- calculate_OE(test_new, lr_seq_predictions1)
  svm_seq_predictions <- svm_seq_predictions.list[[i]]
  svm_seq_ratio[[i]] <- calculate_OE(test_new, svm_seq_predictions)
}

r <- 1
k <- 4
fitnp_list <- vector("list", 2)
rf_fitnp_list <- vector("list", 5)
lr_pc_fitnp_list <- vector("list", 5)
svm_pc_fitnp_list <- vector("list", 5)
lr_seq_fitnp_list <- vector("list", 5)
svm_seq_fitnp_list <- vector("list", 5)
for (i in 1:5) {
  test_new <- test_data_list[[i]]
  
  lp1 <- get(paste("mlr_lp", i, "_1", sep = ""))
  lp2 <- get(paste("mlr_lp", i, "_2", sep = ""))
  lp3 <- get(paste("mlr_lp", i, "_3", sep = ""))
  
  rf_lp1 <- get(paste("rf_lp", i, "_1", sep = ""))
  rf_lp2 <- get(paste("rf_lp", i, "_2", sep = ""))
  rf_lp3 <- get(paste("rf_lp", i, "_3", sep = ""))
  
  lr_pc_lp1 <- get(paste("lr_pc_lp", i, "_1", sep = ""))
  lr_pc_lp2 <- get(paste("lr_pc_lp", i, "_2", sep = ""))
  lr_pc_lp3 <- get(paste("lr_pc_lp", i, "_3", sep = ""))
  
  svm_pc_lp1 <- get(paste("svm_pc_lp", i, "_1", sep = ""))
  svm_pc_lp2 <- get(paste("svm_pc_lp", i, "_2", sep = ""))
  svm_pc_lp3 <- get(paste("svm_pc_lp", i, "_3", sep = ""))
  
  lr_seq_lp1 <- get(paste("lr_seq_lp", i, "_1", sep = ""))
  lr_seq_lp2 <- get(paste("lr_seq_lp", i, "_2", sep = ""))
  lr_seq_lp3 <- get(paste("lr_seq_lp", i, "_3", sep = ""))
  
  svm_seq_lp1 <- get(paste("svm_seq_lp", i, "_1", sep = ""))
  svm_seq_lp2 <- get(paste("svm_seq_lp", i, "_2", sep = ""))
  svm_seq_lp3 <- get(paste("svm_seq_lp", i, "_3", sep = ""))
  
  fitnp_list[[i]] <- vgam(test_new$AKI_num ~ s(lp1, df = 3) + s(lp2, df = 3) + s(lp3, df = 3), 
                          family = multinomial(refLevel = r))
  rf_fitnp_list[[i]] <- vgam(test_new$AKI_num ~ s(rf_lp1, df = 3) + s(rf_lp2, df = 3) + s(rf_lp3, df = 3), 
                          family = multinomial(refLevel = r))
  lr_pc_fitnp_list[[i]] <- vgam(test_new$AKI_num ~ s(lr_pc_lp1, df = 3) + s(lr_pc_lp2, df = 3) + s(lr_pc_lp3, df = 3), 
                          family = multinomial(refLevel = r))
  svm_pc_fitnp_list[[i]] <- vgam(test_new$AKI_num ~ s(svm_pc_lp1, df = 3) + s(svm_pc_lp2, df = 3) + s(svm_pc_lp3, df = 3), 
                          family = multinomial(refLevel = r))
  lr_seq_fitnp_list[[i]] <- vgam(test_new$AKI_num ~ s(lr_seq_lp1, df = 3) + s(lr_seq_lp2, df = 3) + s(lr_seq_lp3, df = 3), 
                          family = multinomial(refLevel = r))
  svm_seq_fitnp_list[[i]] <- vgam(test_new$AKI_num ~ s(svm_seq_lp1, df = 3) + s(svm_seq_lp2, df = 3) + s(svm_seq_lp3, df = 3), 
                          family = multinomial(refLevel = r))
}

#calibration_plot(overall)
calibration_plot <- function(probs, fitnp, k, datapoints = TRUE, smoothing = TRUE) {
  if (isTRUE(datapoints)) {
    for (i in 1:k) {
      p <- unlist(probs[[i]])
      matplot(p, fitted(fitnp)[, i], type = "p", pch = i, col = (1 + i), lwd = 1, ylab = "", xlab = "", xlim = 0:1, ylim = 0:1)
      par(new = TRUE)
    }
  }
  ref <- rbind(c(0, 0), c(1, 1))
  matplot(ref, ref, type = "l", col = 1, lwd = 2, ylab = "Observed proportions", xlab = "Predicted probabilities", xlim = 0:1, ylim = 0:1)
  if (isTRUE(smoothing)) {
    a <- 0.5
    for (i in 1:k) {
      p <- unlist(probs[[i]])
      points(smooth.spline(p, fitted(fitnp)[, i], spar = a), type = "l", col = (1 + i), lwd = 4)
    }
  }
  legende <- c()
  for (i in 1:k) {
    if (i <= 2) {
      legende <- c("AKI 0", "AKI 1")
    }
    if (i > 2) {
      legende <- c(legende, paste("AKI ", i-1, sep = ""))
    }
  }
  legend(x = 0.7, y = (0.20 + (k - 3) * 0.05), col = 2:(k + 1), lty = 1, legend = legende)
  title(main = "Non-parametric calibration plot")
  par(new = FALSE)
}

plot_data_list <- list()
par(mfrow = c(1, 1))
par(mfrow = c(2, 3))
for (i in 1:5) {
  probs <- mlr_probs_list[[i]]
  rf_probs <- rf_probs_list[[i]]
  lr_pc_probs <- lr_pc_probs_list[[i]]
  svm_pc_probs <- svm_pc_probs_list[[i]]
  lr_seq_probs <- lr_seq_probs_list[[i]]
  svm_seq_probs <- svm_seq_probs_list[[i]]
  
  fitnp <- fitnp_list[[i]]
  rf_fitnp <- rf_fitnp_list[[i]]
  lr_pc_fitnp <- lr_pc_fitnp_list[[i]]
  svm_pc_fitnp <- svm_pc_fitnp_list[[i]]
  lr_seq_fitnp <- lr_seq_fitnp_list[[i]]
  svm_seq_fitnp <- svm_seq_fitnp_list[[i]]
  
  calibration_plot(probs, fitnp, k, datapoints = TRUE, smoothing = TRUE)
  calibration_plot(rf_probs, rf_fitnp, k, datapoints = TRUE, smoothing = TRUE)
  calibration_plot(lr_pc_probs, lr_pc_fitnp, k, datapoints = TRUE, smoothing = TRUE)
  calibration_plot(svm_pc_probs, svm_pc_fitnp, k, datapoints = TRUE, smoothing = TRUE)
  calibration_plot(lr_seq_probs, lr_seq_fitnp, k, datapoints = TRUE, smoothing = TRUE)
  calibration_plot(svm_seq_probs, svm_seq_fitnp, k, datapoints = TRUE, smoothing = TRUE)
}

#category calibration plot
calibration_plot_catagory <- function(probs_list, fitnp_list, k, datapoints = TRUE, smoothing = TRUE, plot_title) {
  layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), ncol = 2, byrow = TRUE), heights = c(0.2, 0.8, 0.2, 0.8))
  par(mar = c(2, 4, 1, 1), mgp = c(2, 1, 0))
  for (i in c(0, 1)) {  
    p_combined <- unlist(lapply(probs_list, function(x) x[[i + 1]]))  
    if (i == 0) {
      plot(density(p_combined), 
           main = paste("AKI", i), 
           xlab = "", 
           ylab = "", 
           xaxt = 'n',  
           xlim = c(0.7, 1.0), 
           ylim = c(0, max(density(p_combined)$y) * 1.5))  
    } else if (i == 1) {
      plot(density(p_combined), 
           main = paste("AKI", i), 
           xlab = "", 
           ylab = "", 
           xaxt = 'n',  
           xlim = c(0, 0.3),  
           ylim = c(0, max(density(p_combined)$y) * 1.5))  
    }
    axis(1)  
  }
  
  par(mar = c(4, 4, 1, 1), mgp = c(2, 1, 0))  
  
  for (i in c(0, 1)) { 
    if (i == 0) {
      plot(NA, 
           xlim = c(0.7, 1.0),
           ylim = c(0.7, 1.0), 
           xlab = "Predicted probabilities", 
           ylab = "Observed proportions", 
           main = NULL) 
    } else if (i == 1) {
      plot(NA, 
           xlim = c(0, 0.3), 
           ylim = c(0, 0.3), 
           xlab = "Predicted probabilities", 
           ylab = "Observed proportions", 
           main = NULL) 
    }
    abline(0, 1, col = "black", lwd = 2) 
    
    colors <- rainbow(length(probs_list))
    legends <- vector()
    
    if (isTRUE(smoothing)) {
      for (j in 1:length(probs_list)) {
        p <- unlist(probs_list[[j]][[i + 1]])
        lines(smooth.spline(p, fitted(fitnp_list[[j]])[, i + 1], spar = 0.5), col = colors[j], lwd = 2)
        legends <- c(legends, paste("Fold", j)) 
      }
    }
    
    if (length(legends) > 0) {
      legend("bottomright", legend = legends, col = colors, lwd = 2)
    }
  }
  
  par(mar = c(2, 4, 1, 1), mgp = c(2, 1, 0))  
  
  for (i in c(2, 3)) {  
    p_combined <- unlist(lapply(probs_list, function(x) x[[i + 1]])) 
    plot(density(p_combined), 
         main = paste("AKI", i), 
         xlab = "", 
         ylab = "", 
         xaxt = 'n',  
         xlim = c(0, 0.1),  
         ylim = c(0, max(density(p_combined)$y) * 1.5)) 
    axis(1) 
  }
  
  par(mar = c(4, 4, 1, 1), mgp = c(2, 1, 0))  
  
  for (i in c(2, 3)) {  
    plot(NA, 
         xlim = c(0, 0.1),  
         ylim = c(0, 0.1),  
         xlab = "Predicted probabilities", 
         ylab = "Observed proportions", 
         main = NULL)  
    abline(0, 1, col = "black", lwd = 2) 
    
    colors <- rainbow(length(probs_list))
    legends <- vector()
    
    if (isTRUE(smoothing)) {
      for (j in 1:length(probs_list)) {
        p <- unlist(probs_list[[j]][[i + 1]])
        lines(smooth.spline(p, fitted(fitnp_list[[j]])[, i + 1], spar = 0.5), col = colors[j], lwd = 2)
        legends <- c(legends, paste("Fold", j))
      }
    }
    
    if (length(legends) > 0) {
      legend("bottomright", legend = legends, col = colors, lwd = 2)
    }
  }
  layout(1) 
}

calibration_plot_catagory(mlr_probs_list, fitnp_list, k, datapoints = FALSE, smoothing = TRUE, plot_title = "MLR Model")
calibration_plot_catagory(rf_probs_list, rf_fitnp_list, k, datapoints = FALSE, smoothing = TRUE, plot_title = "RF Model")
calibration_plot_catagory(lr_pc_probs_list, lr_pc_fitnp_list, k, datapoints = FALSE, smoothing = TRUE, plot_title = "LR-OvO Model")
calibration_plot_catagory(svm_pc_probs_list, svm_pc_fitnp_list, k, datapoints = FALSE, smoothing = TRUE, plot_title = "SVM-OvO Model")
calibration_plot_catagory(lr_seq_probs_list, lr_seq_fitnp_list, k, datapoints = FALSE, smoothing = TRUE, plot_title = "LR-Seq Model")
calibration_plot_catagory(svm_seq_probs_list, svm_seq_fitnp_list, k, datapoints = FALSE, smoothing = TRUE, plot_title = "SVM- Model")

#calibration intercept
calculate_calibration_intercept <- function(data, offset_matrix, r, alpha = 0.05) {
  int <- vgam(data$AKI_num ~ 1, offset = offset_matrix, family = multinomial(refLevel = r))
  coeffint <- coefficients(int)
  se <- sqrt(diag(vcov(int)))
  ci <- matrix(NA, nrow = length(coeffint), ncol = 2)
  for (i in seq_along(coeffint)) {
    ci[i, ] <- c(coeffint[i] - qnorm(1 - alpha / 2) * se[i], 
                 coeffint[i] + qnorm(1 - alpha / 2) * se[i])
  }
  estint <- c(coeffint, ci)
  names(estint) <- paste('CALINT', c(paste0('int', seq_along(coeffint)), 
                                     paste0('LLint', seq_along(coeffint)), 
                                     paste0('ULint', seq_along(coeffint))), sep = '.')
  return(estint)
}

mlr_intercept_list <- vector("list", length = 5)
rf_intercept_list <- vector("list", length = 5)
lr_pc_intercept_list <- vector("list", length = 5)
svm_pc_intercept_list <- vector("list", length = 5)
lr_seq_intercept_list <- vector("list", length = 5)
svm_seq_intercept_list <- vector("list", length = 5)

for (i in 1:5) {
  test_new <- test_data_list[[i]]
  
  mlr_intercept_list[[i]] <- calculate_calibration_intercept(test_new, mlr_log[[i]], r)
  print(mlr_intercept_list[[i]])
  
  valid_rows <- apply(rf_log[[i]], 1, function(row) all(is.finite(row) & !is.na(row)))
  rf_log_clean <- rf_log[[i]][valid_rows, , drop = FALSE]
  t_clean <- test_new[valid_rows, , drop = FALSE]
  
  rf_intercept_list[[i]] <- calculate_calibration_intercept(t_clean, rf_log_clean, r)
  print(rf_intercept_list[[i]])
  
  lr_pc_intercept_list[[i]] <- calculate_calibration_intercept(test_new, lr_pc_log.list[[i]], r)
  print(lr_pc_intercept_list[[i]])
  
  svm_pc_intercept_list[[i]] <- calculate_calibration_intercept(test_new, svm_pc_log.list[[i]], r)
  print(svm_pc_intercept_list[[i]])
  
  lr_seq_intercept_list[[i]] <- calculate_calibration_intercept(test_new, lr_seq_log.list[[i]], r)
  print(lr_seq_intercept_list[[i]])
  
  svm_seq_intercept_list[[i]] <- calculate_calibration_intercept(test_new, svm_seq_log.list[[i]], r)
  print(svm_seq_intercept_list[[i]])
}


#calibration_slopes
calculate_calibration_slopes <- function(data, lp1, lp2, lp3, k, r, alpha = 0.05) {
  i <- diag(k - 1)
  i2 <- rbind(1, 0, 0)
  i3 <- rbind(0, 1, 0)
  i4 <- rbind(0, 0, 1)
  clist <- list("(Intercept)" = i, "lp1" = i2, "lp2" = i3, "lp3" = i4)
  slopes <- vgam(data$AKI_num ~ lp1 + lp2 + lp3, family = multinomial(refLevel = r), constraints = clist)
  coeffslopes <- coefficients(slopes)[k:length(coefficients(slopes))]
  se <- sqrt(diag(vcov(slopes)))
  ci <- matrix(NA, nrow = length(coeffslopes), ncol = 2)
  for (i in seq_along(coeffslopes)) {
    ci[i, ] <- c(coeffslopes[i] - qnorm(1 - alpha / 2) * se[i + 3], 
                 coeffslopes[i] + qnorm(1 - alpha / 2) * se[i + 3])
  }
  estslopes <- c(coeffslopes, ci)
  names(estslopes) <- paste('CALSLOPES', c(paste0('lp', seq_along(coeffslopes)), 
                                           paste0('LLlp', seq_along(coeffslopes)), 
                                           paste0('ULlp', seq_along(coeffslopes))), sep = '.')
  return(estslopes)
}

mlr_estslopes_list <- vector("list", length = 5)
rf_estslopes_list <- vector("list", length = 5)
lr_pc_estslopes_list <- vector("list", length = 5)
svm_pc_estslopes_list <- vector("list", length = 5)
lr_seq_estslopes_list <- vector("list", length = 5)
svm_seq_estslopes_list <- vector("list", length = 5)

for (i in 1:5) {
  test_new <- test_data_list[[i]]
  mlr_estslopes_list[[i]] <- calculate_calibration_slopes(test_new, 
                                                          get(paste("mlr_lp", i, "_1", sep = "")), 
                                                          get(paste("mlr_lp", i, "_2", sep = "")), 
                                                          get(paste("mlr_lp", i, "_3", sep = "")), 
                                                          k, r)
  print(mlr_estslopes_list[[i]])
  
  rf_lp1 <- get(paste("rf_lp", i, "_1", sep = ""))
  rf_lp2 <- get(paste("rf_lp", i, "_2", sep = ""))
  rf_lp3 <- get(paste("rf_lp", i, "_3", sep = ""))

  valid_indices_lp1 <- which(is.finite(rf_lp1) & !is.na(rf_lp1))
  valid_indices_lp2 <- which(is.finite(rf_lp2) & !is.na(rf_lp2))
  valid_indices_lp3 <- which(is.finite(rf_lp3) & !is.na(rf_lp3))
  
  valid_indices_all <- Reduce(intersect, list(valid_indices_lp1, valid_indices_lp2, valid_indices_lp3))
  
  test_new_clean <- test_new[valid_indices_all, , drop = FALSE]
  rf_lp1_clean <- rf_lp1[valid_indices_all]
  rf_lp2_clean <- rf_lp2[valid_indices_all]
  rf_lp3_clean <- rf_lp3[valid_indices_all]
  rf_estslopes_list[[i]] <- calculate_calibration_slopes(test_new_clean, 
                                                         rf_lp1_clean, 
                                                         rf_lp2_clean, 
                                                         rf_lp3_clean, 
                                                         k, r)
  print(rf_estslopes_list[[i]])
  
  lr_pc_estslopes_list[[i]] <- calculate_calibration_slopes(test_new, 
                                                            get(paste("lr_pc_lp", i, "_1", sep = "")), 
                                                            get(paste("lr_pc_lp", i, "_2", sep = "")), 
                                                            get(paste("lr_pc_lp", i, "_3", sep = "")), 
                                                            k, r)
  print(lr_pc_estslopes_list[[i]])
  
  svm_pc_estslopes_list[[i]] <- calculate_calibration_slopes(test_new, 
                                                             get(paste("svm_pc_lp", i, "_1", sep = "")), 
                                                             get(paste("svm_pc_lp", i, "_2", sep = "")), 
                                                             get(paste("svm_pc_lp", i, "_3", sep = "")), 
                                                             k, r)
  print(svm_pc_estslopes_list[[i]])
  
  lr_seq_estslopes_list[[i]] <- calculate_calibration_slopes(test_new, 
                                                             get(paste("lr_seq_lp", i, "_1", sep = "")), 
                                                             get(paste("lr_seq_lp", i, "_2", sep = "")), 
                                                             get(paste("lr_seq_lp", i, "_3", sep = "")), 
                                                             k, r)
  print(lr_seq_estslopes_list[[i]])
  
  svm_seq_estslopes_list[[i]] <- calculate_calibration_slopes(test_new, 
                                                              get(paste("svm_seq_lp", i, "_1", sep = "")), 
                                                              get(paste("svm_seq_lp", i, "_2", sep = "")), 
                                                              get(paste("svm_seq_lp", i, "_3", sep = "")), 
                                                              k, r)
  print(svm_seq_estslopes_list[[i]])
}

#ECI
calculate_ECI <- function(dat.valid.pred, prediction_matrix) {
  dat.valid.pred$p1 <- prediction_matrix[, 1]
  dat.valid.pred$p2 <- prediction_matrix[, 2]
  dat.valid.pred$p3 <- prediction_matrix[, 3]
  dat.valid.pred$p4 <- prediction_matrix[, 4]
  dat.valid.pred$Y.fac1 <- ifelse(dat.valid.pred$AKI_num == 1, 1, 0)
  dat.valid.pred$Y.fac2 <- ifelse(dat.valid.pred$AKI_num == 2, 1, 0)
  dat.valid.pred$Y.fac3 <- ifelse(dat.valid.pred$AKI_num == 3, 1, 0)
  dat.valid.pred$Y.fac4 <- ifelse(dat.valid.pred$AKI_num == 4, 1, 0)
  loess.outcome1 <- loess(Y.fac1 ~ p1, data = dat.valid.pred, method = 'loess')
  loess.outcome2 <- loess(Y.fac2 ~ p2, data = dat.valid.pred, method = 'loess')
  loess.outcome3 <- loess(Y.fac3 ~ p3, data = dat.valid.pred, method = 'loess')
  loess.outcome4 <- loess(Y.fac4 ~ p4, data = dat.valid.pred, method = 'loess')
  dat.valid.pred$p.obs1 <- loess.outcome1$fitted
  dat.valid.pred$p.obs2 <- loess.outcome2$fitted
  dat.valid.pred$p.obs3 <- loess.outcome3$fitted
  dat.valid.pred$p.obs4 <- loess.outcome4$fitted
  ECI.numer <- sum((dat.valid.pred$p1 - dat.valid.pred$p.obs1)^2) + 
    sum((dat.valid.pred$p2 - dat.valid.pred$p.obs2)^2) +
    sum((dat.valid.pred$p3 - dat.valid.pred$p.obs3)^2) + 
    sum((dat.valid.pred$p4 - dat.valid.pred$p.obs4)^2)
  ECI.denom <- nrow(dat.valid.pred) * 4
  ECI <- ECI.numer / ECI.denom * (100 * 4/2)
  return(ECI)
}

mlr_eci_list <- vector("list", length = 5)
rf_eci_list <- vector("list", length = 5)
lr_pc_eci_list <- vector("list", length = 5)
svm_pc_eci_list <- vector("list", length = 5)
lr_seq_eci_list <- vector("list", length = 5)
svm_seq_eci_list <- vector("list", length = 5)

for (i in 1:5) {
  test_new <- test_data_list[[i]]
  
  mlr_predictions <- mlr_pred[[i]]
  rf_predictions <- rf_pred[[i]]
  lr_pc_predictions <- lr_pc_predictions.list[[i]]
  svm_pc_predictions <- svm_pc_predictions.list[[i]]
  lr_seq_predictions <- lr_seq_predictions1.list[[i]]
  svm_seq_predictions <- svm_seq_predictions.list[[i]]
  
  mlr_eci_list[[i]] <- calculate_ECI(test_new, mlr_predictions)
  print(mlr_eci_list[[i]])
  
  rf_eci_list[[i]] <- calculate_ECI(test_new, rf_predictions)
  print(rf_eci_list[[i]])
  
  lr_pc_eci_list[[i]] <- calculate_ECI(test_new, lr_pc_predictions)
  print(lr_pc_eci_list[[i]])
  
  svm_pc_eci_list[[i]] <- calculate_ECI(test_new, svm_pc_predictions)
  print(svm_pc_eci_list[[i]])
  
  lr_seq_eci_list[[i]] <- calculate_ECI(test_new, lr_seq_predictions)
  print(lr_seq_eci_list[[i]])
  
  svm_seq_eci_list[[i]] <- calculate_ECI(test_new, svm_seq_predictions)
  print(svm_seq_eci_list[[i]])
}

average_mlr_eci <- mean(unlist(mlr_eci_list), na.rm = TRUE)
average_rf_eci <- mean(unlist(rf_eci_list), na.rm = TRUE)
average_lr_pc_eci <- mean(unlist(lr_pc_eci_list), na.rm = TRUE)
average_svm_pc_eci <- mean(unlist(svm_pc_eci_list), na.rm = TRUE)
average_lr_seq_eci <- mean(unlist(lr_seq_eci_list), na.rm = TRUE)
average_svm_seq_eci <- mean(unlist(svm_seq_eci_list), na.rm = TRUE)


#PDI
calculate_PDI <- function(dat.valid.pred, prediction_matrix) {
dat.valid.pred$outcome <- dat.valid.pred$AKI_num
pred_df <- as.data.frame(prediction_matrix)
colnames(pred_df) <- paste0("p", 1:ncol(pred_df))
dat.valid.pred <- cbind(dat.valid.pred, pred_df)
data <- dplyr::select(dat.valid.pred, outcome, p1, p2, p3, p4)

y <- data$outcome
ymin <- min(y)
ymax <- max(y)
noutcome <- ymax - ymin
p <- prod(table(y))
pdi <- c()
for (i in 1:(noutcome + 1)) {
  predprob <- data[, (i + 1)]  
  t0 <- table(predprob, y)  
  
  dim1 <- dim(t0)[1]
  dim2 <- dim(t0)[2]
  t <- cbind(t0[, i], t0[, -i]) 
  
  restrictt <- if (noutcome == 1) { matrix(t[, 2:(noutcome + 1)], ncol = 1) } else { t[, 2:(noutcome + 1)] }  

  c <- apply(restrictt, 2, cumsum)  
  cnew <- if (noutcome == 1) { rbind(rep(0, noutcome), matrix(c[1:(dim(c)[1] - 1), ], ncol = )) } else { rbind(rep(0, noutcome), c[1:(dim(c)[1] - 1), ]) }  
  
  mat <- c()  
  for (j in 1:noutcome) {
    mat0 <- cbind(mat, 0)
    mat1 <- cbind(mat, 1)
    mat <- rbind(mat0, mat1)
  }
  
  r <- 0
  for (k in 1:dim(mat)[1]) {
    dt <- t(apply(restrictt, 1, function(x) mat[k,] * x))
    dcnew <- t(apply(cnew, 1, function(x) (1 - mat[k,]) * x))
    dfinal <- if (noutcome == 1) { cbind(t[, 1], t(dt + dcnew)) } else { cbind(t[, 1], dt + dcnew) }
    r <- r + sum(apply(dfinal, 1, prod)) / (1 + sum(mat[k,]))  
  }
  
  r <- r / p
  pdi <- rbind(pdi, r)
}
print(pdi)
pdi <- rbind(mean(pdi), pdi)
return(as.numeric(pdi[1]))}

mlr_pdi_list <- vector("list", length = 5)
rf_pdi_list <- vector("list", length = 5)
lr_pc_pdi_list <- vector("list", length = 5)
svm_pc_pdi_list <- vector("list", length = 5)
lr_seq_pdi_list <- vector("list", length = 5)
svm_seq_pdi_list <- vector("list", length = 5)

for (i in 1:5) {
  test_new <- test_data_list[[i]]
  
  mlr_predictions <- mlr_pred[[i]]
  rf_predictions <- rf_pred[[i]]
  lr_pc_predictions <- lr_pc_predictions.list[[i]]
  svm_pc_predictions <- svm_pc_predictions.list[[i]]
  lr_seq_predictions <- lr_seq_predictions1.list[[i]]
  svm_seq_predictions <- svm_seq_predictions.list[[i]]
  
  mlr_pdi_list[[i]] <- calculate_PDI(test_new, mlr_predictions)
  print(mlr_pdi_list[[i]])
  
  rf_pdi_list[[i]] <- calculate_PDI(test_new, rf_predictions)
  print(rf_pdi_list[[i]])
  
  lr_pc_pdi_list[[i]] <- calculate_PDI(test_new, lr_pc_predictions)
  print(lr_pc_pdi_list[[i]])
  
  svm_pc_pdi_list[[i]] <- calculate_PDI(test_new, svm_pc_predictions)
  print(svm_pc_pdi_list[[i]])
  
  lr_seq_pdi_list[[i]] <- calculate_PDI(test_new, lr_seq_predictions)
  print(lr_seq_pdi_list[[i]])
  
  svm_seq_pdi_list[[i]] <- calculate_PDI(test_new, svm_seq_predictions)
  print(svm_seq_pdi_list[[i]])
}

average_mlr_pdi <- mean(unlist(mlr_pdi_list), na.rm = TRUE)
average_rf_pdi <- mean(unlist(rf_pdi_list), na.rm = TRUE)
average_lr_pc_pdi <- mean(unlist(lr_pc_pdi_list), na.rm = TRUE)
average_svm_pc_pdi <- mean(unlist(svm_pc_pdi_list), na.rm = TRUE)
average_lr_seq_pdi <- mean(unlist(lr_seq_pdi_list), na.rm = TRUE)
average_svm_seq_pdi <- mean(unlist(svm_seq_pdi_list), na.rm = TRUE)