############################################################
# Supplementary R Script
# Phenotypic trait relationships in mamey sapote and star apple
# revealed by Bayesian networks and structural equation modelling
############################################################

############################
# 1. Reproducibility
############################
set.seed(123)

############################
# 2. Required packages
############################
required_packages <- c("bnlearn", "MVN", "sem", "psych")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(bnlearn)
library(MVN)
library(sem)
library(psych)

############################
# 3. File paths
############################
# Adjust these file names/paths if necessary
file_star_apple <- "caimito.txt"
file_mamey      <- "sapote.txt"

############################
# 4. Data import
############################

# --- Star apple dataset ---
dataC_raw <- read.table(file_star_apple, header = TRUE)

# Keep the same variables used in the original analysis
dataC <- data.frame(dataC_raw[-c(1, 7, 11:14)])

# --- Mamey sapote dataset ---
dataS0 <- read.table(file_mamey, header = TRUE)

colnames(dataS0) <- c(
  "ID","FRW","FRL","FRD","PUT","PET","PEW",
  "PUP","NSE","SEL","SED","SEW","LEL","LEW",
  "TRH","TTD","TCD","PRO","SAC","GLU","FRU"
)

# Keep the same variables used in the original analysis
dataS <- data.frame(dataS0[-c(1, 5, 8, 9, 13:21)])

cat("\nDimensions of star apple dataset:\n")
print(dim(dataC))

cat("\nDimensions of mamey sapote dataset:\n")
print(dim(dataS))

cat("\nColumn names of star apple dataset:\n")
print(names(dataC))

cat("\nColumn names of mamey sapote dataset:\n")
print(names(dataS))

############################
# 5. Descriptive statistics
############################
make_desc <- function(dat) {
  d <- psych::describe(dat)
  out <- data.frame(
    VARIABLE = names(dat),
    MIN      = d$min,
    MAX      = round(d$max, 1),
    MEAN     = round(d$mean, 2),
    SE       = round(d$se, 2),
    CV       = round((d$sd / d$mean) * 100, 2)
  )
  rownames(out) <- NULL
  return(out)
}

descS <- make_desc(dataS)
descC <- make_desc(dataC)

cat("\nDescriptive statistics - Mamey sapote:\n")
print(descS)

cat("\nDescriptive statistics - Star apple:\n")
print(descC)

############################
# 6. Normality assessment
############################

# Log transformation for selected star apple variables,
# following the original analysis
dataC_t <- data.frame(
  FRW = log(dataC$FRW),
  FRL = dataC$FRL,
  FRD = log(dataC$FRD),
  PET = log(dataC$PET),
  PEW = log(dataC$PEW),
  SEL = dataC$SEL,
  SED = log(dataC$SED),
  SEW = log(dataC$SEW)
)

cat("\nUnivariate and multivariate normality - Mamey sapote:\n")
normS <- MVN::mvn(
  data = dataS,
  mvn_test = "royston",
  univariate_test = "SW"
)
print(normS$multivariate_normality)
print(normS$univariate_normality)

cat("\nUnivariate and multivariate normality - Star apple (transformed):\n")
normC <- MVN::mvn(
  data = dataC_t,
  mvn_test = "royston",
  univariate_test = "SW"
)
print(normC$multivariate_normality)
print(normC$univariate_normality)

############################
# 7. Structure learning
############################

# IAMB
bnS_iamb_005 <- iamb(dataS,   test = "cor", alpha = 0.05)
bnC_iamb_005 <- iamb(dataC_t, test = "cor", alpha = 0.05)

bnS_iamb_020 <- iamb(dataS,   test = "cor", alpha = 0.20)
bnC_iamb_020 <- iamb(dataC_t, test = "cor", alpha = 0.20)

# Tabu
bnS_tabu <- tabu(dataS,   score = "bge")
bnC_tabu <- tabu(dataC_t, score = "bge")

cat("\nLearned network - Mamey sapote - IAMB alpha=0.05:\n")
print(arcs(bnS_iamb_005))

cat("\nLearned network - Mamey sapote - IAMB alpha=0.20:\n")
print(arcs(bnS_iamb_020))

cat("\nLearned network - Mamey sapote - Tabu:\n")
print(arcs(bnS_tabu))

cat("\nLearned network - Star apple - IAMB alpha=0.05:\n")
print(arcs(bnC_iamb_005))

cat("\nLearned network - Star apple - IAMB alpha=0.20:\n")
print(arcs(bnC_iamb_020))

cat("\nLearned network - Star apple - Tabu:\n")
print(arcs(bnC_tabu))

############################
# 8. Plot learned structures
############################
par(mfrow = c(2, 3))

plot(
  bnS_iamb_005,
  main = expression(paste("Mamey sapote - IAMB (", alpha, "=0.05)"))
)

plot(
  bnS_iamb_020,
  main = expression(paste("Mamey sapote - IAMB (", alpha, "=0.20)"))
)

plot(
  bnS_tabu,
  main = "Mamey sapote - Tabu search"
)

plot(
  bnC_iamb_005,
  main = expression(paste("Star apple - IAMB (", alpha, "=0.05)"))
)

plot(
  bnC_iamb_020,
  main = expression(paste("Star apple - IAMB (", alpha, "=0.20)"))
)

plot(
  bnC_tabu,
  main = "Star apple - Tabu search"
)

############################
# 9. Jackknife stability
############################

jackknife_bn <- function(dat, learner = c("iamb", "tabu"), alpha = 0.05) {
  learner <- match.arg(learner)
  
  p <- ncol(dat)
  edge_count <- matrix(0, nrow = p, ncol = p)
  dir_count  <- matrix(0, nrow = p, ncol = p)
  
  colnames(edge_count) <- colnames(dat)
  rownames(edge_count) <- colnames(dat)
  colnames(dir_count)  <- colnames(dat)
  rownames(dir_count)  <- colnames(dat)
  
  for (i in seq_len(nrow(dat))) {
    f <- dat[-i, , drop = FALSE]
    
    bn_fit <- switch(
      learner,
      iamb = iamb(f, test = "cor", alpha = alpha),
      tabu = tabu(f, score = "bge")
    )
    
    amat_i <- amat(bn_fit)
    edge_count <- edge_count + amat_i
    
    dir_mat <- (amat_i - t(amat_i)) > 0
    dir_count <- dir_count + dir_mat
  }
  
  edge_percent <- 100 * edge_count / nrow(dat)
  dir_percent  <- 100 * dir_count  / nrow(dat)
  
  return(list(
    edge_count   = edge_count,
    dir_count    = dir_count,
    edge_percent = edge_percent,
    dir_percent  = dir_percent
  ))
}

# Run jackknife for IAMB (alpha = 0.05)
jkS_iamb_005 <- jackknife_bn(dataS,   learner = "iamb", alpha = 0.05)
jkC_iamb_005 <- jackknife_bn(dataC_t, learner = "iamb", alpha = 0.05)

# Run jackknife for IAMB (alpha = 0.20)
jkS_iamb_020 <- jackknife_bn(dataS,   learner = "iamb", alpha = 0.20)
jkC_iamb_020 <- jackknife_bn(dataC_t, learner = "iamb", alpha = 0.20)

# Run jackknife for Tabu
jkS_tabu <- jackknife_bn(dataS,   learner = "tabu")
jkC_tabu <- jackknife_bn(dataC_t, learner = "tabu")

# Print IAMB alpha = 0.05
cat("\nJackknife edge occurrence (%) - Mamey sapote - IAMB alpha=0.05:\n")
print(round(jkS_iamb_005$edge_percent, 2))

cat("\nJackknife edge direction recovery (%) - Mamey sapote - IAMB alpha=0.05:\n")
print(round(jkS_iamb_005$dir_percent, 2))

cat("\nJackknife edge occurrence (%) - Star apple - IAMB alpha=0.05:\n")
print(round(jkC_iamb_005$edge_percent, 2))

cat("\nJackknife edge direction recovery (%) - Star apple - IAMB alpha=0.05:\n")
print(round(jkC_iamb_005$dir_percent, 2))

# Print IAMB alpha = 0.20
cat("\nJackknife edge occurrence (%) - Mamey sapote - IAMB alpha=0.20:\n")
print(round(jkS_iamb_020$edge_percent, 2))

cat("\nJackknife edge direction recovery (%) - Mamey sapote - IAMB alpha=0.20:\n")
print(round(jkS_iamb_020$dir_percent, 2))

cat("\nJackknife edge occurrence (%) - Star apple - IAMB alpha=0.20:\n")
print(round(jkC_iamb_020$edge_percent, 2))

cat("\nJackknife edge direction recovery (%) - Star apple - IAMB alpha=0.20:\n")
print(round(jkC_iamb_020$dir_percent, 2))

# Print Tabu
cat("\nJackknife edge occurrence (%) - Mamey sapote - Tabu:\n")
print(round(jkS_tabu$edge_percent, 2))

cat("\nJackknife edge direction recovery (%) - Mamey sapote - Tabu:\n")
print(round(jkS_tabu$dir_percent, 2))

cat("\nJackknife edge occurrence (%) - Star apple - Tabu:\n")
print(round(jkC_tabu$edge_percent, 2))

cat("\nJackknife edge direction recovery (%) - Star apple - Tabu:\n")
print(round(jkC_tabu$dir_percent, 2))

############################
# 10. Equivalence classes (CPDAG)
############################

cpdagS <- cpdag(bnS_tabu)
cpdagC <- cpdag(bnC_tabu)

cat("\nCPDAG arcs - Mamey sapote:\n")
print(arcs(cpdagS))

cat("\nCPDAG arcs - Star apple:\n")
print(arcs(cpdagC))

par(mfrow = c(1, 2))
plot(cpdagS, main = "CPDAG - Mamey sapote")
plot(cpdagC, main = "CPDAG - Star apple")

############################
# 11. Helper functions for SEM
############################
run_sem_model <- function(model_text, data_mat, model_name = "SEM model") {
  model_obj <- specifyModel(textConnection(model_text))
  fit <- sem(model_obj, cor(data_mat), N = nrow(data_mat))
  
  cat("\n====================================\n")
  cat(model_name, "\n")
  cat("====================================\n")
  print(summary(fit, fit.indices = c("GFI", "AGFI", "RMSEA")))
  
  return(fit)
}

extract_fit_indices <- function(fit, model_name) {
  txt <- capture.output(summary(fit, fit.indices = c("GFI", "AGFI", "RMSEA")))
  
  line_chi <- txt[grep("Model Chisquare", txt)]
  line_gfi <- txt[grep("Goodness-of-fit index", txt)]
  line_agfi <- txt[grep("Adjusted goodness-of-fit index", txt)]
  line_rmsea <- txt[grep("RMSEA index", txt)]
  
  chisq  <- as.numeric(sub(".*Model Chisquare = *([0-9\\.Ee+-]+).*", "\\1", line_chi))
  df     <- as.numeric(sub(".*Df = *([0-9\\.Ee+-]+).*", "\\1", line_chi))
  pvalue <- as.numeric(sub(".*Pr\\(>Chisq\\) = *([0-9\\.Ee+-]+).*", "\\1", line_chi))
  gfi    <- as.numeric(sub(".*Goodness-of-fit index = *([0-9\\.Ee+-]+).*", "\\1", line_gfi))
  agfi   <- as.numeric(sub(".*Adjusted goodness-of-fit index = *([0-9\\.Ee+-]+).*", "\\1", line_agfi))
  rmsea  <- as.numeric(sub(".*RMSEA index = *([0-9\\.Ee+-]+).*", "\\1", line_rmsea))
  
  out <- data.frame(
    Model = model_name,
    Chisq = chisq,
    Df = df,
    PValue = pvalue,
    GFI = gfi,
    AGFI = agfi,
    RMSEA = rmsea
  )
  
  return(out)
}

############################
# 12. SEM candidate models
############################

# 12A. Mamey sapote candidate models

sapote_M0 <- "
FRL -> PEW, lam1, NA
PEW -> PET, lam2, NA
PEW -> FRW, lam3, NA
FRW -> FRD, lam4, NA
FRL -> FRD, lam5, NA
FRL -> SEL, lam6, NA
SEL -> SEW, lam7, NA
SED -> SEW, lam8, NA

FRW <-> FRW, alp1, NA
FRL <-> FRL, alp2, NA
FRD <-> FRD, alp3, NA
PET <-> PET, alp4, NA
PEW <-> PEW, alp5, NA
SEL <-> SEL, alp6, NA
SED <-> SED, alp7, NA
SEW <-> SEW, alp8, NA
"

sapote_M1 <- "
FRL -> PEW, lam1, NA
PEW -> PET, lam2, NA
PEW -> FRW, lam3, NA
FRW -> FRD, lam4, NA
FRL -> FRD, lam5, NA
FRL -> SEL, lam6, NA
SEL -> SEW, lam7, NA
SED -> SEW, lam8, NA
PEW -> SED, lam9, NA

FRW <-> FRW, alp1, NA
FRL <-> FRL, alp2, NA
FRD <-> FRD, alp3, NA
PET <-> PET, alp4, NA
PEW <-> PEW, alp5, NA
SEL <-> SEL, alp6, NA
SED <-> SED, alp7, NA
SEW <-> SEW, alp8, NA
"

sapote_M2 <- "
FRL -> PEW, lam1, NA
PEW -> PET, lam2, NA
PEW -> FRW, lam3, NA
FRW -> FRD, lam4, NA
FRL -> FRD, lam5, NA
FRL -> SEL, lam6, NA
FRD -> SEL, lam7, NA
SEL -> SEW, lam8, NA
SED -> SEW, lam9, NA

FRW <-> FRW, alp1, NA
FRL <-> FRL, alp2, NA
FRD <-> FRD, alp3, NA
PET <-> PET, alp4, NA
PEW <-> PEW, alp5, NA
SEL <-> SEL, alp6, NA
SED <-> SED, alp7, NA
SEW <-> SEW, alp8, NA
"

sapote_M3 <- "
FRL -> PEW, lam1, NA
PEW -> PET, lam2, NA
PEW -> FRW, lam3, NA
PET -> FRW, lam4, NA
FRW -> FRD, lam5, NA
FRL -> FRD, lam6, NA
FRL -> SEL, lam7, NA
SEL -> SEW, lam8, NA
SED -> SEW, lam9, NA

FRW <-> FRW, alp1, NA
FRL <-> FRL, alp2, NA
FRD <-> FRD, alp3, NA
PET <-> PET, alp4, NA
PEW <-> PEW, alp5, NA
SEL <-> SEL, alp6, NA
SED <-> SED, alp7, NA
SEW <-> SEW, alp8, NA
"

sapote_M4 <- "
FRL -> PEW, lam1, NA
PEW -> PET, lam2, NA
PEW -> FRW, lam3, NA
FRW -> FRD, lam4, NA
FRL -> FRD, lam5, NA
FRL -> SEL, lam6, NA
FRD -> SEL, lam7, NA
SEL -> SEW, lam8, NA
SED -> SEW, lam9, NA
PEW -> SED, lam10, NA

FRW <-> FRW, alp1, NA
FRL <-> FRL, alp2, NA
FRD <-> FRD, alp3, NA
PET <-> PET, alp4, NA
PEW <-> PEW, alp5, NA
SEL <-> SEL, alp6, NA
SED <-> SED, alp7, NA
SEW <-> SEW, alp8, NA
"

sapote_final <- "
FRW <- PET, lam11, NA
FRW <- PEW, lam12, NA
FRW -> FRD, lam13, NA
FRL -> FRD, lam14, NA
FRL -> SEL, lam15, NA
FRD -> SEL, lam16, NA
SEL -> SEW, lam17, NA
SED -> SEW, lam18, NA
PEW -> PET, lam19, NA
PEW <- FRL, lam20, NA
PEW -> SED, lam21, NA
SED <-> SED, alp1, NA
SEL <-> SEL, alp2, NA
SEW <-> SEW, alp3, NA
FRW <-> FRW, alp4, NA
FRL <-> FRL, alp5, NA
FRD <-> FRD, alp6, NA
PET <-> PET, alp7, NA
PEW <-> PEW, alp8, NA
"

# 12B. Star apple candidate models

caimito_M0 <- "
FRL -> PEW, lam1, NA
PEW -> PET, lam2, NA
PEW -> FRW, lam3, NA
FRL -> FRW, lam4, NA
FRW -> FRD, lam5, NA
FRW -> SED, lam6, NA
SED -> SEL, lam7, NA
SED -> SEW, lam8, NA

FRW <-> FRW, alp1, NA
FRL <-> FRL, alp2, NA
FRD <-> FRD, alp3, NA
PET <-> PET, alp4, NA
PEW <-> PEW, alp5, NA
SEL <-> SEL, alp6, NA
SED <-> SED, alp7, NA
SEW <-> SEW, alp8, NA
"

caimito_M1 <- "
FRL -> PEW, lam1, NA
PEW -> PET, lam2, NA
PEW -> FRW, lam3, NA
FRL -> FRW, lam4, NA
PET -> FRW, lam5, NA
FRW -> FRD, lam6, NA
FRW -> SED, lam7, NA
SED -> SEL, lam8, NA
SED -> SEW, lam9, NA

FRW <-> FRW, alp1, NA
FRL <-> FRL, alp2, NA
FRD <-> FRD, alp3, NA
PET <-> PET, alp4, NA
PEW <-> PEW, alp5, NA
SEL <-> SEL, alp6, NA
SED <-> SED, alp7, NA
SEW <-> SEW, alp8, NA
"

caimito_M2 <- "
FRL -> PEW, lam1, NA
PEW -> PET, lam2, NA
FRL -> PET, lam3, NA
PEW -> FRW, lam4, NA
FRL -> FRW, lam5, NA
FRW -> FRD, lam6, NA
FRW -> SED, lam7, NA
SED -> SEL, lam8, NA
SED -> SEW, lam9, NA

FRW <-> FRW, alp1, NA
FRL <-> FRL, alp2, NA
FRD <-> FRD, alp3, NA
PET <-> PET, alp4, NA
PEW <-> PEW, alp5, NA
SEL <-> SEL, alp6, NA
SED <-> SED, alp7, NA
SEW <-> SEW, alp8, NA
"

caimito_M3 <- "
FRL -> PEW, lam1, NA
PEW -> PET, lam2, NA
FRL -> PET, lam3, NA
PEW -> FRW, lam4, NA
FRL -> FRW, lam5, NA
PET -> FRW, lam6, NA
FRW -> FRD, lam7, NA
FRW -> SED, lam8, NA
SED -> SEL, lam9, NA
SED -> SEW, lam10, NA

FRW <-> FRW, alp1, NA
FRL <-> FRL, alp2, NA
FRD <-> FRD, alp3, NA
PET <-> PET, alp4, NA
PEW <-> PEW, alp5, NA
SEL <-> SEL, alp6, NA
SED <-> SED, alp7, NA
SEW <-> SEW, alp8, NA
"

caimito_M4 <- "
FRL -> PEW, lam1, NA
PEW -> PET, lam2, NA
FRL -> PET, lam3, NA
PEW -> FRW, lam4, NA
FRL -> FRW, lam5, NA
PET -> FRW, lam6, NA
FRW -> FRD, lam7, NA
FRW -> SED, lam8, NA
FRL -> SEW, lam9, NA
SED -> SEL, lam10, NA
SED -> SEW, lam11, NA

FRW <-> FRW, alp1, NA
FRL <-> FRL, alp2, NA
FRD <-> FRD, alp3, NA
PET <-> PET, alp4, NA
PEW <-> PEW, alp5, NA
SEL <-> SEL, alp6, NA
SED <-> SED, alp7, NA
SEW <-> SEW, alp8, NA
"

caimito_final <- "
FRW <- PET, lam1, NA
FRW <- PEW, lam2, NA
FRW -> FRD, lam3, NA
FRW <- FRL, lam4, NA
FRW -> SED, lam5, NA
PEW -> PET, lam6, NA
PEW <- FRL, lam7, NA
FRL -> PET, lam8, NA
FRL -> SEW, lam9, NA
SED -> SEL, lam10, NA
SED -> SEW, lam11, NA
SED <-> SED, alp1, NA
SEL <-> SEL, alp2, NA
SEW <-> SEW, alp3, NA
FRW <-> FRW, alp4, NA
FRL <-> FRL, alp5, NA
FRD <-> FRD, alp6, NA
PET <-> PET, alp7, NA
PEW <-> PEW, alp8, NA
"

############################
# 13. Run candidate SEM models
############################

fit_sapote_M0 <- run_sem_model(sapote_M0, dataS, "Sapote M0")
fit_sapote_M1 <- run_sem_model(sapote_M1, dataS, "Sapote M1 (+ PEW -> SED)")
fit_sapote_M2 <- run_sem_model(sapote_M2, dataS, "Sapote M2 (+ FRD -> SEL)")
fit_sapote_M3 <- run_sem_model(sapote_M3, dataS, "Sapote M3 (+ PET -> FRW)")
fit_sapote_M4 <- run_sem_model(sapote_M4, dataS, "Sapote M4 (+ PEW -> SED, + FRD -> SEL)")

fit_caimito_M0 <- run_sem_model(caimito_M0, dataC_t, "Star apple M0")
fit_caimito_M1 <- run_sem_model(caimito_M1, dataC_t, "Star apple M1 (+ PET -> FRW)")
fit_caimito_M2 <- run_sem_model(caimito_M2, dataC_t, "Star apple M2 (+ FRL -> PET)")
fit_caimito_M3 <- run_sem_model(caimito_M3, dataC_t, "Star apple M3 (+ PET -> FRW, + FRL -> PET)")
fit_caimito_M4 <- run_sem_model(caimito_M4, dataC_t, "Star apple M4 (+ PET -> FRW, + FRL -> PET, + FRL -> SEW)")

############################
# 14. Summarise candidate model fit
############################
fit_summary <- rbind(
  extract_fit_indices(fit_sapote_M0, "Sapote M0"),
  extract_fit_indices(fit_sapote_M1, "Sapote M1"),
  extract_fit_indices(fit_sapote_M2, "Sapote M2"),
  extract_fit_indices(fit_sapote_M3, "Sapote M3"),
  extract_fit_indices(fit_sapote_M4, "Sapote M4"),
  extract_fit_indices(fit_caimito_M0, "StarApple M0"),
  extract_fit_indices(fit_caimito_M1, "StarApple M1"),
  extract_fit_indices(fit_caimito_M2, "StarApple M2"),
  extract_fit_indices(fit_caimito_M3, "StarApple M3"),
  extract_fit_indices(fit_caimito_M4, "StarApple M4")
)

fit_summary$Chisq  <- round(fit_summary$Chisq, 3)
fit_summary$PValue <- round(fit_summary$PValue, 6)
fit_summary$GFI    <- round(fit_summary$GFI, 3)
fit_summary$AGFI   <- round(fit_summary$AGFI, 3)
fit_summary$RMSEA  <- round(fit_summary$RMSEA, 3)

cat("\nCandidate SEM model fit summary:\n")
print(fit_summary)

############################
# 15. Fit only the excellent final models
############################

modelC_final <- specifyModel(textConnection(caimito_final))
semC_final <- sem(modelC_final, cor(dataC_t), N = nrow(dataC_t))

cat("\nSEM summary - Final excellent model - Star apple:\n")
print(summary(semC_final, fit.indices = c("GFI", "AGFI", "RMSEA")))

cat("\nModification indices - Final excellent model - Star apple:\n")
print(modIndices(semC_final))

modelS_final <- specifyModel(textConnection(sapote_final))
semS_final <- sem(modelS_final, cor(dataS), N = nrow(dataS))

cat("\nSEM summary - Final excellent model - Mamey sapote:\n")
print(summary(semS_final, fit.indices = c("GFI", "AGFI", "RMSEA")))

cat("\nModification indices - Final excellent model - Mamey sapote:\n")
print(modIndices(semS_final))

############################
# 16. Optional SEM path diagrams
############################
if (!requireNamespace("semPlot", quietly = TRUE)) {
  install.packages("semPlot", dependencies = TRUE)
}
library(semPlot)

par(mfrow = c(2, 1))

semPaths(
  semS_final,
  what = "paths",
  whatLabels = "est",
  residuals = FALSE,
  style = "ram",
  layout = "tree2",
  edge.color = "dimgray",
  edge.label.cex = 0.9
)
title("Mamey sapote - final excellent model", line = 3)

semPaths(
  semC_final,
  what = "paths",
  whatLabels = "est",
  residuals = FALSE,
  style = "ram",
  layout = "tree2",
  edge.color = "dimgray",
  edge.label.cex = 0.9
)
title("Star apple - final excellent model", line = 3)

############################
# 17. Export results
############################

write.csv(descS, "descriptive_mamey_sapote.csv", row.names = FALSE)
write.csv(descC, "descriptive_star_apple.csv", row.names = FALSE)

# Jackknife summaries - IAMB alpha = 0.05
write.csv(round(jkS_iamb_005$edge_percent, 2), "jk_mamey_iamb005_edge_percent.csv")
write.csv(round(jkS_iamb_005$dir_percent, 2),  "jk_mamey_iamb005_dir_percent.csv")

write.csv(round(jkC_iamb_005$edge_percent, 2), "jk_star_iamb005_edge_percent.csv")
write.csv(round(jkC_iamb_005$dir_percent, 2),  "jk_star_iamb005_dir_percent.csv")

# Jackknife summaries - IAMB alpha = 0.20
write.csv(round(jkS_iamb_020$edge_percent, 2), "jk_mamey_iamb020_edge_percent.csv")
write.csv(round(jkS_iamb_020$dir_percent, 2),  "jk_mamey_iamb020_dir_percent.csv")

write.csv(round(jkC_iamb_020$edge_percent, 2), "jk_star_iamb020_edge_percent.csv")
write.csv(round(jkC_iamb_020$dir_percent, 2),  "jk_star_iamb020_dir_percent.csv")

# Jackknife summaries - Tabu
write.csv(round(jkS_tabu$edge_percent, 2), "jk_mamey_tabu_edge_percent.csv")
write.csv(round(jkS_tabu$dir_percent, 2),  "jk_mamey_tabu_dir_percent.csv")

write.csv(round(jkC_tabu$edge_percent, 2), "jk_star_tabu_edge_percent.csv")
write.csv(round(jkC_tabu$dir_percent, 2),  "jk_star_tabu_dir_percent.csv")

# Candidate model fit
write.csv(fit_summary, "candidate_sem_model_fit_summary.csv", row.names = FALSE)

############################
# 18. Session information
############################
cat("\nSession information:\n")
print(sessionInfo())

cat("\nAnalysis completed successfully.\n")
