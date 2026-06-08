############################################################
# Supplementary R script
# Phenotypic dependency networks of fruit traits in
# mamey sapote and star apple germplasm
#
# Journal target: Genetic Resources and Crop Evolution, Springer
#
# This script reproduces the complete analytical workflow:
# 1. data import and preparation;
# 2. descriptive statistics;
# 3. normality assessment and transformation;
# 4. Bayesian network structure learning;
# 5. leave-one-out jackknife stability analysis;
# 6. inspection of Markov equivalence classes through CPDAGs;
# 7. fitting of the final retained structural equation models;
# 8. extraction of path coefficients and Wald tests.
#
# The SEMs fitted at the end of this script are the final retained
# models defined after evaluating the Bayesian network structures,
# jackknife stability, CPDAG information, and biological plausibility.
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
# 3. File paths and output folder
############################

# The input files must be available in the working directory.
# Adjust these file names or paths if necessary.
file_star_apple <- "caimito.txt"
file_mamey      <- "sapote.txt"

output_dir <- "outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

############################
# 4. Helper functions
############################

required_traits <- c("FRW", "FRL", "FRD", "PET", "PEW", "SEL", "SED", "SEW")

read_star_apple_data <- function(file) {
  dat <- read.table(file, header = TRUE)
  
  if (!all(required_traits %in% names(dat))) {
    stop(
      "The star apple dataset must contain the variables: ",
      paste(required_traits, collapse = ", ")
    )
  }
  
  dat <- dat[, required_traits]
  dat <- as.data.frame(dat)
  return(dat)
}

read_mamey_data <- function(file) {
  dat <- read.table(file, header = TRUE)
  
  # If the file already contains the standardized variable names,
  # select the eight variables directly.
  if (all(required_traits %in% names(dat))) {
    dat <- dat[, required_traits]
    dat <- as.data.frame(dat)
    return(dat)
  }
  
  # The original mamey sapote file contains 21 columns.
  # The following names reproduce the variable coding used in the article.
  if (ncol(dat) == 21) {
    colnames(dat) <- c(
      "ID", "FRW", "FRL", "FRD", "PUT", "PET", "PEW",
      "PUP", "NSE", "SEL", "SED", "SEW", "LEL", "LEW",
      "TRH", "TTD", "TCD", "PRO", "SAC", "GLU", "FRU"
    )
  }
  
  if (!all(required_traits %in% names(dat))) {
    stop(
      "The mamey sapote dataset could not be mapped to the variables: ",
      paste(required_traits, collapse = ", ")
    )
  }
  
  dat <- dat[, required_traits]
  dat <- as.data.frame(dat)
  return(dat)
}

make_desc <- function(dat) {
  d <- psych::describe(dat)
  
  out <- data.frame(
    Variable = names(dat),
    N        = d$n,
    Mean     = round(d$mean, 2),
    SE       = round(d$se, 2),
    SD       = round(d$sd, 2),
    Min      = d$min,
    Max      = d$max,
    CV       = round((d$sd / d$mean) * 100, 2)
  )
  
  rownames(out) <- NULL
  return(out)
}

export_arcs <- function(bn_obj, file_name) {
  arc_df <- as.data.frame(arcs(bn_obj))
  write.csv(arc_df, file.path(output_dir, file_name), row.names = FALSE)
}

write_matrix_csv <- function(mat, file_name) {
  write.csv(round(mat, 2), file.path(output_dir, file_name))
}

############################
# 5. Data import and preparation
############################

# Star apple is analyzed using the file containing the standardized
# variable names FRW, FRL, FRD, PET, PEW, SEL, SED, and SEW.
dataC <- read_star_apple_data(file_star_apple)

# Mamey sapote may be supplied either with standardized variable names
# or with the original 21-column coding from the source dataset.
dataS <- read_mamey_data(file_mamey)

cat("\nDimensions of star apple dataset:\n")
print(dim(dataC))

cat("\nDimensions of mamey sapote dataset:\n")
print(dim(dataS))

cat("\nColumn names of star apple dataset:\n")
print(names(dataC))

cat("\nColumn names of mamey sapote dataset:\n")
print(names(dataS))

# Natural logarithmic transformations for star apple variables.
# These transformations were used to improve agreement with the
# Gaussian modeling framework.
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

############################
# 6. Descriptive statistics
############################

# Descriptive statistics are computed on the original measurement scale.
descS <- make_desc(dataS)
descC <- make_desc(dataC)

cat("\nDescriptive statistics - Mamey sapote:\n")
print(descS)

cat("\nDescriptive statistics - Star apple:\n")
print(descC)

write.csv(descS, file.path(output_dir, "descriptive_mamey_sapote.csv"), row.names = FALSE)
write.csv(descC, file.path(output_dir, "descriptive_star_apple.csv"), row.names = FALSE)

############################
# 7. Normality assessment
############################

# Mamey sapote was analyzed on the original measurement scale.
cat("\nUnivariate and multivariate normality - Mamey sapote:\n")
normS <- MVN::mvn(
  data = dataS,
  mvn_test = "royston",
  univariate_test = "SW"
)
print(normS$multivariate_normality)
print(normS$univariate_normality)

# For documentation, normality is also assessed for star apple before
# and after transformation.
cat("\nUnivariate and multivariate normality - Star apple original scale:\n")
normC_original <- MVN::mvn(
  data = dataC,
  mvn_test = "royston",
  univariate_test = "SW"
)
print(normC_original$multivariate_normality)
print(normC_original$univariate_normality)

cat("\nUnivariate and multivariate normality - Star apple after transformation:\n")
normC <- MVN::mvn(
  data = dataC_t,
  mvn_test = "royston",
  univariate_test = "SW"
)
print(normC$multivariate_normality)
print(normC$univariate_normality)

write.csv(normS$multivariate_normality,
          file.path(output_dir, "normality_mamey_multivariate.csv"),
          row.names = FALSE)
write.csv(normS$univariate_normality,
          file.path(output_dir, "normality_mamey_univariate.csv"),
          row.names = FALSE)
write.csv(normC_original$multivariate_normality,
          file.path(output_dir, "normality_star_original_multivariate.csv"),
          row.names = FALSE)
write.csv(normC_original$univariate_normality,
          file.path(output_dir, "normality_star_original_univariate.csv"),
          row.names = FALSE)
write.csv(normC$multivariate_normality,
          file.path(output_dir, "normality_star_transformed_multivariate.csv"),
          row.names = FALSE)
write.csv(normC$univariate_normality,
          file.path(output_dir, "normality_star_transformed_univariate.csv"),
          row.names = FALSE)

############################
# 8. Bayesian network structure learning
############################

# IAMB is a constraint based structure learning algorithm. The test = "cor"
# option applies conditional independence tests based on partial correlations.
bnS_iamb_005 <- iamb(dataS,   test = "cor", alpha = 0.05)
bnC_iamb_005 <- iamb(dataC_t, test = "cor", alpha = 0.05)

bnS_iamb_020 <- iamb(dataS,   test = "cor", alpha = 0.20)
bnC_iamb_020 <- iamb(dataC_t, test = "cor", alpha = 0.20)

# Tabu search is a score based structure learning algorithm. The BGe score
# is appropriate for continuous variables under a Gaussian modeling framework.
bnS_tabu <- tabu(dataS,   score = "bge")
bnC_tabu <- tabu(dataC_t, score = "bge")

cat("\nLearned network - Mamey sapote - IAMB alpha = 0.05:\n")
print(arcs(bnS_iamb_005))

cat("\nLearned network - Mamey sapote - IAMB alpha = 0.20:\n")
print(arcs(bnS_iamb_020))

cat("\nLearned network - Mamey sapote - Tabu search:\n")
print(arcs(bnS_tabu))

cat("\nLearned network - Star apple - IAMB alpha = 0.05:\n")
print(arcs(bnC_iamb_005))

cat("\nLearned network - Star apple - IAMB alpha = 0.20:\n")
print(arcs(bnC_iamb_020))

cat("\nLearned network - Star apple - Tabu search:\n")
print(arcs(bnC_tabu))

export_arcs(bnS_iamb_005, "arcs_mamey_iamb_alpha005.csv")
export_arcs(bnS_iamb_020, "arcs_mamey_iamb_alpha020.csv")
export_arcs(bnS_tabu,     "arcs_mamey_tabu_bge.csv")

export_arcs(bnC_iamb_005, "arcs_star_iamb_alpha005.csv")
export_arcs(bnC_iamb_020, "arcs_star_iamb_alpha020.csv")
export_arcs(bnC_tabu,     "arcs_star_tabu_bge.csv")

############################
# 9. Optional plots of learned structures
############################

# Plotting Bayesian networks with graphviz.plot requires Rgraphviz.
# If Rgraphviz is not available, this step is skipped.
if (requireNamespace("Rgraphviz", quietly = TRUE)) {
  
  pdf(file.path(output_dir, "learned_networks.pdf"), width = 12, height = 8)
  par(mfrow = c(2, 3))
  
  graphviz.plot(bnS_iamb_005,
                main = expression(paste("Mamey sapote - IAMB (", alpha, " = 0.05)")))
  graphviz.plot(bnS_iamb_020,
                main = expression(paste("Mamey sapote - IAMB (", alpha, " = 0.20)")))
  graphviz.plot(bnS_tabu,
                main = "Mamey sapote - Tabu search")
  
  graphviz.plot(bnC_iamb_005,
                main = expression(paste("Star apple - IAMB (", alpha, " = 0.05)")))
  graphviz.plot(bnC_iamb_020,
                main = expression(paste("Star apple - IAMB (", alpha, " = 0.20)")))
  graphviz.plot(bnC_tabu,
                main = "Star apple - Tabu search")
  
  dev.off()
  
} else {
  message("Rgraphviz is not available. Bayesian network plots were skipped.")
}

############################
# 10. Leave-one-out jackknife stability analysis
############################

jackknife_bn <- function(dat, learner = c("iamb", "tabu"), alpha = 0.05) {
  
  learner <- match.arg(learner)
  p <- ncol(dat)
  
  arc_count        <- matrix(0, nrow = p, ncol = p)
  connection_count <- matrix(0, nrow = p, ncol = p)
  direction_count  <- matrix(0, nrow = p, ncol = p)
  
  colnames(arc_count)        <- colnames(dat)
  rownames(arc_count)        <- colnames(dat)
  colnames(connection_count) <- colnames(dat)
  rownames(connection_count) <- colnames(dat)
  colnames(direction_count)  <- colnames(dat)
  rownames(direction_count)  <- colnames(dat)
  
  for (i in seq_len(nrow(dat))) {
    
    dat_i <- dat[-i, , drop = FALSE]
    
    bn_i <- switch(
      learner,
      iamb = iamb(dat_i, test = "cor", alpha = alpha),
      tabu = tabu(dat_i, score = "bge")
    )
    
    amat_i <- amat(bn_i)
    
    # The raw arc matrix is useful to reproduce the orientation-specific
    # output returned by bnlearn.
    arc_count <- arc_count + amat_i
    
    # A connection is counted whenever two nodes are adjacent, regardless of
    # whether the edge is directed or undirected in that resample.
    connection_i <- (amat_i + t(amat_i)) > 0
    connection_count <- connection_count + connection_i
    
    # A direction is counted only when an arc is strictly oriented in one
    # direction and not represented as an undirected edge.
    direction_i <- (amat_i == 1) & (t(amat_i) == 0)
    direction_count <- direction_count + direction_i
  }
  
  arc_percent        <- 100 * arc_count / nrow(dat)
  connection_percent <- 100 * connection_count / nrow(dat)
  direction_percent  <- 100 * direction_count / nrow(dat)
  
  diag(arc_percent)        <- 0
  diag(connection_percent) <- 0
  diag(direction_percent)  <- 0
  
  return(list(
    arc_count          = arc_count,
    connection_count   = connection_count,
    direction_count    = direction_count,
    arc_percent        = arc_percent,
    connection_percent = connection_percent,
    direction_percent  = direction_percent
  ))
}

extract_stable_components <- function(jk, species, method, threshold = 80,
                                      mode = c("connection_or_direction", "direction_only")) {
  
  mode <- match.arg(mode)
  vars <- colnames(jk$connection_percent)
  p <- length(vars)
  
  out <- list()
  k <- 1
  
  for (i in seq_len(p - 1)) {
    for (j in (i + 1):p) {
      
      connection_support <- jk$connection_percent[i, j]
      direction_ij <- jk$direction_percent[i, j]
      direction_ji <- jk$direction_percent[j, i]
      direction_support <- max(direction_ij, direction_ji)
      
      if (direction_support >= threshold) {
        
        if (direction_ij >= direction_ji) {
          component <- paste(vars[i], "->", vars[j])
        } else {
          component <- paste(vars[j], "->", vars[i])
        }
        
        out[[k]] <- data.frame(
          Species = species,
          Method = method,
          Component = component,
          Support = round(direction_support, 2),
          SupportType = "direction recovery",
          stringsAsFactors = FALSE
        )
        k <- k + 1
        
      } else if (mode == "connection_or_direction" &&
                 connection_support >= threshold) {
        
        component <- paste(vars[i], "--", vars[j])
        
        out[[k]] <- data.frame(
          Species = species,
          Method = method,
          Component = component,
          Support = round(connection_support, 2),
          SupportType = "connection occurrence",
          stringsAsFactors = FALSE
        )
        k <- k + 1
      }
    }
  }
  
  if (length(out) == 0) {
    return(data.frame(
      Species = character(0),
      Method = character(0),
      Component = character(0),
      Support = numeric(0),
      SupportType = character(0)
    ))
  }
  
  out <- do.call(rbind, out)
  out <- out[order(out$Species, out$Method, -out$Support), ]
  rownames(out) <- NULL
  return(out)
}

# Run jackknife for IAMB with alpha = 0.05.
jkS_iamb_005 <- jackknife_bn(dataS,   learner = "iamb", alpha = 0.05)
jkC_iamb_005 <- jackknife_bn(dataC_t, learner = "iamb", alpha = 0.05)

# Run jackknife for IAMB with alpha = 0.20.
jkS_iamb_020 <- jackknife_bn(dataS,   learner = "iamb", alpha = 0.20)
jkC_iamb_020 <- jackknife_bn(dataC_t, learner = "iamb", alpha = 0.20)

# Run jackknife for Tabu search.
jkS_tabu <- jackknife_bn(dataS,   learner = "tabu")
jkC_tabu <- jackknife_bn(dataC_t, learner = "tabu")

cat("\nJackknife connection occurrence (%) - Mamey sapote - IAMB alpha = 0.05:\n")
print(round(jkS_iamb_005$connection_percent, 2))
cat("\nJackknife direction recovery (%) - Mamey sapote - IAMB alpha = 0.05:\n")
print(round(jkS_iamb_005$direction_percent, 2))

cat("\nJackknife connection occurrence (%) - Mamey sapote - IAMB alpha = 0.20:\n")
print(round(jkS_iamb_020$connection_percent, 2))
cat("\nJackknife direction recovery (%) - Mamey sapote - IAMB alpha = 0.20:\n")
print(round(jkS_iamb_020$direction_percent, 2))

cat("\nJackknife connection occurrence (%) - Mamey sapote - Tabu search:\n")
print(round(jkS_tabu$connection_percent, 2))
cat("\nJackknife direction recovery (%) - Mamey sapote - Tabu search:\n")
print(round(jkS_tabu$direction_percent, 2))

cat("\nJackknife connection occurrence (%) - Star apple - IAMB alpha = 0.05:\n")
print(round(jkC_iamb_005$connection_percent, 2))
cat("\nJackknife direction recovery (%) - Star apple - IAMB alpha = 0.05:\n")
print(round(jkC_iamb_005$direction_percent, 2))

cat("\nJackknife connection occurrence (%) - Star apple - IAMB alpha = 0.20:\n")
print(round(jkC_iamb_020$connection_percent, 2))
cat("\nJackknife direction recovery (%) - Star apple - IAMB alpha = 0.20:\n")
print(round(jkC_iamb_020$direction_percent, 2))

cat("\nJackknife connection occurrence (%) - Star apple - Tabu search:\n")
print(round(jkC_tabu$connection_percent, 2))
cat("\nJackknife direction recovery (%) - Star apple - Tabu search:\n")
print(round(jkC_tabu$direction_percent, 2))

# Export complete jackknife matrices.
write_matrix_csv(jkS_iamb_005$arc_percent,        "jk_mamey_iamb005_raw_arc_percent.csv")
write_matrix_csv(jkS_iamb_005$connection_percent, "jk_mamey_iamb005_connection_percent.csv")
write_matrix_csv(jkS_iamb_005$direction_percent,  "jk_mamey_iamb005_direction_percent.csv")

write_matrix_csv(jkS_iamb_020$arc_percent,        "jk_mamey_iamb020_raw_arc_percent.csv")
write_matrix_csv(jkS_iamb_020$connection_percent, "jk_mamey_iamb020_connection_percent.csv")
write_matrix_csv(jkS_iamb_020$direction_percent,  "jk_mamey_iamb020_direction_percent.csv")

write_matrix_csv(jkS_tabu$arc_percent,        "jk_mamey_tabu_raw_arc_percent.csv")
write_matrix_csv(jkS_tabu$connection_percent, "jk_mamey_tabu_connection_percent.csv")
write_matrix_csv(jkS_tabu$direction_percent,  "jk_mamey_tabu_direction_percent.csv")

write_matrix_csv(jkC_iamb_005$arc_percent,        "jk_star_iamb005_raw_arc_percent.csv")
write_matrix_csv(jkC_iamb_005$connection_percent, "jk_star_iamb005_connection_percent.csv")
write_matrix_csv(jkC_iamb_005$direction_percent,  "jk_star_iamb005_direction_percent.csv")

write_matrix_csv(jkC_iamb_020$arc_percent,        "jk_star_iamb020_raw_arc_percent.csv")
write_matrix_csv(jkC_iamb_020$connection_percent, "jk_star_iamb020_connection_percent.csv")
write_matrix_csv(jkC_iamb_020$direction_percent,  "jk_star_iamb020_direction_percent.csv")

write_matrix_csv(jkC_tabu$arc_percent,        "jk_star_tabu_raw_arc_percent.csv")
write_matrix_csv(jkC_tabu$connection_percent, "jk_star_tabu_connection_percent.csv")
write_matrix_csv(jkC_tabu$direction_percent,  "jk_star_tabu_direction_percent.csv")

# Components reported in the manuscript table.
# For IAMB, connection occurrence is emphasized when direction recovery is limited.
# For Tabu search, only directions with recovery of at least 80% are reported.
stable_components <- rbind(
  extract_stable_components(jkS_iamb_005, "Mamey sapote", "IAMB alpha = 0.05",
                            threshold = 80, mode = "connection_or_direction"),
  extract_stable_components(jkS_iamb_020, "Mamey sapote", "IAMB alpha = 0.20",
                            threshold = 80, mode = "connection_or_direction"),
  extract_stable_components(jkS_tabu, "Mamey sapote", "Tabu search",
                            threshold = 80, mode = "direction_only"),
  extract_stable_components(jkC_iamb_005, "Star apple", "IAMB alpha = 0.05",
                            threshold = 80, mode = "connection_or_direction"),
  extract_stable_components(jkC_iamb_020, "Star apple", "IAMB alpha = 0.20",
                            threshold = 80, mode = "connection_or_direction"),
  extract_stable_components(jkC_tabu, "Star apple", "Tabu search",
                            threshold = 80, mode = "direction_only")
)

cat("\nStable jackknife components with support >= 80%:\n")
print(stable_components)

write.csv(stable_components,
          file.path(output_dir, "stable_jackknife_components_threshold80.csv"),
          row.names = FALSE)

############################
# 11. Markov equivalence classes through CPDAGs
############################

cpdagS <- cpdag(bnS_tabu)
cpdagC <- cpdag(bnC_tabu)

cat("\nCPDAG arcs - Mamey sapote:\n")
print(arcs(cpdagS))

cat("\nCPDAG arcs - Star apple:\n")
print(arcs(cpdagC))

export_arcs(cpdagS, "cpdag_arcs_mamey_tabu_bge.csv")
export_arcs(cpdagC, "cpdag_arcs_star_tabu_bge.csv")

if (requireNamespace("Rgraphviz", quietly = TRUE)) {
  
  pdf(file.path(output_dir, "cpdag_structures.pdf"), width = 10, height = 5)
  par(mfrow = c(1, 2))
  
  graphviz.plot(cpdagS, main = "CPDAG - Mamey sapote")
  graphviz.plot(cpdagC, main = "CPDAG - Star apple")
  
  dev.off()
  
} else {
  message("Rgraphviz is not available. CPDAG plots were skipped.")
}

############################
# 12. Helper functions for SEM
############################

build_sem_model <- function(paths, variables) {
  
  path_lines <- paste0(
    paths, ", b", sprintf("%02d", seq_along(paths)), ", NA"
  )
  
  variance_lines <- paste0(
    variables, " <-> ", variables, ", v", sprintf("%02d", seq_along(variables)), ", NA"
  )
  
  model_text <- paste(
    c(path_lines, "", variance_lines),
    collapse = "\n"
  )
  
  return(model_text)
}

run_sem_model <- function(model_text, data_mat, model_name = "SEM model") {
  
  model_obj <- specifyModel(textConnection(model_text))
  
  # SEMs are fitted to the correlation matrix. Therefore, the path
  # coefficients among observed variables are interpreted as standardized
  # coefficients.
  fit <- sem(model_obj, cor(data_mat), N = nrow(data_mat))
  
  cat("\n====================================\n")
  cat(model_name, "\n")
  cat("====================================\n")
  print(summary(fit, fit.indices = c("GFI", "AGFI", "RMSEA")))
  
  return(fit)
}

extract_fit_indices <- function(fit, model_name) {
  
  txt <- capture.output(summary(fit, fit.indices = c("GFI", "AGFI", "RMSEA")))
  
  line_chi   <- txt[grep("Model Chisquare", txt)]
  line_gfi   <- txt[grep("Goodness-of-fit index", txt)]
  line_agfi  <- txt[grep("Adjusted goodness-of-fit index", txt)]
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

extract_sem_parameters <- function(fit, paths, species) {
  
  estimates <- coef(fit)
  vc <- tryCatch(vcov(fit), error = function(e) NULL)
  
  if (is.null(vc)) {
    stop("The variance-covariance matrix of the SEM parameters could not be extracted.")
  }
  
  se <- sqrt(diag(vc))
  z_value <- estimates / se
  p_value <- 2 * (1 - pnorm(abs(z_value)))
  
  parameter_table <- data.frame(
    Parameter = names(estimates),
    Estimate = as.numeric(estimates),
    SE = as.numeric(se),
    WaldZ = as.numeric(z_value),
    PValue = as.numeric(p_value),
    stringsAsFactors = FALSE
  )
  
  path_map <- data.frame(
    Parameter = paste0("b", sprintf("%02d", seq_along(paths))),
    Path = paths,
    stringsAsFactors = FALSE
  )
  
  path_table <- merge(path_map, parameter_table, by = "Parameter", all.x = TRUE)
  path_table$Species <- species
  path_table <- path_table[, c("Species", "Path", "Parameter", "Estimate", "SE", "WaldZ", "PValue")]
  
  path_table$Estimate <- round(path_table$Estimate, 4)
  path_table$SE       <- round(path_table$SE, 4)
  path_table$WaldZ    <- round(path_table$WaldZ, 4)
  path_table$PValue   <- signif(path_table$PValue, 4)
  
  return(path_table)
}

save_sem_summary <- function(fit, file_name) {
  txt <- capture.output(summary(fit, fit.indices = c("GFI", "AGFI", "RMSEA")))
  writeLines(txt, file.path(output_dir, file_name))
}

save_modification_indices <- function(fit, file_name) {
  txt <- capture.output(modIndices(fit))
  writeLines(txt, file.path(output_dir, file_name))
}

############################
# 13. Final retained SEM specifications
############################

sem_vars <- required_traits

# Final retained SEM for mamey sapote.
# This specification was defined after evaluating the updated Bayesian
# network structures, jackknife stability, CPDAG information, and
# biological plausibility.
mamey_final_paths <- c(
  "FRW -> PET",
  "FRW -> PEW",
  "FRW -> FRD",
  "PET -> PEW",
  "FRL -> FRD",
  "PEW -> FRL",
  "FRL -> SEL",
  "SEL -> SEW",
  "SEW -> SED"
)

# Final retained SEM for star apple.
# The FRW -> SED path is retained because it belongs to the stable
# score based backbone under jackknife resampling and is coherent with
# the final SEM diagnostics.
star_final_paths <- c(
  "FRW -> PET",
  "FRW -> PEW",
  "FRW -> FRD",
  "FRW -> FRL",
  "PET -> PEW",
  "FRL -> PEW",
  "PEW -> FRD",
  "FRD -> SEW",
  "SEW -> SED",
  "SED -> SEL",
  "FRW -> SED"
)

mamey_final_model <- build_sem_model(
  paths = mamey_final_paths,
  variables = sem_vars
)

star_final_model <- build_sem_model(
  paths = star_final_paths,
  variables = sem_vars
)

cat("\nFinal retained SEM model - Mamey sapote:\n")
cat(mamey_final_model, "\n")

cat("\nFinal retained SEM model - Star apple:\n")
cat(star_final_model, "\n")

writeLines(mamey_final_model, file.path(output_dir, "sem_model_mamey_final.txt"))
writeLines(star_final_model,  file.path(output_dir, "sem_model_star_final.txt"))

############################
# 14. Fit final retained SEMs
############################

fit_mamey_final <- run_sem_model(
  model_text = mamey_final_model,
  data_mat = dataS,
  model_name = "Mamey sapote final retained SEM"
)

fit_star_final <- run_sem_model(
  model_text = star_final_model,
  data_mat = dataC_t,
  model_name = "Star apple final retained SEM"
)

fit_summary_final <- rbind(
  extract_fit_indices(fit_mamey_final, "Mamey sapote final retained SEM"),
  extract_fit_indices(fit_star_final,  "Star apple final retained SEM")
)

fit_summary_final$Chisq  <- round(fit_summary_final$Chisq, 3)
fit_summary_final$PValue <- round(fit_summary_final$PValue, 6)
fit_summary_final$GFI    <- round(fit_summary_final$GFI, 3)
fit_summary_final$AGFI   <- round(fit_summary_final$AGFI, 3)
fit_summary_final$RMSEA  <- round(fit_summary_final$RMSEA, 3)

cat("\nFinal SEM model fit summary:\n")
print(fit_summary_final)

write.csv(fit_summary_final,
          file.path(output_dir, "final_sem_model_fit_summary.csv"),
          row.names = FALSE)

save_sem_summary(fit_mamey_final, "sem_summary_mamey_final.txt")
save_sem_summary(fit_star_final,  "sem_summary_star_final.txt")

############################
# 15. Path coefficients and Wald tests
############################

mamey_path_tests <- extract_sem_parameters(
  fit = fit_mamey_final,
  paths = mamey_final_paths,
  species = "Mamey sapote"
)

star_path_tests <- extract_sem_parameters(
  fit = fit_star_final,
  paths = star_final_paths,
  species = "Star apple"
)

path_tests <- rbind(mamey_path_tests, star_path_tests)

cat("\nPath coefficients and Wald tests - Final retained SEMs:\n")
print(path_tests)

write.csv(mamey_path_tests,
          file.path(output_dir, "path_coefficients_wald_mamey_final.csv"),
          row.names = FALSE)
write.csv(star_path_tests,
          file.path(output_dir, "path_coefficients_wald_star_final.csv"),
          row.names = FALSE)
write.csv(path_tests,
          file.path(output_dir, "path_coefficients_wald_final_models.csv"),
          row.names = FALSE)

############################
# 16. Modification indices for final retained SEMs
############################

cat("\nModification indices - Mamey sapote final retained SEM:\n")
mi_mamey_final <- modIndices(fit_mamey_final)
print(mi_mamey_final)

cat("\nModification indices - Star apple final retained SEM:\n")
mi_star_final <- modIndices(fit_star_final)
print(mi_star_final)

save_modification_indices(fit_mamey_final, "modification_indices_mamey_final.txt")
save_modification_indices(fit_star_final,  "modification_indices_star_final.txt")

############################
# 17. Optional SEM path diagrams
############################

# SEM path diagrams are optional and require the semPlot package.
# If semPlot is not available, this step is skipped.
if (requireNamespace("semPlot", quietly = TRUE)) {
  
  library(semPlot)
  
  pdf(file.path(output_dir, "final_sem_path_diagrams.pdf"), width = 8, height = 10)
  par(mfrow = c(2, 1))
  
  semPaths(
    fit_mamey_final,
    what = "paths",
    whatLabels = "est",
    residuals = FALSE,
    style = "ram",
    layout = "tree2",
    edge.label.cex = 0.8
  )
  title("Mamey sapote - final retained SEM", line = 2)
  
  semPaths(
    fit_star_final,
    what = "paths",
    whatLabels = "est",
    residuals = FALSE,
    style = "ram",
    layout = "tree2",
    edge.label.cex = 0.8
  )
  title("Star apple - final retained SEM", line = 2)
  
  dev.off()
  
} else {
  message("semPlot is not available. SEM path diagrams were skipped.")
}

############################
# 18. Session information
############################

cat("\nSession information:\n")
print(sessionInfo())

cat("\nPackage version - bnlearn:\n")
print(packageVersion("bnlearn"))

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "session_info.txt"))
writeLines(capture.output(packageVersion("bnlearn")),
           file.path(output_dir, "bnlearn_version.txt"))

cat("\nAnalysis completed successfully.\n")
cat("Results were saved in the 'outputs' folder.\n")
