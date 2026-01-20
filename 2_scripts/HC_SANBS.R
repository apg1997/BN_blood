#!/usr/bin/env Rscript

# -------------------------------
# HC_SANBS Analysis Script
# -------------------------------

# Load libraries
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load(tidyverse, bnlearn, dtplyr, parallel, Rgraphviz, htmltools)

# -------------------------------
# Load dataset
# -------------------------------
cat("Loading dataset...\n")
bn_df <- readRDS("bn_df.rds")
n <- nrow(bn_df)
k_bic <- log(n)
cat("Sample size:", n, "BIC penalty:", k_bic, "\n\n")

# -------------------------------
# Build blacklist
# -------------------------------
fixed_vars <- c("Race", "ABO_RH", "Visit_age", "sex")
all_vars <- names(bn_df)

blacklist <- expand.grid(from = all_vars, to = fixed_vars, stringsAsFactors = FALSE)
blacklist <- blacklist[blacklist$from != blacklist$to, ]

cat("Total blacklisted edges:", nrow(blacklist), "\n")
cat("Edges blocked into Race:", sum(blacklist$to == "Race"), "\n")
cat("Edges blocked into ABO_RH:", sum(blacklist$to == "ABO_RH"), "\n")
cat("Edges blocked into Visit_age:", sum(blacklist$to == "Visit_age"), "\n")
cat("Edges blocked into sex:", sum(blacklist$to == "sex"), "\n\n")

# -------------------------------
# Learn BN structure
# -------------------------------
cat("Learning BN structure (Hill-Climbing)...\n")
bn_model <- hc(bn_df, score = "pnal-cg", k = k_bic, blacklist = blacklist)
saveRDS(bn_model, "bn_model_pnal_bic_blacklist.rds")
cat("BN structure saved to 'bn_model_pnal_bic_blacklist.rds'\n\n")

# -------------------------------
# Fit parameters using hard-EM CG
# -------------------------------
cat("Fitting BN with hard-EM CG...\n")
n_cores <- 13
cl <- makeCluster(n_cores)
bn_fitted_hard_em <- bn.fit(bn_model, bn_df, method = "hard-em-cg", cluster = cl)
stopCluster(cl)
saveRDS(bn_fitted_hard_em, "bn_fitted_hard_em.rds")
cat("Hard-EM fitted BN saved to 'bn_fitted_hard_em.rds'\n\n")

# -------------------------------
# Visualize network
# -------------------------------
cat("Plotting network...\n")
g <- as.graphNEL(bn_model)
attrs <- list(
  node = list(
    fontsize = 20,
    shape = "ellipse",
    fillcolor = "lightblue",
    style = "filled"
  ),
  edge = list(color = "gray40")
)
# Save plot as PNG
png("bn_network.png", width = 1200, height = 800)
plot(g, attrs = attrs, "fdp")
dev.off()
cat("Network plot saved as 'bn_network.png'\n\n")

# -------------------------------
# Prediction using cpdist
# -------------------------------
cat("Generating predictions for HB_Value...\n")
evidence_list <- list(
  Race = factor("African Black", levels = levels(bn_df$Race)),
  Visit_age = 35,
  sex = factor("Male", levels = levels(bn_df$sex)),
  SITE_Type = factor("Mobile", levels = levels(bn_df$SITE_Type)),
  donor_type = factor("Repeat", levels = levels(bn_df$donor_type)),
  donation_time = 1245,
  DonProc = factor("WHOLE BLOOD", levels = levels(bn_df$DonProc))
)

n_samples <- 10000
samples <- cpdist(
  bn_fitted_hard_em,
  nodes = "HB_Value",
  evidence = evidence_list,
  method = "lw",
  n = n_samples
)

mean_hb <- mean(samples$HB_Value)
ci_95 <- quantile(samples$HB_Value, probs = c(0.025, 0.975))

cat("\n=== Prediction Results ===\n")
cat("Mean HB_Value:", round(mean_hb, 2), "\n")
cat("95% CrI:[", round(ci_95[1], 2), ",", round(ci_95[2], 2), "]\n")

# -------------------------------
# Save HTML report
# -------------------------------
cat("Saving HTML report...\n")
html_report <- tags$html(
  tags$head(tags$title("HC_SANBS Report")),
  tags$body(
    tags$h1("HC_SANBS Bayesian Network Report"),
    tags$h2("Network Plot"),
    tags$img(src = "bn_network.png", width = "800px"),
    tags$h2("Prediction Results"),
    tags$p(paste("Mean HB_Value:", round(mean_hb, 2))),
    tags$p(paste0("95% CrI: [", round(ci_95[1], 2), ", ", round(ci_95[2], 2), "]"))
  )
)

save_html(html_report, "HC_SANBS_report.html")
cat("HTML report saved to 'HC_SANBS_report.html'\n")
