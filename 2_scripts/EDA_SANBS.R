library(pacman)
p_load(tidyverse, bnlearn, dtplyr, Rgraphviz)

# Load and prepare data
sanbs_df <- read.csv("/d3mod/donor_data/SANBS_2019-2022_2025May.csv")

sanbs_dt <- lazy_dt(sanbs_df)
sanbs_dt <- sanbs_dt %>%
  mutate(Visit_Date = as.Date(Visit_Date),
         next_visit = as.Date(next_visit),
         Def_start = as.Date(Def_start),
         Def_end = as.Date(Def_end))

# Prepare data for Bayesian network
bn_df <- sanbs_dt %>%
  select(-Visit_Date, -next_visit, -Def_start, -Def_end, 
         -DonorID, -MobileID, -Visit_yr, -Visit_Mo, -EffectiveYr, -EffectiveMo,
         -donation_product, -DEF_LENGTH, -DefCode, -don_seq, -col_Int) %>%
  filter(Deferral_permanent != "Perm") %>%
  as.data.frame() %>%
  mutate(across(where(is.character), as.factor)) %>%
  mutate(across(where(is.integer), as.numeric))

# Learn Bayesian networks
bn_model <- hc(bn_df, score = "nal-cg")
bn_model_sparse_10 <- hc(bn_df, score = "pnal-cg", k = 10)

# Save models
saveRDS(bn_model, "bn_model_nal.rds")
saveRDS(bn_model_sparse_10, "bn_model_pnal_10.rds")

# Plot models
g_nal <- as.graphNEL(bn_model)
g_sparse <- as.graphNEL(bn_model_sparse_10)

attrs <- list(
  node = list(
    fontsize = 30,           
    shape = "ellipse",
    fillcolor = "lightblue",
    style = "filled"
  ),
  edge = list(
    color = "gray40"
  )
)

# Plot NAL model
plot(g_nal, attrs = attrs, "fdp")

# Plot PNAL k=10 model
plot(g_sparse, attrs = attrs, "fdp")