library(reshape2)
library(mice)
set.seed(1000)

# 1) Data and variables
library(readxl)
myData_raw <- read_excel("T13-EQ.xlsx")
View(myData_raw)
# Variables by wave
vars_T1 <- c("T1fPEX", "T1fNEX", "T1mPEX", "T1mNEX", "T1SEA", "T1ROE", "T1UOE", "T1OEA", "T1DEA", "T1LPA", "T1SOM", "T1INT")
vars_T2 <- c("T2fPEX", "T2fNEX", "T2mPEX", "T2mNEX", "T2SEA", "T2ROE", "T2UOE", "T2OEA", "T2DEA", "T2LPA", "T2SOM", "T2INT")
vars_T3 <- c("T3fPEX", "T3fNEX", "T3mPEX", "T3mNEX", "T3SEA", "T3ROE", "T3UOE", "T3OEA", "T3DEA", "T3LPA", "T3SOM", "T3INT")
vars_all <- c(vars_T1, vars_T2, vars_T3)
# Covariates
covars <- c("gender", "age", "place_SES")
# Nodes
k <- length(vars_T1)     # 12 nodes
node <- c("fPEX", "fNEX", "mPEX", "mNEX", "SEA", "ROE", "UOE", "OEA", "DEA", "LPA", "SOM", "INT")
# Groups
groups <- c(rep("Parental emotional expressivity",4), rep("Emotional intelligence",4), rep("Depressive symptoms",4))

# 2) Analysis data
analysis_data <- myData_raw[, c(vars_all, covars)]
# Convert to numeric
for(j in seq_len(ncol(analysis_data))){
  if(!is.numeric(analysis_data[[j]])){
    analysis_data[[j]] <-
      as.numeric(analysis_data[[j]])
  }
}
summary(analysis_data)

# 3) Multiple imputation
# 20 PMM imputations
imp <- mice(
  analysis_data,
  m = 20,
  method = "pmm",
  maxit = 20,
  seed = 1000
)
print(imp)
md.pattern(analysis_data)

# 4) Correlation heatmap
#    T1-T2-T3
#    Excluding covariates
# Correlations: imputation 1
cor_data <- complete(imp, 1)
cor_data <- cor_data[, vars_all]
corMat <- cor(cor_data, use = "pairwise.complete.obs")
cor_df <- melt(corMat)
png(
  filename =
    "CorrelationMatrix_T1_T2_T3_MI.png",
  width = 15,
  height = 13,
  units = "in",
  res = 600
)
ggplot(
  cor_df,
  aes(x = Var1, y = Var2, fill = value)
) +
  geom_tile(color="white") +
  geom_text(aes(label=sprintf("%.2f", value)), size=2.6) +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0, limits=c(-1,1), name="r") +
  theme_minimal(base_size=11) +
  theme(
    axis.text.x =
      element_text(angle=45, vjust=1, hjust=1),
    axis.text.y =
      element_text(size=9),
    axis.title =
      element_blank(),
    panel.grid =
      element_blank()
  ) +
  coord_fixed() +
  # Wave separators
  geom_hline(yintercept = c(length(vars_T1)+0.5, length(vars_T1)+length(vars_T2)+0.5), linewidth=0.6) +
  geom_vline(xintercept = c(length(vars_T1)+0.5, length(vars_T1)+length(vars_T2)+0.5), linewidth=0.6) +
  labs(title= "Correlation Matrix (T1-T3)")
dev.off()

# Part 2/3
# CLPN estimation + Multiple Imputation pooling

# 5) CLPN estimator
estimate_CLPN <- function(
    data,
    Xvars,
    Yvars,
    covars,
    alpha = 1
){
  k_local <- length(Xvars)
  adjMat <- matrix(0, nrow = k_local, ncol = k_local)
  # Predictors: prior-wave nodes + covariates
  X_time <- data[,Xvars]
  X_cov <- data[,covars]
  X <- as.matrix(cbind(X_time, X_cov))
  # Regress each next-wave node
  for(i in 1:k_local){
    y <- as.matrix(data[,Yvars[i]])
    # Elastic Net
    cvfit <- cv.glmnet(
      X,
      y,
      alpha = alpha,
      family = "gaussian",
      nfolds = 10,
      standardize = TRUE
    )
    # Use lambda.1se
    lambda_use <- cvfit$lambda.1se
    coef_all <- as.numeric(coef(cvfit, s=lambda_use))
    # Remove intercept
    coef_all <- coef_all[-1]
    # Keep node effects only
    adjMat[,i] <-
      coef_all[1:k_local]
  }
  # Remove autoregressive edges
  diag(adjMat) <- 0
  rownames(adjMat) <- node
  colnames(adjMat) <- node
  return(adjMat)
}

# 6) Estimate CLPNs across 20 imputations
adj_T1_T2_list <- list()
adj_T2_T3_list <- list()
for(m in 1:20){
  cat("Running imputed dataset:", m, "\n")
  data_m <- complete(imp,m )
  # Standardize
  # Note:
  # Standardize covariates as well
  data_m <- as.data.frame(scale(data_m))
  # T1 → T2
  adj_T1_T2_list[[m]] <-
    estimate_CLPN(data = data_m, Xvars = vars_T1, Yvars = vars_T2, covars = covars, alpha = 1)
  # T2 → T3
  adj_T2_T3_list[[m]] <-
    estimate_CLPN(data = data_m, Xvars = vars_T2, Yvars = vars_T3, covars = covars, alpha = 1)
}

# 7) Pool 20 imputed networks
adj_T1_T2 <-
  Reduce("+", adj_T1_T2_list) / length(adj_T1_T2_list)
adj_T2_T3 <-
  Reduce("+", adj_T2_T3_list) / length(adj_T2_T3_list)

# 8) Final network matrices
print(round(adj_T1_T2, 3))
print(round(adj_T2_T3, 3))

# 9) Save network matrices
write.csv(adj_T1_T2, "CLPN_T1_T2_MI20_matrix.csv")
write.csv(adj_T2_T3, "CLPN_T2_T3_MI20_matrix.csv")
# Export directed edges in long format
# T1 -> T2 long format
df_T1_T2 <- melt(adj_T1_T2)
colnames(df_T1_T2) <- c("From_T1", "To_T2", "Weight")
df_T1_T2$Direction <- "T1_to_T2"
# T2 -> T3 long format
df_T2_T3 <- melt(adj_T2_T3)
colnames(df_T2_T3) <- c("From_T2", "To_T3", "Weight")
df_T2_T3$Direction <- "T2_to_T3"
# Optional: remove zero-weight edges
# df_T1_T2_nonzero <- df_T1_T2[df_T1_T2$Weight != 0, ]
# df_T2_T3_nonzero <- df_T2_T3[df_T2_T3$Weight != 0, ]
# Export each edge table
write.csv(df_T1_T2, "CLPN_T1_to_T2_edges_long.csv", row.names = FALSE)
write.csv(df_T2_T3, "CLPN_T2_to_T3_edges_long.csv", row.names = FALSE)
# Combine edge tables
# Standardized column names
colnames(df_T1_T2) <- c("Source", "Target", "Weight", "Direction")
colnames(df_T2_T3) <- c("Source", "Target", "Weight", "Direction")
df_all_edges <- rbind(df_T1_T2, df_T2_T3)
write.csv(df_all_edges, "CLPN_All_Directions_Edges_long.csv", row.names = FALSE)
cat("Directed edge tables exported.\n")

# Part 3/3
# CLPN Network Visualization

# 10) Fixed node layout
# Based on T1 -> T2
layout_fixed <- qgraph(
  adj_T1_T2,
  layout = "spring",
  DoNotPlot = TRUE
)$layout

# 11) Node colors
group_colors <- c(
  "Parental emotional expressivity" = "#33A02C", # green
  "Emotional intelligence" = "#FF7F00",          # orange
  "Depressive symptoms" = "#1F78B4"              # blue
)
node_colors <-
  group_colors[groups]

# 12) T1 -> T2 network
png(
  filename =
    "CLPN_T1_T2_MI20_fixedLayout.png",
  width = 12,
  height = 10,
  units = "in",
  res = 600
)
qgraph(
  adj_T1_T2,
  layout = layout_fixed,
  labels = node,
  groups = groups,
  color = node_colors,
  # Positive edges
  posCol = "blue",
  # Negative edges
  negCol = "red",
  vsize = 10,
  esize = 5,
  minimum = 0.03,
  fade = FALSE,
  legend = TRUE,
  legend.cex = 0.8,
  title =
    "CLPN: T1 → T2 (MI pooled 20)"
)
dev.off()

# 13) T2 -> T3 network
# Same layout
png(
  filename =
    "CLPN_T2_T3_MI20_fixedLayout.png",
  width = 12,
  height = 10,
  units = "in",
  res = 600
)
qgraph(
  adj_T2_T3,
  layout = layout_fixed,
  labels = node,
  groups = groups,
  color = node_colors,
  posCol = "blue",
  negCol = "red",
  vsize = 10,
  esize = 5,
  minimum = 0.03,
  fade = FALSE,
  legend = TRUE,
  legend.cex = 0.8,
  title =
    "CLPN: T2 → T3 (MI pooled 20)"
)
dev.off()

# 14) Combined network plot
# Left: T1 -> T2
# Right: T2 -> T3
# Shared legend
png(
  filename =
    "CLPN_MI20_combined_network.png",
  width = 24,
  height = 10,
  units = "in",
  res = 600
)
layout(
  matrix(c(1,2,3, 1,2,3), nrow = 2, ncol = 3, byrow = TRUE),
  widths =
    c(5,5,3),
  heights =
    c(1,1)
)
# Panel A
# T1 → T2
qgraph(
  adj_T1_T2,
  layout = layout_fixed,
  labels = node,
  groups = groups,
  color = node_colors,
  posCol = "blue",
  negCol = "red",
  vsize = 16,
  esize = 5,
  minimum = 0.03,
  fade = FALSE,
  legend = FALSE,
  title =
    "A. CLPN T1 → T2"
)
# Panel B
# T2 → T3
qgraph(
  adj_T2_T3,
  layout = layout_fixed,
  labels = node,
  groups = groups,
  color = node_colors,
  posCol = "blue",
  negCol = "red",
  vsize = 16,
  esize = 5,
  minimum = 0.03,
  fade = FALSE,
  legend = FALSE,
  title =
    "B. CLPN T2 → T3"
)
# Panel C
# Legend
par(mar=c(0,0,0,0))
plot.new()
legend(
  "center",
  legend =
    names(group_colors),
  fill =
    group_colors,
  horiz = FALSE,
  cex = 2.5,
  bty = "n",
  title =
    "Group"
)
dev.off()

# 15) Save network objects
save(adj_T1_T2, adj_T2_T3, layout_fixed, file = "CLPN_MI20_final_network_objects.RData")

# 16) Prepare strength objects
CLPN_list <- list("T1–T2" = list(adj_noAR = adj_T1_T2), "T2–T3" = list(adj_noAR = adj_T2_T3))
tp_names <- names(CLPN_list)

# 17) Compute Out-/In-strength
strength_colors <- c("Out-strength" = "#F8766D", "In-strength"  = "#00BFC4")
strength_all <- list()
for (tp in tp_names) {
  adj_noAR <- CLPN_list[[tp]]$adj_noAR
  # Out-strength: sum of absolute outgoing weights
  OutStrength <- rowSums(abs(adj_noAR), na.rm = TRUE)
  # In-strength: sum of absolute incoming weights
  InStrength <- colSums(abs(adj_noAR), na.rm = TRUE)
  strength_df <- data.frame(Node = node, OutStrength = OutStrength, InStrength = InStrength)
  strength_long <- melt(strength_df, id.vars = "Node", variable.name = "Type", value.name = "Strength")
  strength_long$Type <- factor(
    strength_long$Type,
    levels = c("OutStrength", "InStrength"),
    labels = c("Out-strength", "In-strength")
  )
  strength_long$Time <- tp
  strength_all[[tp]] <- strength_long
}
strength_all_df <- do.call(rbind, strength_all)

# 18) Export strength values
write.csv(strength_all_df, "CLPN_Out_In_Strength_values.csv", row.names = FALSE)
print(strength_all_df)

# 19) Plot Out-/In-strength
max_s <- max(strength_all_df$Strength, na.rm = TRUE)
pad <- max(0.05 * max_s, 0.1)
png(
  "Strength_T1T2_vs_T2T3.png",
  width = 14,
  height = 8,
  units = "in",
  res = 600
)
ggplot(
  strength_all_df,
  aes(x = Strength, y = reorder(Node, Strength), fill = Type)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", Strength)), position = position_dodge(width = 0.8), hjust = -0.15, size = 3) +
  facet_wrap(~ Time, nrow = 1) +
  scale_fill_manual(values = strength_colors) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, max_s + pad)) +
  labs(x = "Strength (absolute sum)", y = "Node", title = "Node Strength (Out vs. In)") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text = element_text(size = 10),
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )
dev.off()
# Prepare imputed data for bootnet
cat("Preparing imputed data...\n")
# Use imputed dataset 10
myData_df_imputed <- complete(imp, 10)
# Apply the same standardization
myData_df_imputed <- as.data.frame(scale(myData_df_imputed))
# Store data for downstream use
CLPN_list[["T1–T2"]]$data <- myData_df_imputed
CLPN_list[["T2–T3"]]$data <- myData_df_imputed

# 20) Bootnet accuracy and stability
library(bootnet)
library(glmnet)
NBOOT <- 1000
NCORE <- 1
# Use imputed dataset 15
boot_data <- complete(imp, 15)
# 20.1 Custom CLPN estimator
boot_CLPN <- function(data, Xvars, Yvars, covars, node_names){
  # Re-standardize each bootstrap sample
  data[] <- lapply(data, function(x){
    x <- as.numeric(x)
    if(sd(x) == 0){
      return(rep(0, length(x)))
    }
    as.numeric(scale(x))
  })
  X <- as.matrix(data[, c(Xvars, covars)])
  k <- length(Xvars)
  adj <- matrix(0, nrow = k, ncol = k, dimnames = list(node_names, node_names))
  for(i in seq_len(k)){
    y <- data[[Yvars[i]]]
    fit <- cv.glmnet(
      X,
      y,
      alpha = 1,
      family = "gaussian",
      nfolds = 10,
      standardize = FALSE
    )
    beta <- as.numeric(coef(fit, s = "lambda.1se"))[-1]
    adj[, i] <- beta[seq_len(k)]
  }
  # Remove autoregressive edges
  diag(adj) <- 0
  list(graph = adj)
}
# 20.2 Build bootnet objects
make_boot_network <- function(Xvars, Yvars){
  dat <- boot_data[, c(Xvars, Yvars, covars)]
  estimateNetwork(
    dat,
    default = "none",
    fun = boot_CLPN,
    labels = node,
    directed = TRUE,
    signed = TRUE,
    weighted = TRUE,
    Xvars = Xvars,
    Yvars = Yvars,
    covars = covars,
    node_names = node
  )
}
net_T1T2 <- make_boot_network(vars_T1, vars_T2)
net_T2T3 <- make_boot_network(vars_T2, vars_T3)
# 20.3 Bootstrap analysis
run_bootnet <- function(net, name){
  dir.create(paste0("Bootnet_", name), showWarnings = FALSE)
  outdir <- paste0("Bootnet_", name, "/")
  # 1. Nonparametric Bootstrap
  # Edge accuracy and strength stability
  set.seed(1000)
  boot_nonpara <- bootnet(
    net,
    nBoots = NBOOT,
    type = "nonparametric",
    nCores = NCORE,
    statistics = c("edge", "outStrength", "inStrength")
  )
  # Edge confidence intervals
  png(
    paste0(outdir, "01_Edge_CI.png"),
    width = 14,
    height = 11,
    units = "in",
    res = 600
  )
  print(plot(boot_nonpara, statistics = "edge", order = "sample", labels = TRUE))
  dev.off()
  # 2. Case-dropping Bootstrap
  # Out-/In-strength stability
  set.seed(1000)
  boot_case <- bootnet(
    net,
    nBoots = NBOOT,
    type = "case",
    nCores = NCORE,
    caseMin = 0.05,
    caseMax = 0.75,
    caseN = 10,
    statistics = c("outStrength", "inStrength")
  )
  # Case-dropping stability plot
  png(
    paste0(outdir, "05_Case_Dropping.png"),
    width = 12,
    height = 8,
    units = "in",
    res = 600
  )
  print(plot(boot_case, statistics = c("outStrength", "inStrength")))
  dev.off()
  # CS coefficients
  cs <- corStability(boot_case, statistics = c("outStrength", "inStrength"))
  write.csv(
    data.frame(Statistic = names(cs), CS = as.numeric(cs)),
    paste0(outdir, "CS_Coefficients.csv"),
    row.names = FALSE
  )
  # Save bootstrap objects
  saveRDS(boot_nonpara, paste0(outdir, "Bootstrap_Nonparametric.rds"))
  saveRDS(boot_case, paste0(outdir, "Bootstrap_CaseDropping.rds"))
  return(list(nonparametric = boot_nonpara, case = boot_case, CS = cs))
}
# Run both periods
result_T1T2 <- run_bootnet(
  net_T1T2,
  "T1_T2"
)
result_T2T3 <- run_bootnet(
  net_T2T3,
  "T2_T3"
)
# Inspect CS coefficients
result_T1T2$CS
result_T2T3$CS

# 21. CLPN network similarity
#    Four prespecified metrics:
#    (1) Edge-weight correlation
#    (2) Stable-edge overlap
#    (3) Cross-network Out-strength correlation
#    (4) Direction/sign-consistent Jaccard
SIM_DIR <- "CLPN_Network_Similarity"
dir.create(SIM_DIR, showWarnings = FALSE, recursive = TRUE)
# Stable edge: nonzero in >=50% of imputed networks
EDGE_INCLUSION_CUTOFF <- 0.50
# 21.1 Prepare pooled matrices
A <- as.matrix(adj_T1_T2)
B <- as.matrix(adj_T2_T3)
if(!all(dim(A) == dim(B))){
  stop("Network dimensions do not match.")
}
if(is.null(rownames(A))){
  rownames(A) <- node
  colnames(A) <- node
}
if(is.null(rownames(B))){
  rownames(B) <- node
  colnames(B) <- node
}
B <- B[rownames(A), colnames(A)]
diag(A) <- 0
diag(B) <- 0
node_names <- rownames(A)
# 21.2 Safe Spearman function
safe_spearman <- function(x, y, analysis_name){
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  if(length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2){
    return(data.frame(Analysis = analysis_name, N = length(x), Spearman_rho = NA_real_, P_value = NA_real_))
  }
  test_result <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(
    Analysis = analysis_name,
    N = length(x),
    Spearman_rho = unname(test_result$estimate),
    P_value = test_result$p.value
  )
}
# 21.3 Edge-weight correlation
off_diag <- row(A) != col(A)
edge_cor <- safe_spearman(A[off_diag], B[off_diag], "Edge-weight correlation")
# 21.4 Stable-edge selection
calculate_inclusion <- function(network_list, node_order){
  inclusion_list <- lapply(network_list, function(M){
    M <- as.matrix(M)
    if(!is.null(rownames(M))){
      M <- M[node_order, node_order]
    }
    diag(M) <- 0
    abs(M) > 1e-12
  })
  Reduce("+", inclusion_list) / length(inclusion_list)
}
if(exists("adj_T1_T2_list") && exists("adj_T2_T3_list")){
  inclusion_A <- calculate_inclusion(adj_T1_T2_list, node_names)
  inclusion_B <- calculate_inclusion(adj_T2_T3_list, node_names)
  present_A <- inclusion_A >= EDGE_INCLUSION_CUTOFF
  present_B <- inclusion_B >= EDGE_INCLUSION_CUTOFF
} else {
  stop("Imputed network lists not found.")
}
diag(present_A) <- FALSE
diag(present_B) <- FALSE
# 21.5 Stable-edge overlap
common_edges <- present_A & present_B
union_edges  <- present_A | present_B
n_A      <- sum(present_A)
n_B      <- sum(present_B)
n_common <- sum(common_edges)
n_union  <- sum(union_edges)
# Unsigned Jaccard for edge-set overlap
edge_overlap_jaccard <- if(n_union == 0) NA_real_ else n_common / n_union
# 21.6 Cross-network Out-strength correlation
#      Out-strength = sum of absolute outgoing edges
out_strength_A <- rowSums(abs(A), na.rm = TRUE)
out_strength_B <- rowSums(abs(B), na.rm = TRUE)
out_strength_cor <- safe_spearman(out_strength_A, out_strength_B, "Out-strength correlation")
out_strength_df <- data.frame(
  Node = node_names,
  OutStrength_T1_T2 = out_strength_A,
  OutStrength_T2_T3 = out_strength_B,
  stringsAsFactors = FALSE
)
# 21.7 Direction/sign-consistent Jaccard
#      Numerator: stable in both networks with same sign
#      Denominator: stable in either network
same_sign_edges <-
  common_edges &
  sign(A) != 0 &
  sign(B) != 0 &
  sign(A) == sign(B)
opposite_sign_edges <-
  common_edges &
  sign(A) != 0 &
  sign(B) != 0 &
  sign(A) != sign(B)
signed_jaccard <- if(n_union == 0){
  NA_real_
} else {
  sum(same_sign_edges) / n_union
}
# 21.8 Summary
correlation_results <- rbind(edge_cor, out_strength_cor)
similarity_results <- data.frame(
  Metric = c(
    "Stable edges: T1→T2",
    "Stable edges: T2→T3",
    "Common stable edges",
    "Stable edges in either network",
    "Stable-edge overlap Jaccard",
    "Common stable edges with same sign",
    "Common stable edges with opposite sign",
    "Direction/sign-consistent Jaccard"
  ),
  Value = c(
    n_A,
    n_B,
    n_common,
    n_union,
    edge_overlap_jaccard,
    sum(same_sign_edges),
    sum(opposite_sign_edges),
    signed_jaccard
  ),
  stringsAsFactors = FALSE
)
# 21.9 Print results
cat("\n========================================\n")
cat("CLPN network similarity\n")
cat("========================================\n\n")
cat("1. Edge-weight correlation:\n")
print(edge_cor, row.names = FALSE)
cat("\n2. Stable-edge overlap:\n")
cat("T1->T2 stable edges =", n_A, "\n")
cat("T2->T3 stable edges =", n_B, "\n")
cat("Common stable edges =", n_common, "\n")
cat("Stable-edge union =", n_union, "\n")
cat("Stable-edge Jaccard =", round(edge_overlap_jaccard, 3), "\n")
cat("\n3. Cross-network Out-strength correlation:\n")
print(out_strength_cor, row.names = FALSE)
cat("\n4. Direction/sign-consistent Jaccard:\n")
cat("Common same-sign edges =", sum(same_sign_edges), "\n")
cat("Common opposite-sign edges =", sum(opposite_sign_edges), "\n")
cat("Signed Jaccard =", round(signed_jaccard, 3), "\n")
# 21.10 Save results
write.csv(correlation_results, file.path(SIM_DIR, "01_Correlation_Results.csv"), row.names = FALSE)
write.csv(similarity_results, file.path(SIM_DIR, "02_Stable_Edge_and_Jaccard_Results.csv"), row.names = FALSE)
write.csv(out_strength_df, file.path(SIM_DIR, "03_OutStrength_Values.csv"), row.names = FALSE)
# Edge-level audit table
edge_index <- which(row(A) != col(A), arr.ind = TRUE)
edge_df <- data.frame(
  Source = rownames(A)[edge_index[,1]],
  Target = colnames(A)[edge_index[,2]],
  Weight_T1_T2 = A[edge_index],
  Weight_T2_T3 = B[edge_index],
  Inclusion_T1_T2 = inclusion_A[edge_index],
  Inclusion_T2_T3 = inclusion_B[edge_index],
  Stable_T1_T2 = present_A[edge_index],
  Stable_T2_T3 = present_B[edge_index],
  Common_Stable = common_edges[edge_index],
  Same_Sign = same_sign_edges[edge_index],
  Opposite_Sign = opposite_sign_edges[edge_index]
)
write.csv(edge_df, file.path(SIM_DIR, "04_All_Directed_Edges_Check.csv"), row.names = FALSE)
# 21.11 Edge-weight similarity plot
png(
  file.path(SIM_DIR, "EdgeWeight_Similarity.png"),
  width = 8,
  height = 7,
  units = "in",
  res = 600
)
plot(
  A[off_diag],
  B[off_diag],
  pch = 19,
  cex = 0.8,
  xlab = "Edge weight: T1 → T2",
  ylab = "Edge weight: T2 → T3",
  main = paste0(
    "Directed Edge-Weight Similarity\n",
    "Spearman ρ = ", round(edge_cor$Spearman_rho, 3),
    ", p = ", format.pval(edge_cor$P_value, digits = 3, eps = 0.001)
  )
)
abline(h = 0, v = 0, lty = 3)
abline(a = 0, b = 1, lty = 2)
dev.off()
save(
  edge_df,
  out_strength_df,
  correlation_results,
  similarity_results,
  signed_jaccard,
  edge_overlap_jaccard,
  file = file.path(SIM_DIR, "CLPN_Network_Similarity_Objects.RData")
)
cat("\nSimilarity analysis saved to:\n")
cat(normalizePath(SIM_DIR), "\n")
