###################################### All Elections #####################################################
##########################################################################################################


# WD and clean memory
rm(list = ls())
# Windows -> DO NOT USE mclapply()
# Linux -> You can run the script normally
setwd("")

# Avoiding scientific notation
options(scipen = 999)


# Packages
library(dplyr)
library(magrittr)
library(tidyverse)
library(tidyr)
library(stringr)
library(vote)
library(data.table)
library(MCMCpack)
library(parallel)
 

#### 1. Initial Setup ####

# Nr. of Reps
nb.rep = 10000

# Getting the data 
load("rankings.RData")

# ranking.list should be a list of 226 elections
elections_raw = ranking.list
length(elections_raw)

# excluding Belgium (4 elections)
elections_raw %<>% `[`(-c(12:15))
length(elections_raw)

# Political Weights
load("resW.RData")
Peru_2011 = res.partyW$Peru_2011
Peru_2011$totals

# replacing Peru_2011 w. pol. weights into Election Rankings(elections_raw)
elections_raw[["Peru_2011"]] = as.data.frame(Peru_2011$data)

# Ensure each element is a data.frame (robust coercion)
elections = lapply(elections_raw, function(x) {
  if (is.data.frame(x)) return(x)
  if (is.matrix(x))     return(as.data.frame(x))
  if (is.list(x))       return(as.data.frame(x))
  stop("Element not convertible to data.frame. Class: ", paste(class(x), collapse=", "))
})

# Optional: name them with ids (if lengths match)
if (exists("caseid.party") && length(caseid.party) == length(elections)) {
  names(elections) = caseid.party
}


#### 2. Full pipeline as a function ####
run_one_election = function(
  df,
  N = nb.rep,
  phi = 1,
  seed = 55234,
  pseudocount = 0,
  ord_levels = c(
    "A=B=C",
    "B>C>A",
    "A>C>B",
    "C>B>A",
    "C>A>B",
    "B=C>A",
    "C>A=B",
    "A>B=C",
    "A=C>B",
    "B>A=C",
    "B>A>C",
    "A>B>C",
    "A=B>C"
  )
) {
  df = as.data.frame(df)
  
  # Setting as NAs elections with less than 3 candidates
  if (ncol(df) < 3) {
    return(list(
      overall_cycle_rate = NA_real_,
      cycle_rate_by_triplet = numeric(0),
      baseline_outcome = character(0),
      note = "fewer than 3 candidates"
    ))
  }

  # (optional) if ranks are not numeric yet, coerce:
  df[] = lapply(df, function(col) suppressWarnings(as.numeric(col)))

  # 2.1 Triplets + encoding #
  triplets = combn(colnames(df), 3, simplify = FALSE)

  one_triplet_result = function(triplet) {

    df_three = df %>%
      dplyr::select(all_of(triplet)) %>%
      setNames(LETTERS[1:3])

    encode_ranking = function(x) {
      split(names(x), x) |>
        lapply(sort) |>
        sapply(paste, collapse = "=") |>
        paste(collapse = ">")
    }

    ranking_strings = apply(
      X = df_three,
      MARGIN = 1,
      FUN = encode_ranking
    )

    ranking_table = table(ranking_strings)

    ranking_df = as.data.frame(ranking_table) %>%
      dplyr::rename(ordering = ranking_strings,
                    frequency = Freq) %>%
      dplyr::arrange(dplyr::desc(frequency)) %>%
      dplyr::mutate(triplet = paste(triplet, collapse = ","), .before = 1)

    ranking_df
  }

  all_triplet_dfs = lapply(triplets, one_triplet_result)
  ranking_all = dplyr::bind_rows(all_triplet_dfs)

  # 2.2 Dirichlet Preparation #
  make_vec13 = function(ranking_one_triplet_df, ord_levels) {
    out = numeric(length(ord_levels))
    for (i in seq_along(ord_levels)) {
      ord = ord_levels[i]
      hit = ranking_one_triplet_df$frequency[ranking_one_triplet_df$ordering == ord]
      out[i] = if (length(hit) == 0) 0 else hit[1]
    }
    out
  }

  vec_list = split(ranking_all, ranking_all$triplet) |>
    lapply(make_vec13, ord_levels = ord_levels)

  prop_list = lapply(vec_list, function(v) v / sum(v))

  # 2.3 Dirichlet Perturbation #
  alpha = lapply(vec_list, function(v) v * phi)

  set.seed(seed)

  perturbations = lapply(alpha, function(a) {
    # Dirichlet needs alpha > 0; if there are zeros, uncomment one of:   ---> not the case for us
    # a = pmax(a, 1e-8)         # tiny epsilon
    # a = a + 1                 # Laplace pseudocount
    MCMCpack::rdirichlet(N, a)
  })

  # Preserving the labels

  perturbations = lapply(perturbations, function(m) {
    if (is.null(colnames(m))) {
      if (ncol(m) != length(ord_levels)) stop("Number of columns does not match ord_levels.")
      colnames(m) = ord_levels
    }
    m
  })

  names(perturbations) = names(alpha)


  # 2.4 Cycles #
  pair_share = function(ord, x, y) {
    groups = strsplit(ord, ">", fixed = TRUE)[[1]]
    ranks = list()
    for (i in seq_along(groups)) {
      items = strsplit(groups[i], "=", fixed = TRUE)[[1]]
      for (it in items) ranks[[it]] = i
    }

    rx = ranks[[x]]
    ry = ranks[[y]]
    if (is.null(rx) || is.null(ry)) stop("Missing candidate in ordering: ", ord)

    if (rx < ry) return(1)
    if (rx > ry) return(0)
    0.5
  }

  AB_share = sapply(ord_levels, \(o) pair_share(o, "A","B"))
  AC_share = sapply(ord_levels, \(o) pair_share(o, "A","C"))
  BC_share = sapply(ord_levels, \(o) pair_share(o, "B","C"))

  idx = list(
    APB = which(AB_share == 1),
    BPA = which(AB_share == 0),
    ABt = which(AB_share == 0.5),

    APC = which(AC_share == 1),
    CPA = which(AC_share == 0),
    ACt = which(AC_share == 0.5),

    BPC = which(BC_share == 1),
    CPB = which(BC_share == 0),
    BCt = which(BC_share == 0.5)
  )

  condorcet_outcomes_matrix = function(M, idx, ord_levels = NULL) {

    if (!is.null(ord_levels)) {
      M = M[, ord_levels, drop = FALSE]
    }

    APB = rowSums(M[, idx$APB, drop = FALSE]) + 0.5 * rowSums(M[, idx$ABt, drop = FALSE])
    BPA = rowSums(M[, idx$BPA, drop = FALSE]) + 0.5 * rowSums(M[, idx$ABt, drop = FALSE])

    APC = rowSums(M[, idx$APC, drop = FALSE]) + 0.5 * rowSums(M[, idx$ACt, drop = FALSE])
    CPA = rowSums(M[, idx$CPA, drop = FALSE]) + 0.5 * rowSums(M[, idx$ACt, drop = FALSE])

    BPC = rowSums(M[, idx$BPC, drop = FALSE]) + 0.5 * rowSums(M[, idx$BCt, drop = FALSE])
    CPB = rowSums(M[, idx$CPB, drop = FALSE]) + 0.5 * rowSums(M[, idx$BCt, drop = FALSE])

    ab = ifelse(APB > BPA, "A", ifelse(BPA > APB, "B", "T"))
    ac = ifelse(APC > CPA, "A", ifelse(CPA > APC, "C", "T"))
    bc = ifelse(BPC > CPB, "B", ifelse(CPB > BPC, "C", "T"))

    out = character(nrow(M))

    out[ab=="A" & ac=="A"] = "A"
    out[ab=="B" & bc=="B"] = "B"
    out[ac=="C" & bc=="C"] = "C"

    out[out=="" & ab=="T" & ac=="A" & bc=="B"] = "Tie"
    out[out=="" & ac=="T" & ab=="A" & bc=="C"] = "Tie"
    out[out=="" & bc=="T" & ab=="B" & ac=="C"] = "Tie"

    out[out=="" & ab=="T" & ac=="T" & bc=="T"] = "Tie"
    out[out=="" & (ab=="T" | ac=="T" | bc=="T")] = "Tie"

    out[out==""] = "Cycle"
    out
  }

  outcomes_by_triplet = lapply(perturbations, condorcet_outcomes_matrix, idx = idx, ord_levels = ord_levels)

  overall_table = table(unlist(outcomes_by_triplet, use.names = FALSE))
  cycle_rate_by_triplet = sapply(outcomes_by_triplet, \(x) mean(x == "Cycle"))
  overall_cycle_rate = mean(unlist(outcomes_by_triplet, use.names = FALSE) == "Cycle")

  baseline_outcome = sapply(prop_list, function(p) {
    M0 = matrix(p, nrow = 1)
    colnames(M0) = ord_levels
    condorcet_outcomes_matrix(M0, idx = idx, ord_levels = ord_levels)
  })

  # IMPORTANT: explicit return
  return(list(
    ranking_all = ranking_all,
    vec_list = vec_list,
    prop_list = prop_list,
    perturbations = perturbations,
    outcomes_by_triplet = outcomes_by_triplet,
    overall_table = overall_table,
    cycle_rate_by_triplet = cycle_rate_by_triplet,
    overall_cycle_rate = overall_cycle_rate,
    baseline_outcome = baseline_outcome
  ))

  # deleting not longer needed objects from the RAM!!!
  gc()
}

# total nr of ind-level obs -> 381,770
n_obs_per_election = sapply(elections_raw, nrow)
sum(n_obs_per_election)

# average nr of individuals -> 1,719
n_respondents = sapply(elections_raw, nrow)
mean(n_respondents)

# average nr of candidates/parties -> 6.42
n_parties = sapply(elections_raw, ncol)
mean(n_parties)


# total nr of triplets -> 7,270
triplets_in_election = function(elections_raw) {
  m = ncol(elections_raw)
  if (m < 3) return(0)
  choose(m, 3)
}
triplets_by_election = vapply(elections_raw, triplets_in_election, numeric(1))
total_triplets = sum(triplets_by_election)
total_triplets

# Candidates per elections
candidates_per_election = sapply(elections_raw, ncol)

candidate_groups = cut(
  candidates_per_election,
  breaks = c(-Inf, 5, 6, 7, Inf),
  labels = c("≤5", "6", "7", "8+"),
  right = TRUE
)

table(candidate_groups)

# table data
cand_dist = table(n_parties)

summary_table = data.frame(
  metric = c(
    "Total individual observations",
    "Average respondents per election",
    "Average candidates per election",
    "Total number of triplets",
    "Elections with <5 candidates",
    "Elections with 6 candidates",
    "Elections with 7 candidates",
    "Elections with 8+ candidates"
  ),
  value = c(
    sum(n_obs_per_election),
    mean(n_obs_per_election),
    mean(n_parties),
    sum(triplets_by_election),
    sum(cand_dist[as.numeric(names(cand_dist)) <= 5], na.rm = TRUE),
    unname(cand_dist["6"]),
    unname(cand_dist["7"]),
    sum(cand_dist[as.numeric(names(cand_dist)) >= 8], na.rm = TRUE)
  )
)

summary_table




# Results ---> WITH LINUX!!! 
# Urgent: Adjust the function in case you run this script on Windows; in particular, don't use mclapply!
# NO parallelization -> around 13:5 min
timing1 = system.time({
results_nocores = lapply(
  elections,
  run_one_election,
  phi = 1, seed = 55234
)
})
timing1
as.numeric(timing1["elapsed"]) / 60

cycle_rates1 = vapply(results_nocores, \(x) x$overall_cycle_rate, numeric(1))
print(summary(cycle_rates1))


# Parallelization  -> slightly over 5 min
timing2 = system.time({
results_cores = mclapply(
  elections,
  run_one_election,
  phi = 1, seed = 55234,
  mc.cores = parallel::detectCores(),
  mc.set.seed = TRUE
)
})
timing2
print(as.numeric(timing2["elapsed"]) / 60)

cycle_rates2 = vapply(results_cores, \(x) x$overall_cycle_rate, numeric(1))
print(summary(cycle_rates2))

# how often do the cycles appear? -> new function with parall. results   --> Peru_2011 out once the pol. weights are added
# Where are the "fragile" elections
fragile = which(cycle_rates2 > 0.01)                      # >1% cycle probability (non w. >5%)
#fragile
print(cycle_rates2[fragile])

# and Switzerland_2011? (or any specific election)
#cycle_rates2["Switzerland_2011"]                          # Overall not bad, just some triplet trouble some (see case study, same results)



#### 5. By Triplets #####

thresholds = c(.01, .02, .03, .04, .05)

# 1) Extracting triplet-level cycle probs from each election result
triplets_by_election = lapply(results_cores, \(x) x$cycle_rate)

# 2) Flatten across elections -> one long vector of triplet P(cycle)
p_triplet = as.numeric(unlist(triplets_by_election, use.names = FALSE))
p_triplet = p_triplet[!is.na(p_triplet)]

# checks
max(p_triplet, na.rm = TRUE)
sum(p_triplet >= 0.05, na.rm = TRUE)
summary(p_triplet)
p_triplet

# 3) Threshold counts/shares
df_thresh = tibble(p = p_triplet) %>%
  crossing(threshold = thresholds) %>%
  group_by(threshold) %>%
  summarise(
    count = sum(p >= threshold),
    share = mean(p >= threshold),
    .groups = "drop"
  ) %>%
  mutate(threshold_pct = percent(threshold, accuracy = 1))

# 4) Histogram distribution of triplet cycle probabilities
ggplot(tibble(p = p_triplet), aes(x = p)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = thresholds, linetype = "dashed") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "P(cycle) per triplet",
    y = "Count (triplet × election)",
    title = "Distribution of triplet cycle probabilities"
  ) +
  theme_minimal()

# 5) Bar plot: how many triplets exceed 1–5%
ggplot(df_thresh, aes(x = threshold_pct, y = count)) +
  geom_col() +
  labs(
    x = "Cutoff for P(cycle)",
    y = "Count of triplet observations ≥ cutoff",
    title = "Triplets exceeding cycle-probability thresholds"
  ) +
  theme_minimal()
