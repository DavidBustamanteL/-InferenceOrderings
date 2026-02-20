# Inference Orderings DMCP
## Replication file to Dirichlet Monte Carlo Perturbations - Barbaro & Bustamante Lazo

This repository implements the function ```run_one_election()``` to estimate Condorcet cycle probabilities in ranked-choice elections using Dirichlet Monte Carlo perturbations of empirically observed preference distributions.

The function evaluates all 3-candidate subsets **(triplets)** in an election and computes the probability of observing a majority cycle under sampling uncertainty.

This README focuses exclusively on the main function and its internal subfunctions.

# Main Function

## ```run_one_election(df, N, phi, seed, pseudocount, ord_levels)```

This function executes the complete Condorcet-cycle simulation pipeline for a single election.
