# Taxi-Customer Assignment using PSO

This repository contains the implementation of Particle Swarm Optimization (PSO) for assigning taxis to customers based on minimizing the Euclidean distance.

## Files

- `pso_taxi_customer.m` - Main script implementing the PSO algorithm.
- `EnsureOneToOne.m` - Function to ensure one-to-one assignment of taxis to customers.

## How to Run

1. Clone the repository.
2. Open MATLAB and navigate to the repository directory.
3. Run the script `pso_taxi_customer.m`.

## Parameters

- Number of particles: 30
- Maximum iterations: 100
- Inertia weight: 0.7
- Cognitive weight: 1.5
- Social weight: 1.5

## Results

The script will display the initial cost, minimum cost, and the best allocation of taxis to customers. It will also plot the convergence curve and visualize the results.

## Experiments

To experiment with different values of PSO parameters, modify the values in the script and rerun it.

## Future Work

- Extend the implementation to handle multi-objective optimization.
- Modify the algorithm to handle unequal numbers of taxis and customers.
- Introduce different categories of taxis and implement constraints.
