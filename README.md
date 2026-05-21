# IN104: Monte Carlo Simulation of the Ising Model

## How to Use This Project

This project simulates the 2D Ising model using two Monte Carlo algorithms:

- **Metropolis** (local updates)
- **Swendsen–Wang** (cluster updates)

The user can switch between algorithms, configure physical parameters, run the simulation, and generate plots of the main observables.

---

### 1. Compile the Project

From the root directory, run:

	make all

### 2. Configure the simulation

All parameters are defined in: 

	config.txt

- T — initial temperature
- B — external magnetic field
- Lx, Ly — lattice dimensions
- N_tot — number of Monte Carlo steps for measurements
- DeltaT, dT — temperature sweep range and step
- N_eq_metro — equilibration steps for Metropolis
- N_eq_sw — equilibration steps for Swendsen–Wang
- SwendsenWang — algorithm selector

### 3. Run the Simulation

Once the configuration file is set, run:

	./main

The program will:
- read config.txt
- run the simulation using the selected algorithm
- compute *observables* such as:
	- magnetization $M$
	- energy $E$
	- heat capacity $C_v$
	- binder cumulant $b$
	- autocorrelation function $C_{m^2}$

save all results in the **values/** directory

### 4. Generate Plots
To visualize the results, run:

```bash
python3 plot.py
```

This script:
1. loads the data from **values/**
2. generates plots for each observable
3. saves the figures in the **plots/** directory

### 5. Summary

1. make all
2. Edit config.txt
3. ./main
4. python3 plot.py

## Authors 
Maximilien Schirm, Baptiste Respaud
