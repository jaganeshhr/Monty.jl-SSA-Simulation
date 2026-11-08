This repository contains the code used to estimate and compare the carbon dioxide removal potential of basaltic rock samples under enhanced rock weathering. 
It has two parts:

Monty_jl_SSA_Simulation.ipynb (Julia): simulates a field deployment of one feedstock using Monty.jl and exports per-iteration results as CSV files.

SSA CDR Calculator (Python): reads those CSV files for every sample, computes cumulative CDR, and ranks the samples.

Requirements:
Julia 1.10 with a Jupyter Julia kernel (IJulia)
Packages: Monty.jl (from GitHub), Random, Distributions, Meshes, GeoStatsFunctions, GeoStatsProcesses, DataFrames, Statistics, Plots, DimensionalData, CairoMakie, CSV
