# persistence-detection

generateTopologies.py
Enumerates 2- and 3-node topologies, including AND and OR logic sampling at each node, and removes duplicate topologies.

runSimulations.py
Simulation manager that reads in a list of paramters, calls executable of simulate.cc for various input durations, and calls executable of doseResponseMetrics.cc to calculate metrics of kinetic filtering

simulate.cc
Given a set of parameters and input characteristics, simulates a circuit's response to a single pulse of input

doseResponseMetrics.cc
Given a set of inputs and responses, calculates steepness and EC50 of the dose response curve
