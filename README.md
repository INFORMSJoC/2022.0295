
# 2022.0295

## Required Packages
To run this project, make sure you have the following packages installed. You can install them using:

```julia
using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Graphs")
Pkg.add("MetaGraphs")
Pkg.add("SimpleWeightedGraphs")
Pkg.add("StatsBase")
Pkg.add("Random")
Pkg.add("JuMP")
Pkg.add("GLPK")
Pkg.add("Distances")
Pkg.add("Tables")
Pkg.add("Statistics")
```

## Running the Project

### Step 1: Computing SAIDI and Reconnection Times
To perform all the computations for SAIDI and Reconnection times, run the following script:

```julia
include("Greensboro_example.jl")
```

> **Note:** This step must be completed before running the local search, as it populates the graphs required for further analysis.

### Step 2: Running Local Search
To generate all the local search data for the 18 components in the Greensboro dataset, run:

```julia
include("run_local_search.jl")
```

> **Note:** Ensure there is an `output` folder in the directory where the code is executed, as the results will be saved there.
