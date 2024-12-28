[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

<!-- # 2022.0295 -->

# Fair and Reliable Reconnections for Temporary Disruptions in Electric Distribution Networks

This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported in the paper [Fair and Reliable Reconnections for Temporary Disruptions in Electric Distribution Networks](href) by Swati Gupta, Cyrus Hettle, and Daniel Molzahn.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2022.0295

https://doi.org/10.1287/ijoc.2022.0295.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{GHM2024,
  author       = {Swati Gupta and
                  Cyrus Hettle and
                  Daniel Molzahn},
  publisher =     {INFORMS Journal on Computing},
  title =         {Fair and Reliable Reconnections for Temporary Disruptions in Electric Distribution Networks},
  year =          {2024},
  doi =           {10.1287/ijoc.2022.0295.cd},
  note =          {Available for download at: https://github.com/INFORMSJoC/2022.0295},
}  
```

## Required Packages
We use julia version 1.11.0 for the experiments. To run this project, make sure you have the following packages installed. You can install them using:

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
## Data

We test our methods for improving outage metrics on the Greensboro, NC urban-suburban synthetic network from the NREL SMART-DS project. The relevant data from NREL SMART-DS dataset in an accessible CSV format in the data folder.

## Running the Project

### Step 1: Computing SAIDI and Reconnection Times
To perform all the computations for SAIDI and Reconnection times, run the following script (in the scripts folder):

```julia
include("Greensboro_example.jl")
```

> **Note:** This step must be completed before running the local search, as it populates the graphs required for further analysis.

### Step 2: Running Local Search
To generate all the local search data for the 18 components in the Greensboro dataset, run the following script from the scripts folder:

```julia
include("run_local_search.jl")
```

> **Note:** Ensure there is an `output` folder in the directory where the code is executed, as the results will be saved there.
