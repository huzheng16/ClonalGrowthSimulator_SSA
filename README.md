# Clonal Growth Simulator (SSA)

- [Installation](#installation)
- [Usage](#usage)
  - [Run](#run)
  - [Configuration files](#configuration-files)
- [Docker](#docker)
- [Acknowledgements](#acknowledgements)

## Installation

Dependencies:
- A C++ compiler that supports C++11
- Boost libraries: system, filesystem, iostreams, random
- CMake 2.6 or newer
- Gzip (optional)

On Ubuntu 16.04 all dependencies can be installed with the following command:

```
sudo apt-get install cmake libboost-filesystem-dev libboost-iostreams-dev libboost-random-dev gcc g++ python gzip```

Build:
```bash
mkdir build
cd build
cmake ../
make
```


## Usage

### Run
To run a simulation you need to call the binary with an xml file:

```bash
cd examples
mkdir out
../bin/simulator iteratedgrowth_csc_otl.xml
```

Note that the code will crash if the output folder (`<OutPath>out/</OutPath>`) does not exist!

### Configuration files
The simulation settings are stored in a single xml file containing 3 nodes: `Simulation`, `Solver`, and `Experiment`. An example xml can be found in the `examples` folder.


#### Model
The model node specifies the growth model, which includes the species and transitions, and the related settings. The growth model is specified via the `type` attribute, and may be:

- **Simple**: constant growth with one cell type without cell death
- **CSC**: stem cell driven growth with a specified number of division for differentiated cells, and all cells may go into apoptosis.


##### Parameters for the simple growth model

| Element        | Description                                 |
|----------------|---------------------------------------------|
| DivisionRate   | Number of times cells divide per day        |
| DivisionRateSD | Standard deviation of the division rates    |
| DivisionRateLogNormal | Use lognormal instead of normal distribution for division rates |
| DeathRate      | Death rate of the cells                     |
| Seed           | Seed for generating variable division rates |


##### Parameters for the CSC growth model

| Element                         | Description                                                                 |
|---------------------------------|-----------------------------------------------------------------------------|
| MaxAgeTransitCell               | Number of transitions for DCs                                               |
| DivisionRateStemCell            | Number of times CSCs divide per day                                         |
| DivisionRateTransitCell         | Number of times DCs divide per day                                          |
| DeathRateStemCell               | Number of times CSCs die per day                                            |
| DeathRateTransitCell            | Number of times DCs die per day                                             |
| ProbSymmStemCellDivision        | Probability of a symmetric CSC division                                     |
| ProbSymmStemCellDifferentiation | Probability of a symmetric CSC differentiation                              |
| DivisionRateSD                  | Standard deviation of the division rates                                    |
| Seed                            | Seed for generating variable division rates                                 |
| InitialFractions                | Comma-separated list of MaxAgeTransitCell+1 initial fractions               |
| InitialStemCellFraction         | Initial fraction of CSCs, all other species get identical initial fractions |


#### Solver

The solver node specifies the solver and related settings. The solver is set via the `type` attribute, which, for now, can only be **GillspieOTL**, a solver using ordinary tau leaping gillespie. The only setting for **GillespieOTL** solver is **tau**, which refers to the length (in days) of one tau-step.


#### Experiment

The experiment node specifies the kind of experiments and the experiment specific settings. The experiment is specified via the `type` attribute, and can be:

- **ConstantGrowth**: cells grow for a given time
- **IteratedGrowth**: cells grow untill a target number is reached and then a subset is passed to the next generation.

##### Parameters for ConstantGrowth

| Element               | Description                                                                    |
|-----------------------|--------------------------------------------------------------------------------|
| SimulationTime        | Simulation time in days                                                        |
| Seed                  | Seed used for stochastic growth and passage                                    |
| InitSeed              | Seed used for the initialization (if ommitted, Seed is used)                   |
| InitialPopulationSize | Initial number of cells                                                        |
| InitFile              | File with initial clone distribution                                           |
| InitUniform           | Number of clones over which the initial cells are uniformly distributed        |
| Name                  | Simulation name, used for filenames of the output                              |
| SaveFreq              | Frequency for generating output                                                |
| Outpath               | Path to store output to                                                        |
| SaveXML               | Export all settings (including defaults) in an xml file and save it to Outpath |
| gzip                  | Gzip all output files                                                          |
| UseSimDir             | Create a folder in Outpath that contains all generated output                  |


##### Parameters for IteratedGrowth

| Element                  | Description                                                                    |
|--------------------------|--------------------------------------------------------------------------------|
| CriticalPopulationSize   | Population size after which a growth step is stopped                           |
| NumberOfCellsToKeep      | Number of cells to keep after passage                                          |
| NumberOfPassages         | Number of growth and passage cycles                                            |
| MaxPassTime              | Maximum time a growth step may take                                            |
| StopSimIfNPassNotReached | Stop the simulation if MaxPassTime is reached                                  |
| Seed                     | Seed used for stochastic growth and passage                                    |
| InitSeed                 | Seed used for the initialization (if ommitted, Seed is used)                   |
| InitialPopulationSize    | Initial number of cells                                                        |
| InitFile                 | File with initial clone distribution                                           |
| InitUniform              | Number of clones over which the initial cells are uniformly distributed        |
| Name                     | Simulation name, used for filenames of the output                              |
| SaveFreq                 | Frequency for generating output                                                |
| Outpath                  | Path to store output to                                                        |
| SaveXML                  | Export all settings (including defaults) in an xml file and save it to Outpath |
| gzip                     | Gzip all output files                                                          |
| UseSimDir                | Create a folder in Outpath that contains all generated output                  |


## Docker
Alternatively, you can use Docker to build the software:

```docker build -t ssa .```

Then, you can run the software with:

```
cd examples/
mkdir out
docker run -u $(id -u):$(id -g) -v $PWD:/data/ ssa /data/iteratedgrowth_csc_otl_docker.xml
```

Note that this command mounts the current directory as `/data/` inside the container. Therefore, `OutPath` in the xml must start with `/data/` as well.


## Acknowledgements

We thank the developers of the TinyXML-2 library for providing this library (https://github.com/leethomason/tinyxml2).
