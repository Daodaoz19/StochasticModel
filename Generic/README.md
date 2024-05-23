# Generic Simulation for Reaction Diffusion Model

This directory includes files for a generic simulation of the Reaction Diffusion Model. Below is a detailed description of the included files, output files, and instructions on how to run and modify the simulation.

## Code Files

- **model.h, model.c**: The model class for the Reaction Diffusion model, default is Schnakenberg Model.
- **simulation.c**: The main program for the executable simulation file.
- **Makefile**: Make it and get the executable file.

## Output Files

- **firings**: The number of times each reaction fired during the simulation.
- **Location**: The population distribution of all species across the bins at the end of the simulation.
- **Population**: The final population distribution for each species across the bins at the end of the simulation.
- **Populationfile**: The population distribution for each species across the bins during each OutputPoints.
- **Total**: The final total population count for each species at the end of the simulation.

## Running the Simulation

To run the program, make and run the executable file:
./simulation [optional_input_initial_population_file]
If the `optional_input_initial_population_file` is not specified, the simulation will use the initial population from `initialValue` in `model.c`, where the population distributes equally over the bin.

## Modifying the Model

To change the model used in the simulation, modify the `initializeSimulation`, `prop3`, and `reactionChange` functions in `model.c`.

### Example Setup for the Schnakenberg Model

#### initializeSimulation Function

1. **Species Initialization**:
   - **Species X**:
     - Name: 'X'
     - Initial Value: 1
     - Diffusion Rate: 5e-4
   - **Species Y**:
     - Name: 'Y'
     - Initial Value: 2
     - Diffusion Rate: 5e-4 * 10000

2. **Reactions Initialization**:
   - **Reaction A -> X**:
     - Reaction Type: 0 (Zero Molecular Reaction)
     - Reaction Rate: 10.0
   - **Reaction X -> C**:
     - Reaction Type: 1 (One Molecular Reaction)
     - Reaction Rate: 2.0
     - Reactant Index: 0 (Index of Species X is 0)
   - **Reaction B -> Y**:
     - Reaction Type: 0 (Zero Molecular Reaction)
     - Reaction Rate: 30.0
   - **Reaction 2X + Y -> 3X**:
     - Reaction Type: 2 (Two or more Molecular Reaction)
     - Name of Function used to calculate propensity value: `prop3`

#### prop3 Function

The propensity value for the reaction 2X + Y -> 3X is calculated by:
(population of species X) * (population of species X - 1) * (population of species Y)

#### reactionChange Function

- **Reaction A -> X**:
  - Population of X in the bin that reaction fires increases
  - Total population of X increases
- **Reaction X -> C**:
  - Population of X in the bin that reaction fires decreases
  - Total population of X decreases
- **Reaction B -> Y**:
  - Population of Y in the bin that reaction fires increases
  - Total population of Y increases
- **Reaction 2X + Y -> 3X**:
  - Population of Y in the bin that reaction fires decreases
  - Population of X in the bin that reaction fires increases
  - Total population of X increases
  - Total population of Y decreases

## Mode

Modify the `#DEFINE Mode` in `model.h` to select from Auto and Manual:
- **Auto**: Automatically calculate the propensity value in simulation.
- **Manual**: Calculate the propensity for each reaction with a function in `model.c`.
To use Manual mode, also update the function `reaction_propensityChange` and `diffusionChange` in `model.c`
