# Main parameters

## Genetic algorithm parameters

### Population

- POPULATION_LENGTH:
	- The size of the genetic algorithm population. 
	- Default: `100`. 
	- Range: `>0 [int]`
- POPULATION_ORIGIN
	- determines whether initial population is read from file or initizalized randomly. 
	- Default: `random`.
	- Options: `random`, `file`
- POPULATION_FILL_TYPE
	- if `POPULATION_ORIGIN` is file, this determines whether population is uniformly seeded with read organims, or complemented by random organisms.
	- Default: `random`.
	- Options: `random`, `uniform`
- RECOMBINATION_PROBABILITY :0.5
	- The probability that two parents will be recombined to produce offspring

### Stopping criteria

- END_WHILE_METHOD
	- Termination criteria. Determines whether the algorithm terminates after a fixed number of `iterations`, when the `fitness` meets some upper bound or when the difference in fitness between two consecutive iterations is lower than some `threshold`
	- Default: `iterations`.
	- Options: `iterations`, `fitness`, `threshold`
- MIN_ITERATIONS
	- Number of iterations to compute before stopping (by iterations)
	- Default: `100`
	- Range: `>0 [int]`
- MIN_FITNESS
	- Fitness to be achieved before stopping (by fitness)
	- Default: `10`
	- Range: `>0 [int]`
- THRESHOLD
	- Difference between consecutive iteration fitnesses to be  achieved before stopping (by threshold)
	- Default: `0.05`
	- Range: `>0 [float]`

### Fitness
- MAX_SEQUENCES_TO_FIT_POS
	- Number of sequences used to assess the *positive* dataset component of fitness
	- Default: `50`
	- Range: `>0 [int]`
- MAX_SEQUENCES_TO_FIT_NEG
	- Number of sequences used to assess the *negative* dataset component of fitness
	- Default: `50`
	- Range: `>0 [int]`
- COMPLEXITY_FACTOR
	- Strength of complexity penalty (multiplicative on the complexity term)
	- Default: `0.5`
	- Range: `>0 [float]`

### Directory organization

- DATASET_BASE_PATH_DIR
	- Base directory for datasets.
	- Default: `datasets/`
- RESULT_BASE_PATH_DIR
	- Base directory for results [where best organisms are reported].
	- Default: `results/`
- RESULT_TEST_BASE_PATH_DIR
	- Base directory for test results [where test organisms are assessed].
	- Default: `resultsTEST/`
- POSITIVE_FILENAME
  - Name of the positive dataset
- NEGATIVE_FILENAME
  - Name of the control dataset
- INPUT_FILENAME
  - Used in both search and test organisms. specifies the path to the file that includes the organisms that will be used as input in the program.
- OUTPUT_FILENAME
  - Name for the file where the output of the program will be redirected.

### Others

- PERIODIC_EXPORT
  - Number of iterations used to export periodically the max scored organism. Sometimes it's hard to track how the GA is doing, so every X iterations it exports the max organism on that iteration.

## Organism

- CUMULATIVE_FIT_METHOD
  - This is how we compute the fitness of a dataset when we have the energy on individual sequences.
  - Options: sum, mean
- MUTATE_PROBABILITY_SUBSTITUTE_PSSM
  - Probability to execute the Substitute pssm mutator. It creates a random pssm recognizer and inserts it in a random node position.
  - Range
- MUTATE_PROBABILITY_RISE_CHILD
  - Probability to execute the rise child mutator. It selects a random node in the organism and 
- MUTATE_PROBABILITY_NODE_MUTATION
- MUTATE_PROBABILITY_SUNK_CHILD
- MIN_NODES
- MAX_NODES

## Factory

- INITIAL_CONNECTOR_PROBABILITY
- REDUCER_PROBABILITY_FACTOR
- MIN_MU
- MAX_MU
- MIN_SIGMA
- MAX_SIGMA
- PWM_LENGTH
- PWM_PROBABILITY_STEP
- PWM_PROBABILITY_BASE
- PWM_PROBABILITY_DECIMALS

## Connector

- MUTATE_PROBABILITY_SIGMA
- MUTATE_PROBABILITY_MU
- MUTATE_PROBABILITY_SWAP
- TAU
- MUTATE_VARIANCE_SIGMA
- MUTATE_VARIANCE_MU
- PLACEMENT_OPTIONS

## PSSM recognizer

- MUTATE_PROBABILITY_RANDOM_COL
- MUTATE_PROBABILITY_FLIP_COL
- MUTATE_PROBABILITY_FLIP_ROW
- MUTATE_PROBABILITY_SHIFT_LEFT
- MUTATE_PROBABILITY_SHIFT_RIGHT
- UPPER_PRINT_PROBABILITY
- PSEUDO_COUNT
- PLACEMENT_OPTIONS
- SCAN_REVERSE_COMPLEMENT
