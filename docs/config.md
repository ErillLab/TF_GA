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
  - Probability to execute the rise child mutator. It selects a random node in the organism and raises it to an upper level.
- MUTATE_PROBABILITY_NODE_MUTATION
  - Probability to  execute the node mutation mutator. It selectas a random node in the organism and applies a mutation to that specific node.
- MUTATE_PROBABILITY_SUNK_CHILD
  - Probability to execute the sunk child mutator. It selects a random node in the organism and pulls it down to an lower level.
- MIN_NODES
  - Minimum number of nodes allowed in an organism. If organism has less nodes that specified, it will apply an extra complexity penalty.
- MAX_NODES
  - Maximum number of nodes allowed in an organism. If organism has more nodes that specified, it will apply an extra complexity penalty.

## Factory

- INITIAL_CONNECTOR_PROBABILITY
  - Probability to generate a connector object on the root node of the organismi when created.
- REDUCER_PROBABILITY_FACTOR
  - Every time a we go down a level in the organism (due to connectors), the probability of generating a new connector is lowered by multiplying the actual probability by this factor.
- MIN_MU
  - Minimum starting value of Mu variable for connectors
- MAX_MU
  - Maximum starting value of Mu variable for connectors
- MIN_SIGMA
  - Minimum starting value of Sigma variable for connectors
- MAX_SIGMA
  - Maximum starting value of Sigma variable for connectors
- PWM_LENGTH
  - Number of bases the recognizer can recognize
- PWM_PROBABILITY_STEP
  - Step of the points assigned to every base. This allows to specify if final probabilities go from 0.1 to 0.1 or 0.05 to 0.05.
- PWM_PROBABILITY_BASE
  - Total points to assign to the bases in a specific mutation. When all the points are assigned points are divided by the PROBABILITY_BASE to generate a probability.
- PWM_PROBABILITY_DECIMALS
  - Number of decimals used to compute probabilities

## Connector

- MUTATE_PROBABILITY_SIGMA
  - Probability to execute the Sigma mutator. It modifies the Sigma value of the connector.
- MUTATE_PROBABILITY_MU
  - Probability to execute the Mu mutator. It modifies the Mu value of the connector.
- MUTATE_PROBABILITY_SWAP
  - Probability to execute the swap mutator. It swaps the 2 nodes.
- TAU
  - Modules the value of the connector term when computing the connector energy
- MUTATE_VARIANCE_SIGMA
  - Value to change when the Sigma mutator is applied
- MUTATE_VARIANCE_MU
  - Value to change when Mu mutator is applied
- PLACEMENT_OPTIONS
  - Number of "best" options returned by the connector

## PSSM recognizer

- MUTATE_PROBABILITY_RANDOM_COL
  - Probability to execute the random column mutator. Generates a new random column and insterts it in a random position of the PWM
- MUTATE_PROBABILITY_FLIP_COL
  - Probability to execute the flip column mutator. Selects 2 random columns in the PWM and filps it.
- MUTATE_PROBABILITY_FLIP_ROW
  - Probability to execute the flip row mutator. Selects 2 random rows in the PWM and filps it.
- MUTATE_PROBABILITY_SHIFT_LEFT
  - Probability to execute the shift left mutator. Shifts to the left the PWM with feedback.
- MUTATE_PROBABILITY_SHIFT_RIGHT
  - Probability to execute the shift right mutator. Shifts to the right the PWM with feedback.
- UPPER_PRINT_PROBABILITY
  - Probability to show upper case bases when printing a recognizer
- PSEUDO_COUNT
  - Pseudo-count added to deal with negative infinities when computing the PSSM from PWM.
- PLACEMENT_OPTIONS
  - Number of "best" options returned by the pssm recognizer
- SCAN_REVERSE_COMPLEMENT
  - True if the reverse complement of the sequence should be checked. False otherwise.
