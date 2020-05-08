# Main parameters

## Genetic algorithm parameters

### Population

- POPULATION_LENGTH:
	- The size of the genetic algorithm population. 
	- Default: `100`. 
	- Range: `>0 [int]`
- POPULATION_ORIGIN
	- determines whether initial popilation is read from file or initizalized randomly. 
	- Default: `random`.
	- Options: `random`, `file`
- POPULATION_FILL_TYPE
	- if `POPULATION_ORIGIN` is file, this determines whether population is uniformly seeded with read organims, or complemented by random organisms.
	- Default: `random`.
	- Options: `random`, `uniform`
- RECOMBINATION_PROBABILITY :0.5
	- The probability that two parents will be recombined to produce offspring

### Stopping criteria

- - END_WHILE_METHOD
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

## Directory organization

- DATASET_BASE_PATH_DIR
	- Base directory for datasets.
	- Default: `datasets/`
- RESULT_BASE_PATH_DIR
	- Base directory for results [where best organisms are reported].
	- Default: `results/`
- RESULT_TEST_BASE_PATH_DIR
	- Base directory for test results [where test organisms are assessed].
	- Default: `resultsTEST/`
- POSITIVE_FILENAME":"positive_dataset.fas",
- NEGATIVE_FILENAME":"negative_dataset.fas",
- INPUT_FILENAME":"inputOrganisms2.json",
