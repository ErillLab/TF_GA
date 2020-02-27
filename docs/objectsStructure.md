
## Main Objects

### Organism

Organisms contain basic organismal parameters (fitness, etc.) and root node. The complete tree structure of the organism is specified by the nodes instantiated in the node downstream of the root node.

Variables:
- ID: each organism will have unique identifier
- Root node: this node will support the whole tree structure

- Best fitness score (so we dont have to recompute the fitness if we already computed it)
- ~~Best ubication (Best ubication should be for a specific DNA, so this should be an array, and probably it's not optimal to save an array for every organism)~~

Methods:

- cross over and offspring given a certain organism
- mutation of an organism
- compute fitness of an organism given a set of positive and negative data

### OrganismFactory

The organism factory can create all necessary elements to make an organism(PSSMObjects and ConnectorObjects)

Variables:
- ID (we can identify each organism by ID)
 

Methods:
- get a new Organism (builds the full tree structure)
- get a new PSSM Object (creates a PSSM object)
- get a new Connector objects (creates a connector and the corresponding PSSMObjects)



### Node

Node is an abstract class designed to hold the two main types of elements that will populate the GP tree.
The Node class will contain abstract methods for: mutation, compute energy...


### P objects (PSSM-type recognizers)

Meta-parameters:

- Width: length of DNA recognition region (bp)
- ~~Depth: simulated number of sequences making up the motif (i.e. max counts)~~

Variables:
- numpy array PWM
- numpy array PSSM (not a Biopython motif, but can be reinspected)
- ~~last position (this should be an array, so better to have it in another location)~~
- ~~last score (if negative datasets change over the execution, last score will "almost" never have the same value)~~
Methods:

- Mutation sub-methods
  - Reset counts (to random values between 0 and _depth_, adding up to _depth_)
  - Randomize column counts (for one column, selected at random)
  - Flip two columns
  - Flip two rows
  - Shift columns left/right (adding new random column with uniform base frequencies)
- Operational methods
  - bio-motif methods: pwm, pssm, score...
  - compute energy (return (best) binding energy (i.e. score) and position on sequence)
  

### C objects (Connectors)

Variables:

- mu and sigma: mean "ideal" distance and "variance" (smoothness) for distance among connected element
- Two connected nodes

Methods:

- Mutation sub-methods
  - Point mutation (alter C μ or σ by rand(sign) * int(rand(τMC)))
  - **Swap connecting elements within connector**

- Operational methods
  - compute energy (return binding energy, given connected elements, and position on sequence)



## Non-object functions and main variables

Variables:
  - Unmutable positive DNA dataset
  - Unmutable negative DNA dataset
  - Collection of organisms (population)
  - a matrix to save the best position of a given PSSMObject to a DNA sequence **(To check: this could save some time)**
Methods:
  - generate a random organism (can be called multiple times to get the initial population)
  - recombine two organisms (selecting a random node from each organism and swapping pointers)
  - get the similarity between an organism and parents(return similarity to parent A and parent B)


## Basic loop
  - loop while not a good organism:
    - random sort by pairs
    - for every pair
      - Recombine organisms with a given probability pr (2 additional organisms + 2 parents)
      - Mutate children with a given probability pm
      - Compute similarity and match children with his closest parent
      - Best parent/offspring survives (based on fitness function) (2 winners)






