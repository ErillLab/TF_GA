
## Main Objects

### Organism

Organisms contain basic organismal parameters (fitness, etc.) and root node. The complete tree structure of the organism is specified by the nodes instantiated in the node downstream of the root node.

### Node

Node is an abstract class designed to hold the two main types of elements that will populate the GP tree.
The Node class will contain abstract methods for: mutation, compute energy...

### P objects (PSSM-type recognizers)

Meta-parameters:

- Width: length of DNA recognition region (bp)
- Depth: simulated number of sequences making up the motif (i.e. max counts)

Variables:

- the Bio.motifs object

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
- connected nodes

Methods:

- Mutation sub-methods
  - Point mutation (alter C μ or σ by rand(sign) * int(rand(τMC)))
  - Swap connecting elements within connector
- Operational methods
  - compute energy (return binding energy, given connected elements, and position on sequence)


