# TF_GA: Bacterial Promoter Modeling with Genetic Programming

This project aims at using Genetic Programming to evolve models of regulated bacterial promoters, using datasets of upstream sequences of genes known or suspected to be regulated by a common set of transcription factors, and a control dataset of non-regulated sequences.

A tree will be the data struture to save the models (organism). In this data structure we find two main types of nodes: connectors and PSSM-type recognizers.

- Connectors: node objects that connect up to parent node and down to two nodes. A connector has a specified ideal distance between its downstream connecting elements.
- PSSM-type recognizers: node objects that specify a short (e.g. 4-5 bp) DNA recognition pattern. Based on the BioPython [motifs](https://biopython-cn.readthedocs.io/zh_CN/latest/en/chr14.html) class.

## Install

I personally use virtualenv from python to create the virtual environment.

First of all you need `python3` and `pip3` installed. Use ther version [3.7.3](https://www.python.org/downloads/).
Check everything is right running the following commands:
```bash
python3 -V
# Python 3.7.3

pip3 list
# List of packages installed

which python3
# /usr/bin/python3

which pip3
# /usr/local/bin/pip

```

Once you have everything installed, install virtualenv using pip:
```bash
pip3 install virtualenv

# and check it with 
which virtualenv
``` 

Now we can create the `virtualenv`, it will generate a directory with the new env:
```bash
virtualenv -p python3 TF_GA_env

#Activate the new env with:
source TF_GA_env/bin/activate

#Check the new env with:

python -V
# Python 3.7.3

pip list
# List of packages installed

which python
# /path/to/env/bin/python3

which pip
# /path/to/env/bin/pip

#Deactivate the new env with:
deactivate

```
And to install the necessary packages:
```
# pip install PACKAGE
```

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

## To Do
- [x] Basic documentation
- [ ] Implement GA
- [ ] Fitness of an organism

## Python related

Developing in 
**Python 3.7.3** 

Libraries:
- 

