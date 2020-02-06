# TF_GA: Bacterial Promoter Modeling with Genetic Programming

This project aims at using Genetic Programming to evolve models of regulated bacterial promoters, using datasets of upstream sequences of genes known or suspected to be regulated by a common set of transcription factors, and a control dataset of non-regulated sequences.

A tree will be the data struture to save the models (organism). In this data structure we find two main types of nodes: connectors and PSSM-type recognizers.

- Connectors: node objects that connect up to parent node and down to two nodes. A connector has a specified ideal distance between its downstream connecting elements.
- PSSM-type recognizers: node objects that specify a short (e.g. 4-5 bp) DNA recognition pattern. Based on the BioPython [motifs](https://biopython-cn.readthedocs.io/zh_CN/latest/en/chr14.html) class.

## Main Objects

### Organism

Organisms contain basic organismal parameters (fitness, etc.) and root node. The complete tree structure of the organism is specified by the nodes instantiated in the node downstream of the root node.

### P objects (PSSM-type recognizers)

Meta-parameters:

- Width: length of DNA recognition region (bp)
- Depth: simulated number of sequences making up the motif (i.e. max counts)

Variables:

- the Bio.motifs object

Methods:

- Reset counts (to random values between 0 and _depth_, adding up to _depth_)
- Randomize column counts

### C objects (Connectors)

Variables:

Methods:



## To Do
- [x] Basic documentation
- [ ] Implement GA
- [ ] Fitness of an organism

## Python related

Developing in 
**Python 3.7.3** 

Libraries:
- 

