# TF_GA

The main idea of the project is to create a model of nitrogenous bases that fit in our datasets. We create new models using a genetic algorithm and then evaluate how they fit into our data.

A tree will be the data struture to save the models (organism). In this data structure we find connectors and PSSM-type recognizers:

- Connectors: Is an object that connects a node with other two nodes and an ideal distance between them.
- PSSM-type recognizers: Is an object that saves a number of bases of the DNA.
- ~~Shape-type recognizers~~: Is an object that saves a specific shape of the DNA.

## Main Objects

### Organism

Saves a main root connector and a tree-like structure of the model using connector. Leaves of the tree are PSSM-type and Shape-type recognizers.

### P objects (PSSM-type recognizers)

Variables:

Methods:


### C objects (Connectors)

Variables:

Methods:


### ~~S objects~~

## To Do
- [x] Basic documentation
- [ ] Implement GA
- [ ] Fitness of an organism

## Python related

Developing in 
**Python 3.7.3** 

Libraries:
- 

