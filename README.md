# TF_GA

The main idea of the project is to create a model of nitrogenous bases that fit in our datasets. We create new models using a genetic algorithm and then evaluate how they fit into our data.

A tree will be the data struture to save the models (organism). In this data structure we find connectors and PSSM-type recognizers:

- Connectors: Is an object that connects a node with other two nodes and an ideal distance between them.
- PSSM-type recognizers: Is an object that saves a number of bases of the DNA.
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

Saves a main root connector and a tree-like structure of the model using connector. Leaves of the tree are PSSM-type and Shape-type recognizers.

### P objects (PSSM-type recognizers)

Variables:

Methods:


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

