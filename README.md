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

## Python related

Developing in 
**Python 3.7.3** 

Libraries:
- 

