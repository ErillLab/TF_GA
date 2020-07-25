# TF_GA: Bacterial Promoter Modeling with Genetic Programming

This project aims at using Genetic Programming to evolve models of regulated bacterial promoters, using datasets of upstream sequences of genes known or suspected to be regulated by a common set of transcription factors, and a control dataset of non-regulated sequences.

A tree will be the data struture to save the models (organism). In this data structure we find two main types of nodes: connectors and PSSM-type recognizers.

- Connectors: node objects that connect up to parent node and down to two nodes. A connector has a specified ideal distance between its downstream connecting elements.
- PSSM-type recognizers: node objects that specify a short (e.g. 4-5 bp) DNA recognition pattern. Based on the BioPython [motifs](https://biopython-cn.readthedocs.io/zh_CN/latest/en/chr14.html) class.

## Install with virtualenv

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


## Install with conda

If you are more confortable using conda you can import the env using the following YAML:
```yaml
name: TF_GA_env_conda
channels:
  - defaults
dependencies:
  - _libgcc_mutex=0.1=main
  - biopython=1.76=py37h7b6447c_0
  - blas=1.0=mkl
  - ca-certificates=2020.1.1=0
  - certifi=2019.11.28=py37_0
  - intel-openmp=2020.0=166
  - ld_impl_linux-64=2.33.1=h53a641e_7
  - libedit=3.1.20181209=hc058e9b_0
  - libffi=3.2.1=hd88cf55_4
  - libgcc-ng=9.1.0=hdf63c60_0
  - libgfortran-ng=7.3.0=hdf63c60_0
  - libstdcxx-ng=9.1.0=hdf63c60_0
  - mkl=2020.0=166
  - ncurses=6.1=he6710b0_1
  - numpy=1.11.3=py37h7e9f1db_12
  - numpy-base=1.11.3=py37hde5b4d6_12
  - openssl=1.1.1d=h7b6447c_4
  - pip=20.0.2=py37_1
  - python=3.7.3=h0371630_0
  - readline=7.0=h7b6447c_5
  - setuptools=45.2.0=py37_0
  - sqlite=3.31.1=h7b6447c_0
  - tk=8.6.8=hbc83047_0
  - wheel=0.34.2=py37_0
  - xz=5.2.4=h14c3975_4
  - zlib=1.2.11=h7b6447c_3

```

You can create it using the following command:
```bash
conda env create -f environment.yml
```

## Generate documentation
To generate documentation we use the `pdoc` command (not installed by default).

```bash
# Move to the source dir
cd /path/to/project/src

# Generate the new docmentation
pdoc --html .

# Remove old documentation
rm -rf ../docs/html

# Move the new documentation to the docs dir
mv html/ ../docs/

```


## Python related

Developing in 
**Python 3.7.3** 

Libraries:
- BioPython 1.76
- numpy 1.11.3 

