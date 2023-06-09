# concrete-biopython

***concrete-biopython*** is a FHE library based on python [***biopython***](https://biopython.org/) library. It implements the same class and functions when they are compatible with FHE

## Main features
- `FheSeq` class following `Bio.Seq`
- `FheMutableSeq` class following `Bio.MutableSeq`

## Installation

#### Install dependencies
Use poetry to install dependencies in a local virtual environment:
```
cd concrete-biopython
poetry init
```

#### Sett up the kernel for the jupyter tutorial

1. Activate the poetry project's virtual environment:
```
poetry shell
```
2. Create a jupyter kernel to run the tutorial with this virtual env:
```
python -m ipykernel install --user --name=concrete-biopython
```
3. Then run `jupyter notebook`, open the `tutorial/tutorial.ipynb` and select `Kernel > Change kernel > concrete-biopython`

## Usage
Please refer to the tutorial in `tutorial/tutorial.ipynb`

## Dev

For testing:
```
poetry run python tests/test_FheSeq.py
```