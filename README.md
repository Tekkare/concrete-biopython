# concrete-biopython

**_concrete-biopython_** is a FHE library based on python [**_biopython_**](https://biopython.org/) library. It implements the same class and functions when they are compatible with FHE

## Main features

- `FheSeq` class following `Bio.Seq`
- `FheMutableSeq` class following `Bio.MutableSeq`
- `SeqInterface` class interfacing `Bio.Seq`
  with `FheSeq`, and vice versa.
- `BioCircuit` class that wraps a concrete circuit and takes sequences as input

## Installation

#### Install dependencies

Use poetry to install dependencies in a local virtual environment:

```
poetry install
```

#### Set up the kernel for the jupyter tutorial

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

Please refer to the tutorial in `tutorial/tutorial.ipynb` for usage explanations.

**Note**: `FheSeq` objects are created from integer arrays within a FHE circuit. They work on encrypted arrays in the same way `Bio.Seq` objects work on clear string sequences.  
In your app, in order to go from a `Bio.Seq` object to a `FheSeq` object and vice versa, you need to use the interfacing class `SeqInterface` as explained in the tutorial. You can also simply use a `BioCircuit` which will make the interfacing internally.

To run an app `./my_app.py`, use:

```
poetry run python tests/my_app.py
```

or

```
poetry shell
python tests/my_app.py
```

## Dev

For testing:

```
poetry run python tests/test_FheSeq.py
```

or

```
poetry shell
python tests/test_FheSeq.py
```
