# Sparse Boolean Matrix Algebra with SuiteSparse:GraphBLAS
In the folder `include` you can find `bm_baselineGB` class to operate with boolean matrices.
## Application: Regular Path Queries on Compressed Adjacency Matrices

The file  `include/rpq_solver` implements the transformation of Regular Path Queries to algebra operations on boolean matrices.

### Instructions
installation

```Bash
git clone https://github.com/adriangbrandon/rpq-matrix.git
cd rpq-matrix
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```
### Example

From your dataset obtain following values:
- *N* = Triples
- *P* = Predicates
- *V* = Vertices

Change 10th line in `include/rpq_solver.hpp` on your *P* value that was obtained on previos step

Generate all neceasssary files for your dataset and store them into one specific folder. In this folder you should have:
- `dataset.dat.SO` file
- `dataset.dat.P` file
- `dataset.dat.baseline-64` folder with matrices for predicates (`0001.mat` - `P.mat` files)

If your dataset named `dataset.dat` and stored in `dataset` folder, then you can run solver on queries form `queries/paths.tsv` with this command:
```Bash
./build/baseline_queryGB dataset/dataset.dat ./queries/paths.tsv P N