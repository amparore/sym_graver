# Symbolic Computation of the Graver basis of a matrix

The tool `sym_hilbert` is a tool that computes the Graver basis, the Hilbert basis or the set of extremal rays of a system of equations, using a symbolic data structure (Decision Diagrams) to perform the computation.
The tool implements several variations of the algorithm described by Loïc Pottier for the computation of the Graver basis of a matrix.

The tool is similar in scope to [4ti2](https://4ti2.github.io/), but it differs from 4ti2 in several details about the algorithm, the data structures and on how the computation is carried out.

---
### Compile the tool (Unix/macOS)

To compile, the tool, you first need the [MEDDLY library](https://github.com/asminer/meddly), which implements the decision diagram data structure used throughout the tool. To do so, type these commands in a terminal:

```
git clone  https://github.com/asminer/meddly.git meddly
cd meddly
./autogen.sh
./configure --prefix=/usr/local
make
sudo make install
```

Please, refer to the Meddly site for further information on how to download, compile and install it.

To compile and install `sym_hilbert`, type
```
git clone XXXXXXX sym_hilbert
cd sym_hilbert/build
cmake .. && make
```

You should now have the tool `sym_hilbert` in the `build/` folder.

---
### How to use the `sym_hilbert` tool

The tool expects the systems of equations to be represented as a `<input>.mat` file, in the following format
```
<number_of_rows> <number_of_columns> 
<coefficients> ...
```
As an example, the following is a valid input matrix:
```
3 5
1  3  1 -2  0
0 -1  0  0  2
0  0 -2  1 -1
```
Given the coefficient's matrix, the tool can be used to compute one of the following targets:
- the Hilbert basis, saved as `<input>.hilb` (*default*);
- the Graver basis, saved as `<input>.grav` (use option **-g**);
- the set of extremal rays of the cone, saved as `<input>.xray` (use option **-x**).

The tool is called from the command line as follows:
```
./build/sym_hilbert  <options>  <input>
```

---
The tool performs several stages of the computation, and each one can be modified by options.

#### Input loading

Loads the input matrix `<input>.mat` from the disk as the input matrix **A**.

- Use **-t** to transpose the input matrix.

#### Compute the $\mathbb{Z}$-basis of the kernel of **A**

Obtain the $\mathbb{Z}$-basis from the Hermite normal form of **A**

- Use **-k** to say that the input file is already the $\mathbb{Z}$-basis.

#### Generate the Decision Diagram variable order and the pivot order

This is an essential step of every system that manipulates Decision Diagrams, as the variables (stored as levels) requires to be ordered. Such order has a huge impact in DD performances, and heuristics are tipically sued to generate *reasonable* variable orders.

- Use **-no** to disable variable ordering heuristic and use the sequential order;
- **-sloan** use Sloan heuristics (*default*) to reorder variables;
- **-sloan2** use Boost::C++ Sloan heuristics to reorder variables (requires boost::C++ when building the tool);
- **-fast** use fast reordering (usually less effective than Sloan);
- **-np** disable pivoting heuristic and use sequential order for pivoting;
- **-po** use pivoting heuristic order for both pivot and DD variable order.

#### Set target

As already said before, the default target is the Hilbert basis of **A** in the first orthant, but the target can be changed to the Graver basis (use **-g**) or the set of extremal rays of the cone of **A** in the first orthant (use **-x**).

#### Pottier algorithm

This is the main computational step of the tool, which performs the Pottier algorithm. Several options are available, divided in these groups>

- Compute by generators:
  - Use **-ne** to start from all generators at once (*default*);
  - Use **-ye** to add one generator at a time, incrementally.

- Compute by variables:
  - Use **-nl** to perform completions/reductions to all variables at once;
  - Use **-yl** to perform completions/reductions one variable at a time, following the pivot order  (*default*).

- Pottier completion/reduction algorithm:
  - Normalize all levels at once (*default*) or just the projected levels (use **-z**);
  - Perform S-Vectors operations in graded order (use **-d**) using the set of candidate vectors **C**;
  - Perform S-Vectors operations dynamically generating the candidate vectors in graded order (use **-s**).

#### Output

Once the target has been computed, it can be saved to the disk.

- Use **-o** to save the computed matrix esplicitly
- Use **-dd** to save also a representation of the Decision Diagram, as a PDF (requires the GraphViz tool);
- Use **-od** to save the computed matrix Decision Diagram in MEDDLY format;
- Use **-p** to show some performance statistics;
- Use **-c** to check the results against a saved file (for regression tests).

#### Other options

- Use **-v** and **-vv** to have verbose outputs;
- **-cs** `n` to set the Meddly cache to `n` entries;