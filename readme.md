# DiffDomain

DiffDomain is a suite of numeric domains tailored for abstract interpretation of partial differential equation solvers.

It is codesigned with the PDEnclose project, but can be used externally.

## Dependencies
* cereal serialization library

## Design
DiffDomain includes three domains:
* Winterval: standard interval arithmetic
* Caffeine: affine interval arithmetic
* DualDomain: affine-interval reduced product

These domains override operations defined in the Numeric concept. 

## Using

### Building
We've tested building with CMAKE >= 3.5 and Unix makefiles. The libraries should compile with `gcc` or `clang`.

* (recommended) in the project root, make a temporary directory `build`. `cd` into `build`.
* Execute `cmake -DCMAKE_BUILD_TYPE=Release ..`
* Execute `make`

### Testing
Test binaries are generated in `project_root/out/tests`. We also provide a helper script to run all tests quickly.

### Importing
The file `release/diffdomain-release.tar.gz` includes all source files. Note that the header-only library `cereal` must be in the buildpath!

Extracting this file will provide `cmakelists` that make the Winterval, Caffeine, and DualDomain libraries avaliable to cmake projects.

