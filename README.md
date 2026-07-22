# morpho-pairwise

Fast pairwise potentials.

## Installation

The `pairwise` package can be installed with `morphopm`. Type the following into a terminal:

    morphopm install pairwise

## Use

The package can be loaded into morpho using the `import` keyword.

    import pairwise

## Manual installation 

To install manually, clone this repository onto your computer in any convenient place:

    git clone https://github.com/morpho-lang/morpho-pairwise.git

then add the location of this repository to your .morphopackages file.

    echo PACKAGEPATH >> ~/.morphopackages 
    where PACKAGEPATH is the location of the git repository.

You need to compile the extension, which you can do by cd'ing to the repository's base folder and typing

    cmake -S . -B build
    cmake --build build --config Release
    cmake --install build --config Release
