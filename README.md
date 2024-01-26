# morpho-pairwise

Fast pairwise potentials.

## Installation 

To install this, clone this repository onto your computer in any convenient place:

    git clone https://github.com/morpho-lang/morpho-pairwise.git

then add the location of this repository to your .morphopackages file.

    echo PACKAGEPATH >> ~/.morphopackages 
    where PACKAGEPATH is the location of the git repository.

You need to compile the extension, which you can do by cd'ing to the repository's base folder and typing

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release .. 
    make install

The package can be loaded into morpho using the `import` keyword.

    import pairwise
