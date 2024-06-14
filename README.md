# world-s-simplest-ECFP
the world's simplest ECFP. needs rdkit (for getting adjacency matrix) and numpy (for fast matrix powers)

it works because squaring an adjacency matrix or laplacian like matrix is equivalent to one morgan algorithm iteration

usage `python ecfp.py c1ccccc1` will print the on bits for benzene, change the smiles to your favorite molecule.
