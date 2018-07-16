##pyMutate

Simple replacement for scwrl to mutate side chains.
-New CB atoms are built from backbone
-New atoms are built using AMBER03 Prep Lib, except for PRO

usage: pyMutate [-h] [--debug] [--useModels USEMODELS] [-i PDB_PATH]
                [-m MUTATIONLIST] [--map MUTATIONMAP] [--rlib RESIDUELIB]
                [-o OUTPUT_PDB_PATH]

arguments:
  -h, --help            show this help message and exit
  --debug, -d           Produce DEBUG output
  --useModels USEMODELS
                        Use all Models [no, auto, force]
  -i PDB_PATH           PDB File | pdb:pdbId
  -m MUTATIONLIST       List of mutations ([chain:]OldIdNumNewId as in
                        A:Arg232Gln, no chain or * for all chains ) |
                        file:file_path
  --map MUTATIONMAP     Mutation rules (default: dat/pyMutateData.json)
  --rlib RESIDUELIB     Residue Lib (amber prep format) (default dat/allamino03.in)
  -o OUTPUT_PDB_PATH    Output PDB File



