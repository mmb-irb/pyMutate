# pyMutate
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/035eb023fd9b4484977ace3dbc73acd1)](https://www.codacy.com/app/jlgelpi/pyMutate?utm_source=mmb.irbbarcelona.org&amp;utm_medium=referral&amp;utm_content=gitlab/BioExcel/pyMutate&amp;utm_campaign=Badge_Grade)
### Simple utility to mutate side chains.

## DEPRECATED. Functionality available on Biobb_structure_checking at [https://github.com/bioexcel/biobb_structure_checking]

Uses a set of rules (/dat/pyMutateData.json) to transform side chains  
New CB atoms are built from backbone  
New atoms are built using AMBER03 Prep Lib (/dat/all_amino03.in), with a modified PRO residue


`usage: pyMutate [-h] [--debug] [--useModels USEMODELS] [-i PDB_PATH]
                [-m MUTATIONLIST] [--map MUTATIONMAP] [--rlib RESIDUELIB]
                [-o OUTPUT_PDB_PATH]`

### arguments:  
 `-h, --help           show this help message and exit  
 --debug, -d           Produce DEBUG output  
 --useModels USEMODELS Use all Models [no, auto, force]  
 -i PDB_PATH           PDB File | pdb:pdbId  
 -m MUTATIONLIST       List of mutations ([chain:]OldIdNumNewId as in
                        A:Arg232Gln, no chain or * for all chains ) |
                        file:file_path  
 -o OUTPUT_PDB_PATH   Output PDB File  
--map MUTATIONMAP    Mutation rules (default: dat/pyMutateData.json)  
 --rlib RESIDUELIB    Residue Lib (amber prep format) (default dat/allamino03.in)  `
 
### Requirements
 python 3.x  
 numpy    
 BioPython (Bio.PDB)  

