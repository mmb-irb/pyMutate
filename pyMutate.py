#! /usr/bin/python

# 
# pyMutate: simple script to mutate one or more side chains
#           compatible with biobb         
# 

import sys
import pyMutateLib

__author__ = "gelpi"
__date__ = "$13-jul-2018 15:52:55$"

if __name__ == "__main__":
    cmdline = pyMutateLib.cmdLine()
    args=cmdline.parse_args()
# ============================================================================
    pdbIo = pyMutateLib.PDBManager(args)
    
    st=pdbIo.loadStructure(args.pdb_path, args.debug)
    
    if not hasattr(args,'id') and hasattr(loader,'id'):
        args.id = pdbIo.id
    args.format = pdbIo.format
    
    print ('#HEADER')
    print ('#HEADER Simple side-chain mutation utility')
    print ('#HEADER J.L Gelpi 2018')
    print ('#HEADER')
    cmdline.printArgs(args)
    
#=============================================================================
# Do Something


#=============================================================================
    pdbIo.saveStructure(st,args.output_pdb_path)

