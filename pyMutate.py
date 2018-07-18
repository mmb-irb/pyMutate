#! /usr/bin/python3

"""
 pyMutate: utility script to mutate one or more side chains
           compatible with biobb. Standalone Version
"""

__author__ = "gelpi"
__date__ = "$13-jul-2018 15:52:55$"

import sys
import os

import pyMutateLib

class pyMutate():
    def __init__(self, input_pdb_path, output_pdb_path, args):
        self.input_pdb_path = input_pdb_path
        self.output_pdb_path = output_pdb_path
        self.mutation_list = args.mutation_list
        self.resLib_path = args.resLib_path
        self.mutMap_path = args.mutMap_path
        self.use_models = args.use_models
        self.debug = args.debug
        self.remove_H = args.remove_H
        if self.remove_H == 'no':
            print ("WARNING: removeH = no is not implemented (yet), using default (mut)")
            self.remove_H = 'mut'

#load data
        self.residue_lib = pyMutateLib.ResidueLib(self.resLib_path)
        self.mutation_map = pyMutateLib.MutationMap(self.mutMap_path)

    def launch(self):
# load structure ==============================================================
        print ("Loading structure from " + self.input_pdb_path)
        pdbdata = pyMutateLib.StructureManager()
        pdbdata.loadStructure(self.input_pdb_path, self.use_models, self.remove_H, self.debug)
# Do Mutations ================================================================
        self.muts = pyMutateLib.MutationManager(self.mutation_list, self.debug)

        self.muts.checkMutations(pdbdata.st, self.debug)

        for mut in self.muts.mutList:
            mut.apply(pdbdata.st, self.mutation_map, self.residue_lib, self.remove_H, self.debug)
# Save output =================================================================
        print ("Saving final structure to " + self.output_pdb_path)
        pdbdata.saveStructure(self.output_pdb_path)
        print ("Done")

def main():
    
    if os.getenv('pyMUTATEDIR') == None:
        print ("WARNING: $pyMUTATEDIR not set")
        defaults={'resLib_path' : '', 'mutMap_path' : ''}
    
    else:
# Default data
        defaults={
            'resLib_path' : os.getenv('pyMUTATEDIR') + '/dat/all_amino03.in',
            'mutMap_path' : os.getenv('pyMUTATEDIR') + '/dat/pyMutateData.json'
        }

    cmd_line = pyMutateLib.CmdLine(defaults)
    args = cmd_line.parse_args()

    print ('==============================================')
    print ('     Simple side-chain mutation utility')
    print ('             J.L Gelpi 2018')
    print ('==============================================')

    pyMutateLib.CmdLine.printArgs(args)

    pyMutate(args.input_pdb_path, args.output_pdb_path, args).launch()

if __name__ == "__main__":
    main()
