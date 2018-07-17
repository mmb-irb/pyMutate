#! /usr/bin/python

#
# pyMutate: simple script to mutate one or more side chains
#           compatible with biobb . Standalone Version
#
__author__ = "gelpi"
__date__ = "$13-jul-2018 15:52:55$"

import sys
import os
import pyMutateLib



class pyMutate():
    def __init__(self, input_pdb_path, output_pdb_path, args):
        self.input_pdb_path = input_pdb_path
        self.output_pdb_path = output_pdb_path
        self.mutationList = args.mutationList
        self.resLibFile = args.resLibFile
        self.mutMapFile = args.mutMapFile
        self.useModels = args.useModels
        self.debug = args.debug
#load data
        self.resLib = pyMutateLib.ResidueLib(self.resLibFile)
        self.mutMap = pyMutateLib.MutationMap(self.mutMapFile)

    def launch(self):
# load structure ==============================================================
        print ("Loading structure from " + self.input_pdb_path)
        pdbdata = pyMutateLib.PDBManager()
        pdbdata.loadStructure(self.input_pdb_path, self.useModels, self.debug)
# Do Mutations ================================================================
        self.muts = pyMutateLib.mutationManager(self.mutationList, self.debug)

        self.muts.checkMutations(pdbdata.st, self.debug)

        for mut in self.muts.mutList:
            mut.apply(pdbdata.st, self.mutMap, self.resLib, self.debug)
# Save output =================================================================
        print ("Saving final structure to " + self.output_pdb_path)
        pdbdata.saveStructure(self.output_pdb_path)
        print ("Done")

def main():
    if os.getenv('pyMUTATEDIR') == None:
        print ("WARNING: $pyMUTATEDIR not set")
        defaults={'resLibFile' : '', 'mutMapFile' : ''}
    else:
# Default data
        defaults={
            'resLibFile' : os.getenv('pyMUTATEDIR') + '/dat/all_amino03.in',
            'mutMapFile' : os.getenv('pyMUTATEDIR') + '/dat/pyMutateData.json'
        }

    cmdline = pyMutateLib.cmdLine(defaults)
    args = cmdline.parse_args()

    print ('==============================================')
    print ('     Simple side-chain mutation utility')
    print ('             J.L Gelpi 2018')
    print ('==============================================')

    pyMutateLib.cmdLine.printArgs(args)

    pyMutate(args.input_pdb_path, args.output_pdb_path, args).launch()

if __name__ == "__main__":
    main()
