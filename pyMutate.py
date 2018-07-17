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

# Default data
resLibFile = os.environ['pyMUTATEDIR'] + '/dat/all_amino03.in'
mutMapFile = os.environ['pyMUTATEDIR'] + '/dat/pyMutateData.json'

class pyMutate():
    def __init__(self, input_pdb_path, output_pdb_path, args):
        self.input_pdb_path = input_pdb_path
        self.output_pdb_path = output_pdb_path
        self.mutationList = args.mutationList
        self.resLibFile = args.resLibFile
        self.mutMapFile = args.mutMapFile
        self.useModels = args.useModels
        self.debug = args.debug

    def launch(self):
#load data
        self.resLib = pyMutateLib.ResidueLib(self.resLibFile)
        self.mutMap = pyMutateLib.MutationMap(self.mutMapFile)

#load structure
        pdbIo = pyMutateLib.PDBManager(self.useModels)

        self.st = pdbIo.loadStructure(self.input_pdb_path)
        self.format = pdbIo.format

    #Checking models
        if len(self.st) > 1:
            if self.useModels == 'no':
                print ("#WARNING: Input Structure contains models, but using only first one due to useModels settings")
                self.useModels = False
            elif self.useModels == 'auto':
                if pdbIo.models == 'alt':
                    print ("#WARNING: Input Structure contains models, but they look as NMR models, using the first one (override with force)")
                    self.useModels = False
                else:
                    self.useModels = True
            elif self.useModels == 'force':
                if pdbIo.models == 'alt':
                    print ('#WARNING: Models found look like NMR models, but using all due to useModels = force')
                self.useModels = True
            else:
                print ("#ERROR: Unknown useModels option", file=sys.stderr)
                sys.exit(1)

            if not self.useModels:
                print ("Removing models")
                ids = []
                for md in self.st.get_models():
                    ids.append(md.id)
                for i in range(1, len(ids)):
                    self.st.detach_child(ids[i])
                self.useModels = False
        else:
            self.useModels = False

#=============================================================================
# Do Mutations
        self.muts = pyMutateLib.mutationManager(self.mutationList, self.debug)

        self.muts.checkMutations(self.st, self.debug)

        for mut in self.muts.mutList:
            mut.apply(self.st, self.mutMap, self.resLib, self.debug)
#=============================================================================
        pyMutateLib.PDBManager.saveStructure(self.st, self.output_pdb_path)
        print ("Done")

def main():
    if os.environ['pyMUTATEDIR'] == '':
        print ("WARNING: $pyMUTATEDIR not set")

    cmdline = pyMutateLib.cmdLine(defaults={'resLibFile':resLibFile, 'mutMapFile':mutMapFile})
    args = cmdline.parse_args()

    print ('==============================================')
    print ('     Simple side-chain mutation utility')
    print ('             J.L Gelpi 2018')
    print ('==============================================')

    pyMutateLib.cmdLine.printArgs(args)

    pyMutate(args.input_pdb_path, args.output_pdb_path, args).launch()

if __name__ == "__main__":
    main()
