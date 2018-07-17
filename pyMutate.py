#! /usr/bin/python

# 
# pyMutate: simple script to mutate one or more side chains
#           compatible with biobb . Standalone Version        
# 
__author__ = "gelpi"
__date__ = "$13-jul-2018 15:52:55$"

import sys
import pyMutateLib

# Setting required for standalone use
# 
# homeDir = "PATH_TO_APPDIR"
#
homeDir = "/home/gelpi/data/DEVEL/BioExcel/pyMutate"

# Default data
resLibFile = homeDir + '/dat/all_amino03.in'
mutMapFile = homeDir + '/dat/pyMutateData.json'

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
        resLib = pyMutateLib.ResidueLib(self.resLibFile)
        mutMap = pyMutateLib.MutationMap(self.mutMapFile)

#load structure
        pdbIo = pyMutateLib.PDBManager(self.useModels)

        st = pdbIo.loadStructure(self.input_pdb_path)
        self.format = pdbIo.format

    #Checking models
        if len(st) > 1:
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
        if not self.useModels:
            print ("Removing models")
            ids = []
            for md in st.get_models():
                ids.append(md.id)
            for i in range(1, len(ids)):
                st.detach_child(ids[i])
            self.useModels = False
        else:
            self.useModels = False
        
#=============================================================================
# Mutations
        muts = pyMutateLib.mutationManager(self.mutationList, self.debug)

        muts.checkMutations(st, self.debug)
    
        for mut in muts.mutList:
            mut.apply(st, mutMap, resLib, self.debug)
#=============================================================================
        pdbIo.saveStructure(st, self.output_pdb_path)
        print ("Done")

def main():

    cmdline = pyMutateLib.cmdLine(defaults={'resLibFile':resLibFile, 'mutMapFile':mutMapFile})
    args = cmdline.parse_args()
    
    print ('==============================================')
    print ('     Simple side-chain mutation utility')
    print ('             J.L Gelpi 2018')
    print ('==============================================')
    cmdline.printArgs(args)
    
 
    pyMutate(args.input_pdb_path, args.output_pdb_path, args).launch()

if __name__ == "__main__":
    main()
