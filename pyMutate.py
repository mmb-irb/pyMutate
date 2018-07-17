#! /usr/bin/python

#
# pyMutate: simple script to mutate one or more side chains
<<<<<<< HEAD
#           compatible with biobb         
# 
__author__ = "gelpi"
__date__ = "$13-jul-2018 15:52:55$"

=======
#           compatible with biobb . Standalone Version
#
__author__ = "gelpi"
__date__ = "$13-jul-2018 15:52:55$"
>>>>>>> master

import sys
import pyMutateLib

<<<<<<< HEAD
from biobb_common.configuration import  settings
from biobb_common.tools import file_utils as fu
from biobb_common.command_wrapper import cmd_wrapper

# Setting required for standalone use
# 
# homeDir = "PATH_TO_APPDIR"
#
homeDir = "/home/gelpi/data/DEVEL/BioExcel/pyMutate"


=======
# Setting required for standalone use
#
# homeDir = "PATH_TO_APPDIR"
#
homeDir = "/home/gelpi/data/DEVEL/BioExcel/pyMutate"

>>>>>>> master
# Default data
resLibFile = homeDir + '/dat/all_amino03.in'
mutMapFile = homeDir + '/dat/pyMutateData.json'

class pyMutate():
<<<<<<< HEAD
    def __init__(self,input_pdb_path, output_pdb_path, properties):
        self.input_pdb_path = input_pdb_path
        self.output_pdb_path = output_pdb_path
        self.properties = properties
        self.global_log=properties.get('global_log', None)
        self.prefix=properties.get('prefix',None)
        self.step=properties.get('step', None)
        self.path=properties.get('path', None)
        self.mutation = properties.get('mutation', None)
        self.resLib = properties.get('resLib')
        self.mutMap = properties.get('mutMap')
        self.useModels = properties.get('useModels', None)
    
    def launch(self):
        out_log, err_log = fu.get_logs(path=self.path, prefix=self.prefix, step=self.step)

        pdbIo = pyMutateLib.PDBManager(self.useModels)
     
        st=pdbIo.loadStructure(self.input_pdb_path)
        self.format = pdbIo.format
    
           
    #Chekcing models
        if len(st) > 1:
            if self.useModels == 'no':
                out_log.info ("#WARNING: Input Structure contains models, but using only first one due to useModels settings")
                self.useModels = False
            elif self.useModels == 'auto':
                if pdbIo.models =='alt':
                    out_log.info ("#WARNING: Input Structure contains models, but they look as NMR models, using the first one (override with force)")
                    self.useModels = False
                else:
                    self.useModels = True
            elif self.useModels=='force':             
                if pdbIo.models=='alt':
                    out_log.info ('#WARNING: Models found look like NMR models, but using all due to useModels = force')
                    self.useModels=True
        else:
            err_log.error ("Unknown useModels option")
        if not self.useModels:
            out_log.info ("Removing models")
            ids=[]
            for md in st.get_models():
                ids.append(md.id)
            for i in range(1,len(ids)):
                st.detach_child(ids[i])
            self.useModels=False
        else:
            self.useModels=False
        
#=============================================================================
        muts = pyMutateLib.mutationManager()
        muts.loadMutationList(self.mutation)
    
        resLib = pyMutateLib.ResidueLib(self.resLib)
    
        muts.checkMutations(st)
        mutMap = pyMutateLib.MutationMap(self.mutMap)
    
        for mut in muts.mutList:
            mut.apply(st, mutMap, resLib)
#=============================================================================
        out_log.info ("Writing modified structure to "+ self.output_pdb_path)
        pdbIo.saveStructure(st,self.output_pdb_path)
        out_log.info ("Done")


def main():

    cmdline = pyMutateLib.cmdLine(defaults={'resLibFile':resLibFile, 'mutMapFile':mutMapFile})
    args=cmdline.parse_args()
 
    properties = {
        'mutation': args.mutationList,
        'resLib': args.resLibFile,
        'mutMap': args.mutMapFile,
        'useModels':  args.useModels
    }
    pyMutate(args.input_pdb_path, args.output_pdb_path, properties).launch()
=======
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

    cmdline = pyMutateLib.cmdLine(defaults={'resLibFile':resLibFile, 'mutMapFile':mutMapFile})
    args = cmdline.parse_args()

    print ('==============================================')
    print ('     Simple side-chain mutation utility')
    print ('             J.L Gelpi 2018')
    print ('==============================================')

    pyMutateLib.cmdLine.printArgs(args)


    pyMutate(args.input_pdb_path, args.output_pdb_path, args).launch()
>>>>>>> master

if __name__ == "__main__":
    main()
