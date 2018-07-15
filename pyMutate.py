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
    
    if not hasattr(args,'id') and hasattr(pdbIo,'id'):
        args.id = pdbIo.id
    args.format = pdbIo.format
    
    print ('')
    print ('     Simple side-chain mutation utility')
    print ('             J.L Gelpi 2018')
    print ('')
    cmdline.printArgs(args)
    
    #Chekcing models
    if len(st) > 1:
        if args.useModels == 'no':
            print ("#WARNING: Input Structure contains models, but using only first one due to useModels settings")
            args.useModels = False
        elif args.useModels == 'auto':
            if pdbIo.models =='alt':
                print ("#WARNING: Input Structure contains models, but they look as NMR models, using the first one (override with force)")
                args.useModels = False
            else:
                args.useModels = True
        elif args.useModels=='force':             
            if pdbIo.models=='alt':
                print ('#WARNING: Models found look like NMR models, but using all due to useModels = force')
                args.useModels=True
        else:
            print ("#ERROR: unknown useModels option")
        if not args.useModels:
            print ("Removing models")
            ids=[]
            for md in st.get_models():
                ids.append(md.id)
            for i in range(1,len(ids)):
                st.detach_child(ids[i])
            args.useModels=False
    else:
        args.useModels=False
        
#=============================================================================
# Do Something
    muts = pyMutateLib.mutationManager()
    muts.loadMutationList(args.mutationList, args.debug)
    
    muts.checkMutations(st, args.debug)
    mutMap = pyMutateLib.MutationMap(args.mutationMap)
    
    for mut in muts.mutList:
        mut.apply(st,mutMap,args.debug)
#=============================================================================
    pdbIo.saveStructure(st,args.output_pdb_path)
    print ("Done")

