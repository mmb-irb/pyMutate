
import argparse
import sys

homeDir="/home/gelpi/data/DEVEL/BioExcel/pyMutate"
resLibFile = homeDir + '/dat/all_amino03.in'
mutMapFile = homeDir + '/dat/pyMutateData.json'

class cmdLine():
    def __init__(self, defaults=[]):
        self.argparser = argparse.ArgumentParser(
            prog='pyMutate',
        description='Basic utility to mutate protein side chains'
        )

        self.argparser.add_argument(
            '--debug', '-d',
            action='store_true',
            dest='debug',
            help='Produce DEBUG output'
        )

        self.argparser.add_argument(
            '--useModels',
            action='store',
            dest='useModels',
            help='Use all Models [no, auto, force]',
            default='auto'
        )
    
        self.argparser.add_argument(
            '-i',
            action='store',
            dest='pdb_path',
            help='PDB File | pdb:pdbId'
        )
        
        self.argparser.add_argument(
            '-m',
            action='store',
            dest = 'mutationList',
            help='List of mutations ([chain:]OldIdNumNewId as in A:Arg232Gln, no chain or * for all chains ) | file:file_path ',
        )
        
        self.argparser.add_argument(
            '--map',
            action='store',
            dest='mutationMap',
            help='Mutation rules',
            default=mutMapFile
        )
        
        self.argparser.add_argument(
            '--rlib',
            action='store',
            dest='residueLib',
            help='Residue Lib (amber prep format)',
            default=resLibFile
        )
        
        self.argparser.add_argument(
            '-o',
            action='store',
            dest='output_pdb_path',
            help='Output PDB File'
        )

    def parse_args(self):    
        args = self.argparser.parse_args()
        print (args)
        if not args.pdb_path:
            self.argparser.print_help()
            sys.exit(2)
        return args

    def printArgs(self,args):
        print ('Arguments list')
        print ('==============')
        print (' pdb_path:       ', args.pdb_path)
        print (' mutations:      ', args.mutationList)
        print (' output_pdb_path:', args.output_pdb_path)
        print (' Use PDB Models: ', args.useModels)
        print (' Mutation Rules: ', args.mutationMap)
        print (' Residue Library:', args.residueLib)
        if args.debug:
            print (' DEBUG mode on')
        print()
