
import argparse

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
    
        self.argparser.add_argument('pdb_path',
            help='PDB File | pdb:pdbId'
        )
        
        self.argparser.add_argument(
            'mutationList',
            help='List of mutations ( [chain:]OldresIdNumNewResID as in A:Arg232Gln, no chain or * for all chains )',
        )
        
        self.argparser.add_argument('mutationMap',
            help='Mutation rules'
        )
        
        self.argparser.add_argument('residueLib',
            help='Residue Lib'
        )
        
        self.argparser.add_argument('output_pdb_path',
            help='Output PDB File'
        )
        


    def parse_args(self):    
        args = self.argparser.parse_args()
        if not args.pdb_path:
            argparser.print_help()
            sys.exit(2)
        return args

    def printArgs(self,args):
        print ('Arguments list')
        print ('==============')
        print (' pdb_path:       ', args.pdb_path)
        print (' mutations:      ', args.mutationList)
        print (' output_pdb_path:', args.output_pdb_path)
        print (' Use PDB Models: ', args.useModels)
        if args.debug:
            print (' DEBUG mode on')
        print()
