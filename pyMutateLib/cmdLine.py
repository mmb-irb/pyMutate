import argparse

class cmdLine():
    def __init__(self, defaults):
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
            dest='input_pdb_path',
            help='PDB/mmCIF File | pdb:pdbId',
            required=True
        )

        self.argparser.add_argument(
            '-m',
            action='store',
            dest = 'mutationList',
            help='List of mutations ([chain:]OldIdNumNewId as in A:Arg232Gln, no chain or * for all chains ) | file:file_path ',
            required=True
        )

        self.argparser.add_argument(
            '--map',
            action='store',
            dest='mutMapFile',
            help='Mutation rules',
            default=defaults['mutMapFile']
        )

        self.argparser.add_argument(
            '--rlib',
            action='store',
            dest='resLibFile',
            help='Residue Lib (amber prep format)',
            default=defaults['resLibFile']
        )

        self.argparser.add_argument(
            '-o',
            action='store',
            dest='output_pdb_path',
            help='Output PDB File',
            required=True
        )
        
        self.argparser.add_argument(
            '--removeH',
            action='store',
            dest='removeH',
            help='Remove H atoms is any before make changes (no, mut, all). Recommended',
            default='mut')

    def parse_args(self):
        args = self.argparser.parse_args()
        return args

    @classmethod
    def printArgs(cls, args):

        print (' input_pdb_path: ', args.input_pdb_path)
        print (' mutations:      ', args.mutationList)
        print (' output_pdb_path:', args.output_pdb_path)
        print (' Use PDB Models: ', args.useModels)
        print (' Remove H:       ', args.removeH)
        print (' Mutation Rules: ', args.mutMapFile)
        print (' Residue Library:', args.resLibFile)
        if args.debug:
            print (' DEBUG mode on')
        print()
