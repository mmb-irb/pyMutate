"""
  Module to manage command line arguments
"""

import argparse

class CmdLine():
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
            '--use_models',
            action='store',
            dest='use_models',
            help='Use structure models [no, auto, force]',
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
            dest = 'mutation_list',
            help='List of mutations ([chain:]OldIdNumNewId as in A:Arg232Gln, no chain or * for all chains ) | file:file_path ',
            required=True
        )

        self.argparser.add_argument(
            '--map',
            action='store',
            dest='mutMap_path',
            help='Mutation rules',
            default=defaults['mutMap_path']
        )

        self.argparser.add_argument(
            '--rlib',
            action='store',
            dest='resLib_path',
            help='Residue Lib (amber prep format)',
            default=defaults['resLib_path']
        )

        self.argparser.add_argument(
            '-o',
            action='store',
            dest='output_pdb_path',
            help='Output PDB File',
            required=True
        )
        
        self.argparser.add_argument(
            '--remove_H',
            action='store',
            dest='remove_H',
            help='Remove H atoms is any before make changes (no, mut, all). Recommended',
            default='mut')

    def parse_args(self):
        args = self.argparser.parse_args()
        return args

    @classmethod
    def printArgs(cls, args):

        print (' input_pdb_path: ', args.input_pdb_path)
        print (' mutations:      ', args.mutation_list)
        print (' output_pdb_path:', args.output_pdb_path)
        print (' Use PDB Models: ', args.use_models)
        print (' Remove H:       ', args.remove_H)
        print (' Mutation Rules: ', args.mutMap_path)
        print (' Residue Library:', args.resLib_path)
        if args.debug:
            print (' DEBUG mode on')
        print()
