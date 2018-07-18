"""
  module to load mutationMap json file
"""

import json

class MutationMap():

    def __init__(self, file_path):
        try:
            fh = open (file_path, "r")
            json_map = json.load(fh)
            self.mutation_map = json_map['mutationMap']

        except IOError:
            print ("ERROR: unable to open mutation_map "+ file_path, file=sys.stderr)
            sys.exit(2)

    def getRules(self,aa_in,aa_out,rule_group):

        return self.mutation_map[aa_in][aa_out][rule_group]



