"""
  Module to manage mutations
"""

import re
import sys

from Bio.PDB.Atom import Atom


import pyMutateLib.util as util

ADD = 'Add'
DEL = 'Del'
MOV = 'Mov'

class MutationManager():

    def __init__(self, id_list, debug=False):
        self.mutation_list = []

        if 'file:' in id_list:
            #Load from file
            id_list = id_list.replace ('file:', '')
            print ('Reading mutation list from file ' + id_list)
            for line in open(id_list, 'r'):
                self.id_list.append(line.replace('\n', '').replace('\r', ''))

        else:
            id_list = id_list.replace(' ', '').split(',')

        #convert to list ids to Mutation instances and make the list
        self.mutation_list = list(map (Mutation, id_list))

    # Check mutation_list and unroll chains/models
    def prepareMutations (self, st, debug=False, stop_on_error=True):
        for mut in self.mutation_list:
            mut.prepare(st, debug, stop_on_error)

    # Perform the modifications according to rules
    def applyMutations (self, st, mutation_map, residue_lib, remove_H, debug):
        for mut in self.mutation_list:
            mut.apply(st, mutation_map, residue_lib, remove_H, debug)

    def __str__(self):
        return ','.join(self.id_list)

#===============================================================================

class Mutation():

    def __init__(self, mut_id):
        if ':' not in mut_id:
            mut_id = '*:' + mut_id

        self.id = mut_id.upper()

        [self.chain, mut] = mut_id.split(':')

        mut_comps = re.match('([A-z]*)([0-9]*)([A-z]*)', mut)

        self.old_id = util.residueCheck(mut_comps.group(1))
        self.new_id = util.residueCheck(mut_comps.group(3))
        self.res_num = mut_comps.group(2)

        self.id = self.chain + ":" + self.old_id + self.res_num + self.new_id

    def prepare (self, st, debug=False, stop_on_error=True): # Check which mutations are possible
        if debug:
            print ('#DEBUG: Checking ' + self.id)

        self.mutations = []
        ok = 0
        for model in st.get_models():
            if self.chain == '*':
                chains = []
                for ch in model.get_list():
                    chains.append(ch.get_id())
            else:
                chains = [self.chain]

            for ch in chains:
                if int(self.res_num) in model[ch]:
                    r = model[ch][int(self.res_num)]
                    if r.get_resname() == self.old_id:
                        self.mutations.append(
                            {
                                'model':model.get_id(),
                                'chain':ch,
                                'residue':r.get_id(),
                                'new_id':self.new_id
                            })
                        ok += 1
                    else:
                        print ('#WARNING: Unknown residue ' + ch + ':' + self.old_id + self.res_num)
                else:
                    print ("#WARNING: Unknown residue " + ch + ':' + self.old_id + self.res_num)

        if ok == 0 and stop_on_error:
            print ("#ERROR: no mutations available for " + self.id)
            sys.exit(1)

    def apply(self, st, mut_map, res_lib, remove_H, debug=False):
        if debug:
            print (self.mutations)
            print ("Mutation Rules")
            print (mut_map.mutation_map[self.old_id][self.new_id])

        for m in self.mutations:
            r = st[m['model']][m['chain']][m['residue']]
            print ("Replacing " + util.residueid(r) + " to " + self.new_id)
# Renaming ats
            for rule in mut_map.getRules(r.get_resname(), self.new_id, MOV):
                [old_at, new_at] = rule.split("-")
                print ("  Renaming " + old_at + " to " + new_at)
                for at in r.get_atoms():
                    if at.id == old_at:
                        r.detach_child(at.id)
                        at.id = new_at
                        at.element = new_at[0:1]
                        at.fullname = ' ' + new_at
                        r.add(at)
# Deleting atoms
            for at_id in mut_map.getRules(r.get_resname(), self.new_id, DEL):
                print ("  Deleting " + at_id)
                r.detach_child(at_id)
# Deleting H
            if remove_H == 'mut':
                util.removeHFromRes(r, verbose=True)

# Adding atoms (new_id required as r.resname is still the original)
            for at_id in mut_map.getRules(r.get_resname(), self.new_id, ADD):
                print ("  Adding new atom " + at_id)
                if at_id == 'CB':
                    coords = util.buildCoordsCB(r)
                else:
                    coords = util.buildCoordsOther(r, res_lib, self.new_id, at_id)

                at = Atom(
                    at_id, 
                    coords, 
                    99.0, 
                    1.0, 
                    ' ', 
                    ' ' + at_id + ' ', 
                    0, 
                    at_id[0:1]
                    )

                r.add(at)

#Renaming residue
            r.resname = self.new_id
        print ("")

    def __str__(self):
      return self.id

