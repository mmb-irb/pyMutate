"""
  Module to manage mutations
"""

from Bio.PDB.Atom import Atom
import numpy as np
from numpy import cos
from numpy import pi
from numpy import sin
from numpy.linalg import norm
from pyMutateLib.structure_manager import StructureManager
import re
import sys

# TODO: replace by Bio.PDB equivalent
one_letter_residue_code = {
    'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G',
    'HIS':'H', 'HID':'H', 'HIE':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
    'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S',
    'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
}

three_letter_residue_code = {
    'A':'ALA', 'C': 'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY',
    'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN',
    'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER',
    'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'
}

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

        self.old_id = _residueCheck(mut_comps.group(1))
        self.new_id = _residueCheck(mut_comps.group(3))
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
            print ("Replacing " + _residueid(r) + " to " + self.new_id)
# Renaming ats
            for r in mut_map.getRules(r.get_resname(), self.new_id, MOV):
                [old_at, new_at] = r.split("-")
                print ("  Renaming " + old_at + " to " + new_at)
                for at in res.get_atoms():
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
                print ("  Deleting H atoms ")
                StructureManager.removeHFromRes(r)

# Adding atoms (new_id required as r.resname is still the original)
            for at_id in mut_map.getRules(r.get_resname(), self.new_id, ADD):
                print ("  Adding new atom " + at_id)
                if at_id == 'CB':
                    coords = _buildCoordsCB(r)
                else:
                    coords = _buildCoordsOther(r, res_lib, self.new_id, at_id)

                at = Atom(at_id, coords, 99.0, 1.0, ' ', ' ' + at_id + ' ', 0, at_id[0:1])

                r.add(at)

#Renaming residue
            r.resname = self.new_id
        print ("")

    def __str__(self):
      return self.id

#==============================================================================
def _residueid(r):
    return r.get_resname() + " " \
        + str(r.get_parent().id) \
        + str(r.id[1]) + "/" \
        + str(r.get_parent().get_parent().id)

def _residueCheck(r):
    r = r.upper()
    id = ''
    if r in three_letter_residue_code.keys():
        id = three_letter_residue_code[r]
    elif r in one_letter_residue_code.keys():
        id = r
    else:
        print ('#ERROR: unknown residue id ' + r)
        sys.exit(1)

    return id

def _buildCoordsOther(r, res_lib, new_res, at_id):

    resid_def = res_lib.residues[new_res]
    i = 1
    while resid_def.ats[i].id != at_id and i < len(resid_def.ats):
        i = i + 1
    if resid_def.ats[i].id == at_id:
        return _buildCoords(
            r[resid_def.ats[resid_def.ats[i].link[0]].id].get_coord(),
            r[resid_def.ats[resid_def.ats[i].link[1]].id].get_coord(),
            r[resid_def.ats[resid_def.ats[i].link[2]].id].get_coord(),
            resid_def.ats[i].geom
            )
    else:
        print ("#ERROR: Unknown target atom")
        sys.exit(1)

def _buildCoordsCB(r): # Get CB from Backbone

    return _buildCoords(
        r['CA'].get_coord(),
        r['N'].get_coord(),
        r['C'].get_coord(),
        [1.5, 115.5, -123.]
    )

def _buildCoords(avec, bvec, cvec, geom):

    dst = geom[0]
    ang = geom[1] * pi / 180.
    tor = geom[2] * pi / 180.0

    v1 = avec-bvec
    v2 = avec-cvec

    n = np.cross(v1, v2)
    nn = np.cross(v1, n)

    n /= norm(n)
    nn /= norm(nn)

    n *= -sin(tor)
    nn *= cos(tor)

    v3 = n + nn
    v3 /= norm(v3)
    v3 *= dst * sin(ang)

    v1 /= norm(v1)
    v1 *= dst * cos(ang)

    return avec + v3 - v1
