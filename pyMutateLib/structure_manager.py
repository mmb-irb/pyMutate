"""
    StructureManager: module to handle structure data.
"""

import math
import sys

from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBList import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser

MODELS_MAXRMS = 5.0    # Threshold value to detect NMR models (angs)

class StructureManager():

    def __init__(self):
        self.model_type = ''
        self.num_ats = 0

    def loadStructure(self, input_pdb_path, use_models, remove_H, debug=False):

        if "pdb:"in input_pdb_path:
            pdbl = PDBList(pdb='tmpPDB')

            try:
                input_pdb_path = input_pdb_path[4:].upper()
                real_pdb_path = pdbl.retrieve_pdb_file(input_pdb_path)
                parser = MMCIFParser()
                input_format = 'cif'

            except IOError:
                print ("#ERROR: fetching PDB " + input_pdb_path, file=sys.stderr)
                sys.exit(2)
        else:
            real_pdb_path = input_pdb_path
            if '.pdb' in real_pdb_path:
                parser = PDBParser(PERMISSIVE=1)
                input_format = 'pdb'
            elif '.cif' in real_pdb_path:
                parser = MMCIFParser()
                input_format = 'cif'
            else:
                print ('#ERROR: unknown filetype', file=sys.stderr)
                sys.exit(2)
        try:
            self.st = parser.get_structure('st', real_pdb_path)

        except OSError:
            print ("#ERROR: parsing PDB", file=sys.stderr)
            sys.exit(2)

        #====== Internal residue renumbering =========================================
        i = 1
        for r in self.st.get_residues():
            r.index = i
            i += 1

        #Atom renumbering for mmCIF,
        if input_format == 'cif':
            i = 1
            for at in self.st.get_atoms(): # Check numbering in models
                at.serial_number = i
                if hasattr(at, 'selected_child'):
                    at.selected_child.serial_number = i
                i += 1

        for at in self.st[0].get_atoms():
            self.num_ats += 1

        #checking models type
        if len(self.st) > 1:

            if use_models == 'no':
                print ("WARNING: Input Structure contains models, but using only first one due to use_models settings")
                remove_models = True

            elif use_models == 'auto':
                if debug:
                    print ("DEBUG: Found " + str(len(self.st)) + " models")
                    print ("DEBUG: RMS " + str(StructureManager.calcRMSdAll(self.st[0], self.st[1])))

                if StructureManager.calcRMSdAll(self.st[0], self.st[1]) < MODELS_MAXRMS:
                    if debug:
                        print ("DEBUG: Models look like alternative conformations, will consider only one")
                    self.model_type = 'alt'
                    print ("WARNING: Input Structure contains models, but they look as NMR models, using the first one (override with force)")
                    remove_models = True

                else:
                    self.model_type = 'traj'
                    remove_models = False

            elif use_models == 'force':
                if StructureManager.calcRMSdAll(self.st[0], self.st[1]) < MODELS_MAXRMS:
                    print ('#WARNING: Models found look like NMR models, but using all due to use_models = force')
                remove_models = False

            else:
                print ("#ERROR: Unknown use_models option", file=sys.stderr)
                sys.exit(1)

            if remove_models:
                print ("Removing models")
                ids = []

                for md in self.st.get_models():
                    ids.append(md.id)

                for i in range(1, len(ids)):
                    self.st.detach_child(ids[i])
        # Hydrogens

        if remove_H == 'all':
            print ("Removing H atoms")

            for r in self.st.get_residues():
                StructureManager.removeHFromRes(r)

    def saveStructure(self, output_pdb_path):

        pdbio = PDBIO()
        pdbio.set_structure(self.st)

        try:
            pdbio.save(output_pdb_path)

        except OSError:
            print ("#ERROR: unable to save PDB data on " + output_path, file=sys.stderr)

    @classmethod
    def calcRMSdAll (cls, st1, st2):
        ats1 = []
        ats2 = []

        for at in st1.get_atoms():
            ats1.append(at)
        for at in st2.get_atoms():
            ats2.append(at)

        rmsd = 0

        i = 0
        while i < len(ats1)and i < len(ats2):
            d = ats1[i]-ats2[i]
            rmsd = rmsd + d * d / len(ats1)
            i = i + 1

        return (math.sqrt(rmsd))

    @classmethod
    def removeHFromRes (cls,r):
        H_list=[]
        for at in r.get_atoms():
            if at.element == 'H':
                H_list.append(at.id)
        for at_id in H_list:
            r.detach_child(at_id)
