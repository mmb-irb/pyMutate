
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBList import PDBList
from Bio.PDB.PDBParser import PDBParser
import math
import sys

MODELS_MAXRMS = 5.0

class PDBManager():
    def __init__(self):
        self.models = 'no'

    def loadStructure(self, pdb_path, useModels, debug=False):
        if "pdb:"in pdb_path:
            pdbl = PDBList(pdb='tmpPDB')
            try:
                pdb_path = pdb_path[4:].upper()
                self.id = pdb_path
                self.pdb_path = pdbl.retrieve_pdb_file(pdb_path)
                self.parser = MMCIFParser()
                self.format = 'cif'
            except IOError:
                print ("#ERROR: fetching PDB " + pdb_path, file=sys.stderr)
                sys.exit(2)
        else:
            self.pdb_path = pdb_path
            if '.pdb' in pdb_path:
                self.parser = PDBParser(PERMISSIVE=1)
                self.format = 'pdb'
            elif '.cif' in pdb_path:
                self.parser = MMCIFParser()
                self.format = 'cif'
            else:
                print ('#ERROR: unknown filetype', file=sys.stderr)
                sys.exit(2)
        try:
            self.st = self.parser.get_structure('st', self.pdb_path)
        except OSError:
            print ("#ERROR: parsing PDB", file=sys.stderr)
            sys.exit(2)
        #====== Internal residue renumbering =========================================
        i = 1
        for r in self.st.get_residues():
            r.index = i
            i = i + 1
        #Atom renumbering for mmCIF,
        if self.format == 'cif':
            i = 1
            for at in self.st.get_atoms(): # Check numbering in models
                at.serial_number = i
                if hasattr(at, 'selected_child'):
                    at.selected_child.serial_number = i
                i = i + 1
        self.numAts = 0

        for at in self.st[0].get_atoms():
            self.numAts = self.numAts + 1

        #checking models type
        if len(self.st) > 1:
            if useModels == 'no':
                print ("#WARNING: Input Structure contains models, but using only first one due to useModels settings")
                self.useModels = False
            elif useModels == 'auto':
                if debug:
                    print ("#DEBUG: Found " + str(len(self.st)) + " models")
                    print ("#DEBUG: RMS " + PDBManager.calcRMSdAll(self.st[0], self.st[1]))
                if PDBManager.calcRMSdAll(self.st[0], self.st[1]) < MODELS_MAXRMS:
                    if debug:
                        print ("#DEBUG: Models look like alternative conformations, will consider only one")
                    self.models = 'alt'
                    print ("#WARNING: Input Structure contains models, but they look as NMR models, using the first one (override with force)")
                    self.useModels = False
                else:
                    self.models = 'traj'
                    self.useModels = True
            elif useModels == 'force':
                if PDBManager.calcRMSdAll(self.st[0], self.st[1]) < MODELS_MAXRMS:
                    print ('#WARNING: Models found look like NMR models, but using all due to useModels = force')
                self.useModels = True
            else:
                print ("#ERROR: Unknown useModels option", file=sys.stderr)
                sys.exit(1)
        else:
            self.useModels = False
        if not self.useModels:
            print ("Removing models")
            ids = []
            for md in self.st.get_models():
                ids.append(md.id)
            for i in range(1, len(ids)):
                self.st.detach_child(ids[i])
            self.useModels = False

    def saveStructure(self, output_path):
        pdbio = PDBIO()
        pdbio.set_structure(self.st)
        try:
            pdbio.save(output_path)
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
