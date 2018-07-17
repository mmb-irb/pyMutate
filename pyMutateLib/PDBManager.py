
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBList import PDBList
from Bio.PDB.PDBIO import PDBIO


import math
import sys

MODELS_MAXRMS=5.0

class PDBManager():
    def __init__(self,useModels):
        self.models=useModels

    def loadStructure(self, pdb_path, debug=False):
        if "pdb:"in pdb_path:
            pdbl= PDBList(pdb='tmpPDB')
            try:
                pdb_path = pdb_path[4:].upper()
                self.id=pdb_path
                self.pdb_path=pdbl.retrieve_pdb_file(pdb_path)
                self.parser = MMCIFParser()
                self.useChains=True
                self.format='cif'
            except IOError:
                print ("#ERROR: fetching PDB "+ pdb_path, file=sys.stderr)
                sys.exit(2)
        else:
            self.pdb_path = pdb_path
            if '.pdb' in pdb_path:
                self.parser = PDBParser(PERMISSIVE=1)
                self.format='pdb'
            elif '.cif' in pdb_path:
                self.parser = MMCIFParser()
                self.format='cif'
            else:
                print ('#ERROR: unknown filetype', file=sys.stderr)
                sys.exit(2)

        try:
            st = self.parser.get_structure('st', self.pdb_path)
        except OSError:
            print ("#ERROR: parsing PDB", file=sys.stderr)
            sys.exit(2)
        #====== Internal residue renumbering =========================================
        i=1
        for r in st.get_residues():
            r.index = i
            i=i+1
        #Atom renumbering for mmCIF,
        if self.format == 'cif':
            i=1
            for at in st.get_atoms(): # Check numbering in models
                at.serial_number = i
                if hasattr(at,'selected_child'):
                    at.selected_child.serial_number=i
                i=i+1
        self.numAts=0

        for at in st[0].get_atoms():
            self.numAts=self.numAts+1
        #checking models
        if len(st) > 1 and self.models=='auto':
            if debug:
                print ("#DEBUG: Found "+str(len(st))+" models")
                print ("#DEBUG: RMS " + PDBManager.calcRMSdAll(st[0],st[1]))
            if PDBManager.calcRMSdAll(st[0],st[1]) < MODELS_MAXRMS:
                if debug:
                    print ("#DEBUG: Models look like alternative conformations, will consider only one")
                self.models='alt'
        return st


    @classmethod
    def saveStructure(cls, st, output_path):
        pdbio=PDBIO()
        pdbio.set_structure(st)
        try:
            pdbio.save(output_path)
        except OSError:
            print ("#ERROR: unable to save PDB data on "+ output_path, file=sys.stderr)

    @classmethod
    def calcRMSdAll (cls, st1,st2):
        ats1=[]
        ats2=[]
        for at in st1.get_atoms():
            ats1.append(at)
        for at in st2.get_atoms():
            ats2.append(at)
        rmsd=0
        i=0
        while i<len(ats1)and i <len(ats2):
            d=ats1[i]-ats2[i]
            rmsd = rmsd + d*d/len(ats1)
            i=i+1
        return (math.sqrt(rmsd))
