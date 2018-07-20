#! /usr/bin/python3

"""
 pyCheckStructure: utility script to check structure before setup
    follows MDWeb. Standalone Version
"""

__author__ = "gelpi"
__date__ = "$13-jul-2018 15:52:55$"

import sys
import os

from pyMutateLib.cmd_line import pyCheckingCmdLine as CmdLine
from pyMutateLib.structure_manager import StructureManager
import pyMutateLib.util as util

SS_DIST = 2.5
CA_DIST = 3.8
CA_DIST_ERR = 1.
STERIC_CLASH_DIST =1.
MAX_DIST = 4.5;
STERIC_APOLAR_CLASH_DIST = 2.9
STERIC_POLAR_CLASH_DIST = 3.1
IONIC_CLASH_DIST = 3.5
apolar_elements = ["C","S"]
polar_elements = ["O", "N"]
ionic_ats = ["OD1", "OD2", "OE1", "OE2", "NZ", "NE", "NH1", "NH2"]

class pyCheckStructure():
    def __init__(self, input_pdb_path, output_pdb_path, args):
        self.input_pdb_path = input_pdb_path
        self.output_pdb_path = output_pdb_path
        self.use_models = args.use_models
        self.remove_H = args.remove_H
        self.debug = args.debug

    def launch(self):
# load structure ==============================================================
        print ("Loading structure from " + self.input_pdb_path)
        pdbdata = StructureManager()
        pdbdata.loadStructure(self.input_pdb_path, self.use_models, self.remove_H, self.debug)

# Distance related Checks ================================================================
        all_CA_distances = pdbdata.get_all_CA_distances()
        num_CA_clash=0
        num_wrong_CA_dist=0
        for pair in all_CA_distances:
            [at1,at2,dist] = pair
            if not util.same_chain(at1.get_parent(), at2.get_parent()):
                continue
            if dist < STERIC_CLASH_DIST:
                print ("CA Clash: ", util.atomid(at1), util.atomid(at2), dist)
                num_CA_clash += 1
                continue
            if util.seq_consecutive(at1.get_parent(), at2.get_parent()) and \
                (abs(dist - CA_DIST) > CA_DIST_ERR):
                print ("Incorrect CA-CA distance: ", 
                    residueid(at1.get_parent()),
                    residueid(at2.get_parent()), dist)
                num_wrong_CA_dist += 1
        if num_CA_clash == 0:
            print ("No CA Clashes detected")
        if num_wrong_CA_dist == 0:
            print ("Consecutive CA Distances ok")
        all_distances = pdbdata.get_all_distances()
#SS Bonds       
        SS_bonds=[]
        steric_clashes=[]
        steric_polar_clashes=[]
        steric_apolar_clashes=[]
        steric_ionic_clashes=[]
        for pair in all_distances:
            [at1,at2,dist] = pair
            if not util.same_model(at1.get_parent(), at2.get_parent()):
                continue
            if at1.id == "SG" and at2.id == "SG":
                if dist < SS_DIST:
                    SS_bonds.append(pair)
            elif not util.same_residue(at1,at2) and not util.seq_consecutive(at1.get_parent(),at2.get_parent()):
                if dist < STERIC_CLASH_DIST:
                    steric_clashes.append(pair)
                if dist < STERIC_APOLAR_CLASH_DIST and at1.element in apolar_elements and at2.element in apolar_elements:
                    steric_apolar_clashes.append(pair)
                if dist < STERIC_POLAR_CLASH_DIST and at1.element in polar_elements and at2.element in polar_elements:
                    steric_polar_clashes.append(pair)
                if dist < IONIC_CLASH_DIST and at1.id in ionic_ats and at2.id in ionic_ats:
                    steric_ionic_clashes.append(pair)
        if len(SS_bonds):            
            print ("Possible SS Bonds")
            for pair in SS_bonds:
                print ("    ", util.atomid(pair[0]), util.atomid(pair[1]), pair[2])
        else:
            print ("No SS bonds detected")
        if len(steric_clashes):
            print ("Severe Steric clashes")
            for pair in steric_clashes:
                print ("   ",util.atomid(pair[0]), util.atomid(pair[1]), pair[2])
        else:
            print ("No severe steric clashes detected")
        if len(steric_apolar_clashes):
            print ("Steric apolar clashes")
            for pair in steric_apolar_clashes:
                print ("   ",util.atomid(pair[0]), util.atomid(pair[1]), pair[2])
        else:
            print ("No steric apolar clashes detected")
        if len(steric_polar_clashes):
            print ("Steric polar clashes or HBonds")
            for pair in steric_polar_clashes:
                print ("   ", util.atomid(pair[0]), util.atomid(pair[1]), pair[2])
        else:
            print ("No steric polar clashes detected")
        if len(steric_ionic_clashes):
            print ("Steric ionic clashes or salt bridges")
            for pair in steric_ionic_clashes:
                print ("   ", util.atomid(pair[0]), util.atomid(pair[1]), pair[2])
        else:
            print ("No steric ionic clashes detected")
            
        
        

# Save output =================================================================
        print ("Saving final structure to " + self.output_pdb_path)
        pdbdata.saveStructure(self.output_pdb_path)
        print ("Done")

def main():

    if os.getenv('pyMUTATEDIR') == None:
        print ("WARNING: $pyMUTATEDIR not set")
        defaults={'resLib_path' : '', 'mutMap_path' : ''}

    else:
# Default data
        defaults={
            'resLib_path' : os.getenv('pyMUTATEDIR') + '/dat/all_amino03.in',
            'mutMap_path' : os.getenv('pyMUTATEDIR') + '/dat/pyMutateData.json'
        }

    cmd_line = CmdLine(defaults)
    args = cmd_line.parse_args()

    print ('==============================================')
    print ('     PDB checking utility for MD Setup')
    print ('           J.L Gelpi 2018')
    print ('==============================================')

    CmdLine.printArgs(args)

    pyCheckStructure(args.input_pdb_path, args.output_pdb_path, args).launch()

if __name__ == "__main__":
    main()
