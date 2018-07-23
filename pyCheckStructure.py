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
CA_CA_THRESHOLD = 15
apolar_elements = ["C","S"]
polar_acceptor = ["O", "S"]
polar_donor = ["N"]
pos_ats = ["NZ", "NE", "NH1", "NH2"]
neg_ats = ["OD1", "OD2", "OE1", "OE2"]

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
        #all_distances = pdbdata.get_all_distances()
#SS Bonds       
        SS_bonds=[]
        steric_clashes=[]
        steric_acceptor_clashes=[]
        steric_donor_clashes=[]
        steric_apolar_clashes=[]
        steric_positive_clashes=[]
        steric_negative_clashes=[]
        for ca_pair in all_CA_distances:
            [ca1,ca2,ca_dist] = ca_pair
            if not util.same_model(ca1.get_parent(), ca2.get_parent()):
                continue
            if ca_dist > CA_CA_THRESHOLD:
                continue
#            if 'CA' in at1.get_parent() and 'CA' in at2.get_parent():
#                print (util.atomid(at1), util.atomid(at2),dist,at1.get_parent()['CA']-at2.get_parent()['CA'])
            for pair_ats in util.get_all_rr_distances(ca1.get_parent(), ca2.get_parent()):
                [at1,at2,dist] = pair_ats
                if at1.id == "SG" and at2.id == "SG":
                    if dist < SS_DIST:
                        SS_bonds.append(pair_ats)
                elif not util.same_residue(at1,at2) and not util.seq_consecutive(at1.get_parent(),at2.get_parent()):
                    if dist < STERIC_CLASH_DIST:
                        steric_clashes.append(pair_ats)
                    if dist < STERIC_APOLAR_CLASH_DIST and (at1.element in apolar_elements or at2.element in apolar_elements):
                        steric_apolar_clashes.append(pair_ats)
                    if dist < STERIC_POLAR_CLASH_DIST and at1.element in polar_donor and at2.element in polar_donor:
                        steric_donor_clashes.append(pair_ats)
                    if dist < STERIC_POLAR_CLASH_DIST and at1.element in polar_acceptor and at2.element in polar_acceptor:
                        steric_acceptor_clashes.append(pair_ats)
                    if dist < IONIC_CLASH_DIST and at1.id in pos_ats and at2.id in pos_ats:
                        steric_positive_clashes.append(pair_ats)
                    if dist < IONIC_CLASH_DIST and at1.id in neg_ats and at2.id in neg_ats:
                        steric_negative_clashes.append(pair_ats)
            
            if len(SS_bonds):            
                print (len(SS_bonds),"Possible SS Bonds")
            for pair in SS_bonds:
                print ("    ", util.atomid(pair[0]), util.atomid(pair[1]), pair[2])
        else:
            print ("No SS bonds detected")
        if len(steric_clashes):
            print (len(steric_clashes),"Severe Steric clashes")
            for pair in steric_clashes:
                print ("   ",util.atomid(pair[0]), util.atomid(pair[1]), pair[2])
        else:
            print ("No severe steric clashes detected")
        if len(steric_apolar_clashes):
            print (len(steric_apolar_clashes) , "Steric apolar clashes")
            for pair in steric_apolar_clashes:
                print ("   ",util.atomid(pair[0]), util.atomid(pair[1]), pair[2])
        else:
            print ("No steric apolar clashes detected")
        if len(steric_donor_clashes):
            print (len(steric_donor_clashes),"steric polar donor clashes")
            for pair in steric_donor_clashes:
                print ("   ", util.atomid(pair[0]), util.atomid(pair[1]), pair[2])
        else:
            print ("No steric polar donor detected")
        if len(steric_acceptor_clashes):
            print (len(steric_acceptor_clashes),"steric polar acceptor clashes")
            for pair in steric_acceptor_clashes:
                print ("   ", util.atomid(pair[0]), util.atomid(pair[1]), pair[2])
        else:
            print ("No steric polar acceptor detected")
        if len(steric_positive_clashes):
            print (len(steric_positive_clashes) , "Steric ionic positive clashes")
            for pair in steric_positive_clashes:
                print ("   ", util.atomid(pair[0]), util.atomid(pair[1]), pair[2])
        else:
            print ("No steric ionic positive clashes detected")
        if len(steric_negative_clashes):
            print (len(steric_negative_clashes) , "Steric ionic negative clashes")
            for pair in steric_negative_clashes:
                print ("   ", util.atomid(pair[0]), util.atomid(pair[1]), pair[2])
        else:
            print ("No steric ionic negative clashes detected")
            
        
        

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
