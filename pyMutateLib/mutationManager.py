#
# Mutation Manager for pyMutate
#

import re
import sys
from Bio.PDB.Atom import Atom
import math
from numpy import sin, cos, pi
from numpy.linalg import norm
import numpy as np

oneLetterResidueCode = {
        'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G',
        'HIS':'H', 'HID':'H', 'HIE':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
        'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S',
        'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y'
    }

threeLetterResidueCode = {
        'A':'ALA', 'C': 'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY',
        'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN', 
        'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER',
        'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'
    }
    

class mutationManager():

    def __init__(self):
        self.idlist=[]
        self.mutList=[]

    def loadMutationList(self,idList,debug=False):
        self.idList=[]
        if 'file:' in idList:
            #Load from file
            idList = idList.replace ('file:','')
            print ('Reading mutation list from file ' + idList)
            for line in open(idList,'r'):
                self.idList.append(line)
            self.idList = list(map (lambda x: x.replace('\n','').replace('\r',''), self.idList))
        else:
            self.idList=idList.replace(' ','').split(',')
        #convert to list of mut objects
        self.idList = list(map (lambda x: Mutation(x) , self.idList))
        
        
    def checkMutations (self, st, debug=False, stop_on_error=True):
        self.mutList=[]
        for mut in self.idList:
            self.mutList.append(mut.check(st, debug, stop_on_error))        
    
    def __str__(self):
       return ','.join(self.list)

        
class Mutation():
    def __init__(self, mutId):
        if ':' not in mutId:
            mutId = '*:'+mutId
        self.id = mutId.upper()
        [self.chain, mut] = mutId.split(':')
        mutcomps = re.match('([A-z]*)([0-9]*)([A-z]*)',mut)
        self.oldid = _residueCheck(mutcomps.group(1))
        self.newid = _residueCheck(mutcomps.group(3))
        self.resNum = mutcomps.group(2)
        self.id = self.chain + ":" + self.oldid + self.resNum + self.newid        
    
    def check (self, st, debug=False, stop_on_error=True): # Check which mutations are possible
        if debug:
            print ('#DEBUG: Checking '+ self.id)
        self.mutList=[]
        ok=0
        for model in st.get_models():
            if self.chain=='*':
                chains=[]
                for ch in model.get_list():
                    chains.append(ch.get_id())
            else:
                chains = [self.chain]
            for ch in chains:
                if int(self.resNum) in model[ch]:
                    residue=model[ch][int(self.resNum)]
                    if residue.get_resname() == self.oldid:
                        self.mutList.append(
                        {
                            'model':model.get_id(),
                            'chain':ch,
                            'residue':residue.get_id(), 
                            'newRes':self.newid
                        })
                        ok = ok + 1
                    else:
                        print ('#WARNING: Unknown residue ' + ch + ':' +self.oldid + self.resNum)
                else:
                    print ("#WARNING: Unknown residue " +  ch + ':' +self.oldid + self.resNum)
        if ok == 0 and stop_on_error:
            print ("#ERROR: no mutations available for " + self.id)
            sys.exit(1)
        return self
    
    def apply(self, st, map, resLib, debug=False):
        if debug:
            print (self.mutList)
            print ("#DEBUG: Mutation Rules")
            print (map.map[self.oldid][self.newid])
        for m in self.mutList:
            res = st[m['model']][m['chain']][m['residue']]
            print ("Replacing " + _residueid(res) + " to " + self.newid)
# Renaming ats
            for r in map.getRules(res.get_resname(),self.newid, 'Mov'):
                [oldat,newat]=r.split("-")
                print ("  Renaming "+ oldat + " to " + newat)
                for at in res.get_atoms():
                    if at.id == oldat:
                        res.detach_child(at.id)
                        at.id = newat
                        at.element = newat[0:1]
                        at.fullname=' ' + newat
                        res.add(at)
                        
# delete atoms  
            for atid in map.getRules(res.get_resname(),self.newid,'Del'):
                print ("  Deleting "+ atid)
                res.detach_child(atid)
#Adding atoms
            for atid in map.getRules(res.get_resname(),self.newid,'Add'):
                print ("  Adding new atom "+ atid)
                if atid == 'CB':
                    at = Atom(atid, _buildCoordsCB(res), 99.0, 1.0, ' ', ' '+atid+' ', 0, atid[0:1])
                else:
                    at = Atom(atid, _buildCoordsOther(res,resLib,self.newid,atid), 99.0, 1.0, ' ', ' '+atid+' ', 0, atid[0:1])
                res.add(at)
                
#Renaming residue
            res.resname=self.newid
            print ("")
        
                
    def __str__(self):
        return self.id

def _residueid(res):
    return res.get_resname() + " " +str(res.get_parent().id) +  str(res.id[1]) +"/"+ str(res.get_parent().get_parent().id)

def _residueCheck(r):
    r=r.upper()
    resid=''
    if r in threeLetterResidueCode.keys():
        resid = threeLetterResidueCode[r]
    elif r in oneLetterResidueCode.keys():
            resid = r
    else:
       print ('#ERROR: unknown residue id ' + r)
       sys.exit(1)
    return resid

def _buildCoordsOther(res, resLib, newres,atid):
#    print (res.get_list())            
#    print (res.child_dict.keys())
    residDef = resLib.residues[newres]
    i=1
    while residDef.ats[i].id != atid and i<len(residDef.ats):
        i=i+1
    if residDef.ats[i].id == atid:
#        print (i)
#        print (vars(residDef.ats[i]))
#        print (residDef.ats[residDef.ats[i].link[0]].id)
#        print (residDef.ats[residDef.ats[i].link[1]].id)
#        print (residDef.ats[residDef.ats[i].link[2]].id)
#        print (residDef.ats[residDef.ats[i].link[0]],residDef.ats[i].geom[0],res[residDef.ats[residDef.ats[i].link[0]].id].get_coord())
#        print (residDef.ats[residDef.ats[i].link[1]],residDef.ats[i].geom[1],res[residDef.ats[residDef.ats[i].link[1]].id].get_coord())
#        print (residDef.ats[residDef.ats[i].link[2]],residDef.ats[i].geom[2],res[residDef.ats[residDef.ats[i].link[2]].id].get_coord())
        return _buildCoords(
            res[residDef.ats[residDef.ats[i].link[0]].id].get_coord(),
            res[residDef.ats[residDef.ats[i].link[1]].id].get_coord(), 
            res[residDef.ats[residDef.ats[i].link[2]].id].get_coord(),
            residDef.ats[i].geom
            )
    else:
        print ("#ERROR: Unknown target atom")
        sys.exit(1)

def _buildCoordsCB(res): # Get CB from Backbone
    return _buildCoords(
        res['CA'].get_coord(),
        res['N'].get_coord(),
        res['C'].get_coord(),
        [1.5, 115.5, -123.]
        )

            
def _buildCoords(avec,bvec,cvec,geom):
    dst = geom[0]
    ang = geom[1] * pi / 180.
    tor = geom[2] * pi / 180.0

    v1= avec-bvec
    v2= avec-cvec

    n = np.cross(v1,v2)
    nn = np.cross(v1,n)

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
    
    
