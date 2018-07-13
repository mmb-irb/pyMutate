#
# Mutation Manager for pyMutate
#

import re
import sys

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
        if 'file:' in idList:
            #Load from file
            idList = idList.replace ('file:','')
            print ('#INFO: reading mutation list from file ' + idList)
            for line in open(idList,'r'):
                self.idList.append(line)
            self.idList = list(map (lambda x: x.replace('\n','').replace('\r',''), self.idList))
        else:
            self.idList=idList.split(',')
        #convert to list of mut objects
            self.idList = list(map (lambda x: Mutation(x) , self.idList))
        
    def checkMutations (self, st, stop_on_error=True):
        self.mutList=[]
        for mut in self.idList:
            self.mutList.append(mut.check(st, stop_on_error))        
    
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
    
    def check (self, st, stop_on_error=True): # Check which mutations are possible
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
    
    def apply(self, st):
        print (self.mutList)
        
        return st
                
    def __str__(self):
        return self.id


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
