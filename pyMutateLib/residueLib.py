#
# Class to read AMBER residue lib
#
import sys
import re

class ResidueLib():
    def __init__(self, file_path):
        try:
            fh=open(file_path,"r")
        except IOError:
            print ("#ERROR: unable to open "+ file_path )
            sys.exit(1)

        line=fh.readline()
        line=fh.readline()

        self.residues={}
        atgroup=False
        chgroup=False
        imgroup=False
        res=ResidueDef()
        for line in fh:
            line = line.replace("\n","").replace("\r","")
            if line=='':
                continue
            elif line == "DONE":
                self.residues[res.id]=res
                atgroup=False
                chgroup=False
                imgroup=False
                res=ResidueDef()
            elif line=='CHARGE':
                atgroup=False
                chgroup=True
            elif line=='IMPROPER':
                chgroup=False
                imgroup=True
            elif re.match(' (...)  INT',line):
                residstr=re.match(' (...)  INT',line)
                res.id=residstr.group(1)
                atgroup=True
            elif atgroup:
                data = line.split()
                #print (len(data))
                if len(data) <11:
                    continue
                at=AtomDef(data)
                res.ats.append(at)
            elif chgroup:
                for c in line.split():
                    res.charges.append(c)
            elif imgroup:
                res.improper.append(line)

class ResidueDef():
    def __init__(self):
        self.id=''
        self.name=''
        self.ats=['']
        self.improper=[]
        self.charges=[]

class AtomDef():
    def __init__(self,data):
        self.id=data[1]
        self.type=data[2]
        self.branch=data[3]
        self.link=[int(data[4]), int(data[5]), int(data[6])]
        self.geom=[float(data[7]), float(data[8]), float(data[9])]
        self.coord=[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
        self.ch=float(data[10])
    def __str__(self):
        return self.id
