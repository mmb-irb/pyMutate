# 
# mutationMap
# 

import json

class MutationMap():
    def __init__(self, file):
        try:
            fh = open (file, "r")
            jsonmap = json.load(fh)
            self.map = jsonmap['mutationMap']
        except IOError:
            print ("#ERROR: fetching MutationMap "+ file, file=sys.stderr)
            sys.exit(2)
    
    def getRules(self,aain,aaout,type):
        return self.map[aain][aaout][type]
    


