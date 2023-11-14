import numpy as np

class tree_node():
    def __init__(self, N):
        self.N = N

    def __repr__(self):
        return str(self.N)

class tree_struct():
    def __init__(self):
        self.origin = tree_node(0)        

    def __repr__(self):
        s="TreeStructure\n"
        s+= repr(self.origin)
        return s

def medial_axis_transform(domain):
    #kill homology
    print("*")

Tree = tree_struct()
print(Tree)