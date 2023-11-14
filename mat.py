import numpy as np

class tree_node():
    def __init__(self, value, level=0):
        self.value = value    
        self.children = []
        self.level=level

    def __repr__(self):
        return str(self.value)

    def traverse(self):
        # moves through each node referenced from self downwards
        nodes_to_visit = [self]
        while len(nodes_to_visit) > 0:
            current_node = nodes_to_visit.pop()
            print(current_node.level*" ", current_node.value)
            nodes_to_visit += current_node.children

    def find_value(self, value):
        # finds the node with the respective value
        nodes_to_visit = [self]
        while len(nodes_to_visit) > 0:
            current_node = nodes_to_visit.pop()
            if current_node.value==value: return current_node
            nodes_to_visit += current_node.children

    def add_child(self, value):
        self.children.append(tree_node(value, level=self.level+1))

def medial_axis_transform(domain):
    v_node_n = 0
    #kill homology
    if len(domain.interior)==0: return 0
    else:
        Tree = tree_node("V"+str(v_node_n))

        for interior in domain.interior:
            highest_y_point = -np.inf
            for curve in interior:
                h = curve.bounding_box()[1].imag
                if h>highest_y_point:
                    highest_y_point=h
            print(highest_y_point)
        # compute the highest point of each boundary;
        # sort the boundary curves according to the y-coordinate
        # of the highest point;
        # for i:=1 to the number of the inner boundary curves do
        # begin
        # find the contact circle B(pi) at qi;
        # attach B(pi) to the virtual node;
        # end;
        # end KILL_HOMOLOGY

# Tree = tree_node("V")

# Tree.add_child("A")
# Tree.add_child("B")
# Tree.add_child("C")
# Tree.children[2].add_child("D")
# Tree.children[2].add_child("E")
# Tree.traverse()
# print("found:", Tree.find_value("E"))