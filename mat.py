import numpy as np
from geometries import *
import time

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

def medial_axis_transform(domain, ax):
    v_node_n = 0

    #kill homology
    while len(domain.interior)>0:
        #Tree = tree_node("V"+str(v_node_n))

        # determine highest points of each hole
        highest_point = 0-np.inf*1j
        for curve in domain.interior.pop():
            p = curve.highest_point()
            if p.imag>highest_point.imag:
                highest_point=p

        r=1
        circle = Arc(highest_point, highest_point, r, angle_to_center=np.pi/2)
        n_exp=-10
        while True:
            #print("highest_point: ", highest_point, r)
            circle = Arc(highest_point, highest_point, r, angle_to_center=np.pi/2)
            #for curve in domain.exterior:
            n_intersections = 0
            for curve in domain.exterior:
                n_intersections += len(find_intersection(curve, circle))    
            if n_intersections>0: 
                r-=2**(-n_exp)
            else: 
                r+=2**(-n_exp)            
            n_exp+=1
            if n_exp==10: break
            print("r: ", r)
                # n_found=len(find_intersection(curve, circle))
            time.sleep(0.1)

        t=np.linspace(0,1,1000)
        ax.plot(circle.at(t).real, circle.at(t).imag, "r-", linewidth=1)
        


            #find contact circle

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