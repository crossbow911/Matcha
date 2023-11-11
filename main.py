import numpy as np
import matplotlib.pyplot as plt
from utilities import *
from geometries import *

if __name__=="__main__":
    # ''' Test cases '''
    # fix, ax = plt.subplots(figsize=(10,4))
    # t=np.linspace(0,1,10000)

    # arc1 = Arc(4+5j, 8+2.5j, 2.4)
    # ax.plot(arc1.at(t).real, arc1.at(t).imag)

    # bez1 = Bezier([4+2.8j, 7+2.j, 0.5+4.4j, 6+4j])
    # ax.plot(bez1.at(t).real, bez1.at(t).imag)

    # line1 = Line([4+5j, 4.8+2j])
    # ax.plot(line1.at(t).real, line1.at(t).imag)
                

    # intersec = find_intersection(arc1, bez1)
    # print("arc1, bez1: ", intersec)
    # for p in intersec:
    #     ax.plot(p.real, p.imag, "ro", markersize=4)

    # intersec = find_intersection(arc1, line1)
    # print("arc1, line1: ", intersec)
    # for p in intersec:
    #     ax.plot(p.real, p.imag, "ro", markersize=4)
        
    # intersec = find_intersection(line1, bez1)
    # print("line1, bez1: ", intersec)
    # for p in intersec:
    #     ax.plot(p.real, p.imag, "ro", markersize=4)
        
    # ax.set_aspect("equal")
    # plt.show()

    fix, ax = plt.subplots(figsize=(10,4))
    t=np.linspace(0,1,10000)

    line1 = Line([1+1j,3+2j])
    point = 2+1j
    ax.plot(line1.at(t).real, line1.at(t).imag)
    ax.plot(point.real, point.imag,"ro")
    ax.plot([point.real,line1.closest_point(point).real], [point.imag,line1.closest_point(point).imag],"r--")

    ax.set_aspect("equal")
    plt.show()
