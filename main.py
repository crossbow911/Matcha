import numpy as np
import matplotlib.pyplot as plt
from utilities import *
from geometries import *
from mat import *
import matplotlib.cm as cm
import sys

'''
ToDo:
svgs werden nicht korrekt eingelesen. Drehrichtungen, absolute Positionen, Skalierungen
'''
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

    # fix, ax = plt.subplots(figsize=(10,4))
    # t=np.linspace(0,1,10000)

    # line1 = Line([1+1j,3+2j])
    # point = 2+1j
    # ax.plot(line1.at(t).real, line1.at(t).imag)
    # ax.plot(point.real, point.imag,"ro")
    # ax.plot([point.real,line1.closest_point(point).real], [point.imag,line1.closest_point(point).imag],"r--")

    # ax.set_aspect("equal")
    # plt.show()













    svgPaths, _ = svgpathtools.svg2paths(r"E:\-.-\Projekte\CNC\Macha\svg\test07.svg")
    list_of_domains = create_list_of_domains(svgPaths)
    #print_domain_structure(list_of_domains)

    fig, ax = plt.subplots(figsize=(14,8))

    t=np.linspace(0,1,1000)
    print(list_of_domains[0].exterior)

    for nn in range(len(list_of_domains)):
        for curve in list_of_domains[nn].exterior:
            cr = curve.at(t).reshape(100,10)
            for n in range(100):
                plt.plot(cr[n].real, cr[n].imag, color=cm.jet(n/100))

        for inr in list_of_domains[nn].interior:
            for curve in inr:
                cr = curve.at(t).reshape(100,10)
                for n in range(100):
                    plt.plot(cr[n].real, cr[n].imag, color=cm.jet(n/100))
            

    print_domain_structure(list_of_domains)
    medial_axis_transform(list_of_domains[0], ax)

    ax.set_aspect("equal")
    plt.show()




    # fix, ax = plt.subplots(figsize=(10,4))
    # t=np.linspace(0,1,1000)
    # arc1 = Arc(1+1j, 1+1j, 2)

    # arc1a, arc1b = split_half(arc1)
    # arc1aa, arc1ab = split_half(arc1a)
    # ax.plot(arc1aa.at(t).real, arc1aa.at(t).imag)
    # ax.plot(arc1ab.at(t).real, arc1ab.at(t).imag)
    # ax.plot(arc1b.at(t).real, arc1b.at(t).imag)
    # ax.set_aspect("equal")
    # plt.show()
