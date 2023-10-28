import matplotlib.pyplot as plt
import numpy as np
    

def P2R(radii, angles):
    return radii * np.exp(1j*angles)

def pole(a,b,c,d):
    return (b and (not d)) or (a and b and (not c) and d) or ((not a) and c and (not b) and (not d))

def plot_bb(tbb, ax=plt.subplots()[1]):
    ax.plot([tbb.real[0], tbb.real[1], tbb.real[1], tbb.real[0], tbb.real[0]], \
            [tbb.imag[0], tbb.imag[0], tbb.imag[1], tbb.imag[1], tbb.imag[0]], \
            "g-", linewidth=1)