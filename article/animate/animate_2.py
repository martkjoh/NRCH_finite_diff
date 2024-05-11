from matplotlib import animation
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
from numpy import pi, sqrt

import os, sys
sys.path.insert(0, os.path.abspath("./"))
from loadfiles import *

plt.rc("font", family="serif", size=16)
plt.rc("mathtext", fontset="cm")
plt.rc("lines", lw=2)

SAVE = False

def plot_vid(anim, path, **kwargs):
    if SAVE: anim.save(path, **kwargs)
    else: plt.show()


def make_anim(folder, filenames):
    filenames = [f[:-4] for f in filenames]
    filename = filenames[0].split("*")[0]

    phitparams =  [load_file(folder, filename) for filename in filenames]
    param = phitparams[0][1]
    u, a, b, phibar1, phibar2, N, L, T, dt = param 
    phits = [phit for (phit, param) in phitparams]
    dx = L / N
    x = np.linspace(0, L, N)

    fig = plt.figure(layout="constrained", figsize=(30, 12), dpi=60)
    fig.suptitle(", ".join(filename_from_param(param).split('_')))

    gs = GridSpec(3, 4, figure=fig)
    ax = []

    for i in range(4):
        axa = fig.add_subplot(gs[0:1, i])
        axb = fig.add_subplot(gs[1:3, i])
        axa.set_xlabel("$x$")
        axa.set_ylabel("$\\varphi$")
        axb.set_xlabel("$\\varphi_2$")
        axb.set_ylabel("$\\varphi_1$")
        ax.append([axa, axb])
    
    ls = []
    ms = []
    t = np.linspace(0, 2*pi)
    prange = 1.2
    frames = len(phits[0])

    n = 100
    for i, axi in enumerate(ax):
        axa, axb = axi

        l1, = axa.plot([], [], 'r-', label='$\\varphi_1$')
        l2, = axa.plot([], [], 'k-', label='$\\varphi_2$')
        axa.plot([0, L], [phibar1, phibar1], 'r--')
        axa.plot([0, L], [0, 0], 'k--')
        axa.set_xlim(0, L) 
        axa.set_ylim(-1.2, 1.2)
        axa.legend(loc=1)

        ls.append([l1, l2])

        m, = axb.plot([], [], 'r--.')
        axb.plot(0, phibar1, 'ro')
        axb.plot(np.cos(t), np.sin(t), 'k--') 
        axb.set_xlim(-prange, prange)
        axb.set_ylim(-prange, prange)

        axa.set_xlabel("$x$")
        axa.set_ylabel("$\\varphi$")
        axb.set_xlabel("$\\varphi_2$")
        axb.set_ylabel("$\\varphi_1$")

        ms.append(m)

    def animate(k):
        k = k*n
        for i, phit in enumerate(phits):
            l1, l2 = ls[i]
            m = ms[i]
            p = phit[k]
            l1.set_data(x, p[:, 0])
            l2.set_data(x, p[:, 1])
            m.set_data([*p[:, 1], p[0, 1]], [*p[:, 0], p[0, 0]])

        n3 = frames//100
        if k//n3 - (k-n)//n3 == 1:
            txt = str((k+1)//n3) + "%"


    anim = animation.FuncAnimation(fig, animate, interval=1, frames=frames//n)
    plot_vid(anim, folder_vid+filename+".mp4", fps=30)

name = "2"
folder = "article/data/" + name + "/"
folder_vid = "article/vid/" + name + "/"
fnames = get_all_filenames_in_folder(folder)
make_anim(folder, fnames)
