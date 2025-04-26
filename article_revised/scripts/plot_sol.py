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

name = "3"
folder = "article/data/" + name + "/"
fnames = get_all_filenames_in_folder(folder)

def plot(filename):
    phit, param = load_file(folder, filename)

    u, a, b, phibar1, phibar2, N, L, T, dt = param
    dx = L / N
    x = np.linspace(0, L, N)

    fig = plt.figure(layout="constrained", figsize=(10,5))
    gs = GridSpec(2, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0]) 
    ax2 = fig.add_subplot(gs[:, 1])
    ax3 = fig.add_subplot(gs[1, 0])

    ax = [ax1, ax2, ax3]
    title = [filename_from_param(param).split('_')[i] for i in (0, 1, 3)]
    title[2] = "\\bar\\varphi" + title[2][3:]
    title[1] = "\\alpha" + title[1][1:]
    fig.suptitle("$" + ", ".join(title) + "$")

    sol1 = (1 + phibar1)*np.cos(2*x/L*2*np.pi) + phibar1
    sol2 = 2*np.sqrt(-phibar1-phibar1**2)*np.cos(x/L*2*np.pi)

    ax[0].plot(x, sol1, '-.', color='gray', alpha=1, label="sol", lw=5)
    ax[0].plot(x, sol2, '-.', color='gray', alpha=1, lw=5)
    

    T = -50

    p = phit[T]
    ax[0].plot(x, p[:, 0], 'r-', label='$\\varphi_1$')
    ax[0].plot(x, p[:, 1], 'k-', label='$\\varphi_2$')
    ax[0].plot([0, L], [phibar1, phibar1], 'r--')
    ax[0].plot([0, L], [0, 0], 'k--')

    ax[0].set_xlim(0, L) 
    ax[0].set_ylim(-1.2, 1.2)
    ax[0].set_ylabel("$\\varphi$")
    ax[0].legend(loc=1)

    d1 = p[:, 0] - sol1
    d2 = p[:, 1] - sol2
    dtot = (np.sum(np.abs(d1)) + np.sum(np.abs(d2))) * dx / L
    dmax = np.max([np.max(np.abs(d1)), np.max(np.abs(d2))])

    ax[2].plot(x, d1, 'r-', label='$\\Delta\\varphi_1$')
    ax[2].plot(x, d2, 'k-', label='$\\Delta\\varphi_2$')
    ax[2].set_xlim(0, L)
    ax[2].set_ylabel("$\\Delta\\varphi$")
    ax[2].set_xlabel("$x$")
    ax[2].legend()

    t = np.linspace(0, 2*pi)
    prange = 1.2
    ax[1].plot(0, phibar1, 'ro')
    ax[1].plot(np.cos(t), np.sin(t), 'k--') 
    ax[1].set_xlim(-prange, prange)
    ax[1].set_ylim(-prange, prange)
    ax[1].plot([*p[:, 1], p[0, 1]], [*p[:, 0], p[0, 0]], 'r--.')
    ax[1].set_xlabel("$\\varphi_2$")
    ax[1].set_ylabel("$\\varphi_1$")

    ax[2].ticklabel_format(axis='y', scilimits=(1,-1))
    plt.show()
    # plt.savefig("article/scripts/fig/sol2_"+str(i)+".pdf")


def plot_poster(filename):
    phit, param = load_file(folder, filename)

    u, a, b, phibar1, phibar2, N, L, T, dt = param
    dx = L / N
    x = np.linspace(0, L, N)

    fig, ax = plt.subplots(1, 2, figsize=(9, 4.5), sharey=True)

    sol1 = (1 + phibar1)*np.cos(2*x/L*2*np.pi) + phibar1
    sol2 = 2*np.sqrt(-phibar1-phibar1**2)*np.cos(x/L*2*np.pi)

    ax[0].plot(x/L, sol1, '-.', color='gray', alpha=1, label="sol", lw=5)
    ax[0].plot(x/L, sol2, '-.', color='gray', alpha=1, lw=5)

    p = phit[-1]
    ax[0].plot(x/L, p[:, 0], 'r-', label='$\\varphi_1$')
    ax[0].plot(x/L, p[:, 1], 'k-', label='$\\varphi_2$')
    ax[0].plot([0, 1], [phibar1, phibar1], 'r--')
    ax[0].plot([0, 1], [0, 0], 'k--')

    ax[0].set_xlim(0, 1.) 
    ax[0].set_ylim(-1.2, 1.2)
    ax[1].set_xlabel("$x/L$")
    ax[0].set_ylabel("$\\varphi/\\varphi^*$")
    ax[0].legend(loc=1)

    d1 = p[:, 0] - sol1
    d2 = p[:, 1] - sol2
    dtot = (np.sum(np.abs(d1)) + np.sum(np.abs(d2))) * dx / L
    dmax = np.max([np.max(np.abs(d1)), np.max(np.abs(d2))])


    t = np.linspace(0, 2*pi)
    prange = 1.2
    ax[1].plot(0, phibar1, 'ro')
    ax[1].plot(np.cos(t), np.sin(t), 'k--') 
    ax[1].set_xlim(-prange, prange)
    ax[1].set_ylim(-prange, prange)
    ax[1].plot([*p[:, 1], p[0, 1]], [*p[:, 0], p[0, 0]], 'r--.')
    ax[1].set_xlabel("$\\varphi_2/\\varphi^*$")

    plt.tight_layout()
    plt.show()
    # plt.savefig("fig/sol2_poster_"+str(i)+".pdf")


for i, filename in enumerate(fnames):
    filename = filename[:-4]
    plot(filename)
    plot_poster(filename) 
    