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


def plot_error(ax, phit, param):
    u, a, b, phibar1, phibar2, N, L, T, dt = param
    ax2 = ax.twinx()
    dx = L / N
    frames = len(phit)
    t = np.linspace(0, frames*dt, frames-1)

    pt = np.einsum('txi->ti', phit) * dx / L
    dpt = (pt[1:] - pt[:-1])/dt
    ax.plot(t, dpt[:,0], label="$\\frac{\\mathrm{d} \\bar \\varphi_1}{\\mathrm{d} t}$")
    ax.plot(t, dpt[:,1], label="$\\frac{\\mathrm{d} \\bar \\varphi_2}{\\mathrm{d} t}$")

    t = np.linspace(0, frames*dt, frames)
    ax2.plot(t, pt[:, 0]-pt[0,0], 'k--', label="$\\varphi_1(t) - \\varphi_1(0)$")
    ax2.plot(t, pt[:, 1]-pt[0,1], 'r--', label="$\\varphi_2(t) - \\varphi_2(0)$")
    
    ax.set_ylabel("$\\dot{\\bar\\varphi}$")
    ax.set_xlabel("$t$")
    ax2.set_ylabel("$\\Delta\\bar\\varphi$")

    legend1 = ax.legend(loc=3)
    legend1.remove()
    ax2.add_artist(legend1)

    ax2.legend(loc=4)

def plot_sol2(ax, param):
    u, a, b, phibar1, phibar2, N, L, T, dt = param
    tt = np.linspace(0, L, 1000)
    ax.plot(tt, (1 + phibar1)*np.cos(2*tt/L*2*np.pi) + phibar1,":r",label="$A\\cos2\\phi + c$")
    ax.plot(tt, 2*np.sqrt(-phibar1-phibar1**2)*np.cos(tt/L*2*np.pi),":k",label="$B\\cos^2\\phi$")


def make_anim(folder, filename):
    filename = filename[:-4]
    print(filename)

    phit, param = load_file(folder, filename)
    u, a, b, phibar1, phibar2, N, L, T, dt = param
    dx = L / N
    x = np.linspace(0, L, N)
    D2 = lambda J : ( np.roll(J, 1, axis=-1) + np.roll(J, -1, axis=-1) - 2 * J ) / (dx)**2 
    D = lambda J : (np.roll(J, 1, axis=-1) - np.roll(J, -1, axis=-1) ) / (2 * dx)

    fig = plt.figure(layout="constrained", figsize=(18, 12))
    gs = GridSpec(3, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, :]) 
    ax2 = fig.add_subplot(gs[1:, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[2, 0])

    ax1.set_xlabel("$x$")
    ax1.set_ylabel("$\\varphi$")
    ax2.set_xlabel("$\\varphi_2$")
    ax2.set_ylabel("$\\varphi_1$")
    ax3.set_xlabel("$x$")
    ax3.set_ylabel("$\\Delta \\varphi$")

    ax = [ax1, ax2, ax3]
    fig.suptitle(", ".join(filename_from_param(param).split('_')))

    plot_error(ax4, phit, param)

    plot_sol2(ax[0], param)

    l1, = ax[0].plot([], [], 'r-', label='$\\varphi_1$')
    l2, = ax[0].plot([], [], 'k-', label='$\\varphi_2$')
    ax[0].plot([0, L], [phibar1, phibar1], 'r--')
    ax[0].plot([0, L], [0, 0], 'k--')

    ax[0].set_xlim(0, L) 
    ax[0].set_ylim(-1.2, 1.2)
    ax[0].legend(loc=1)
    l5 = ax[0].text(L/10, 1, 'progress:')

    l3, = ax[2].plot([], [], 'r-', label='$\\varphi_1$')
    l4, = ax[2].plot([], [], 'k-', label='$\\varphi_2$')
    ax[2].set_xlim(0, L) 
    # ax[2].set_ylim(-0.02, 0.02)
    l6 = ax[2].text(L/10, 0.05, 'average error:\nmax error:')

    tt = np.linspace(0, L, N)
    sol1 = (1 + phibar1)*np.cos(2*tt/L*2*np.pi) + phibar1
    sol2 = 2*np.sqrt(-phibar1-phibar1**2)*np.cos(tt/L*2*np.pi)
    

    t = np.linspace(0, 2*pi)
    prange = 1.2
    m1, = ax[1].plot([], [], 'r--.')
    ax[1].plot(0, phibar1, 'ro')
    ax[1].plot(np.cos(t), np.sin(t), 'k--') 
    ax[1].set_xlim(-prange, prange)
    ax[1].set_ylim(-prange, prange)


    frames = len(phit)

    n = 1
    def animate(m):
        m = m*n
        n2 = frames//10
        txt = 'progress:' + (m+1)//n2*'|'
        l5.set_text(txt)

        p = phit[m]
        l1.set_data(x, p[:, 0])
        l2.set_data(x, p[:, 1])

        d1 = p[:, 0] - sol1
        d2 = p[:, 1] - sol2
        dtot = (np.sum(np.abs(d1)) + np.sum(np.abs(d2))) * dx / L
        dmax = np.max([np.max(np.abs(d1)), np.max(np.abs(d2))])
        l3.set_data(x, d1)
        l4.set_data(x, d2)
        txt = 'average error: ' + str(dtot) + "\nmax error: " + str(dmax)
        l6.set_text(txt)

        m1.set_data([*p[:, 1], p[0, 1]], [*p[:, 0], p[0, 0]])

        n3 = frames//1

        if m//n3 - (m-n)//n3 == 1:
            txt = str((m+1)//n3) + "%"
            print(current_process().name, '\t', txt)

    anim = animation.FuncAnimation(fig, animate, cache_frame_data=False,   interval=1, frames=frames//n, repeat=False)
    # plt.show()
    anim.save(folder_vid+filename+".mp4", fps=30)

name = "sol"
folder = "article_revised/data/" + name + "/"
folder_vid = "article_revised/vid/" + name + "/"

import os, shutil
from multiprocessing import Pool, current_process
if os.path.isdir(folder_vid):
    shutil.rmtree(folder_vid)
os.mkdir(folder_vid)

fnames = get_all_filenames_in_folder(folder)
print(fnames)

import time
startTime = time.time()

[make_anim(folder, fname) for fname in fnames]

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))