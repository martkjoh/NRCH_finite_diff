from matplotlib import animation, cm
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

rgba_to_hex = lambda rgba : '#'+''.join([f'{int(v*255):02x}' for v in rgba])
color = rgba_to_hex(cm.viridis(.25))


def add_phase(ax, phibar1, phibar2, a):
    aa = [0, 0.5, 1, 1.5]
    N = 500
    L = 1.1
    u, v = np.linspace(-L, L, N), np.linspace(-L, L, N) 
    u, v = np.meshgrid(u, v)

    f1 = lambda u, v, a :  (- 1 + 3 * (u**2 + v**2) - sqrt((3 * (u**2 - v**2))**2 - a**2 + 0j)).real
    f2 = lambda u, v, a :  (- 1 + 3 * (u**2 + v**2) + sqrt((3 * (u**2 - v**2))**2 - a**2 + 0j)).real
    f = [f1, f2]

    g = lambda v, a: sqrt(v**2 + a/3)
    v0 = np.linspace(0, L, 500)
    sgn1 = [1, 1, -1, -1]
    sgn2 = [1, -1, 1, -1]

    for i in range(2):
        u0 = g(v0, a)
        for s1, s2 in zip(sgn1, sgn2):
            ax.plot(s1*v0, s2*u0, "purple")
            ax.plot(s1*u0, s2*v0, "purple")

        eig = f[i](u, v, a)
        if i==0:
            ax.contour(u, v, eig, levels=[0], colors="black", linestyles="-")
            ax.contourf(u, v, eig, levels=[np.min(f[i](u, v, a)), 0],  colors=color, alpha=0.3)
        elif i==1:
            ax.contour(u, v, eig, levels=[0], colors="blue", linestyles="--")
            ax.contourf(u, v, eig, levels=[np.min(f[i](u, v, a)), 0],  colors=color, alpha=0.6, hatches=["///"])

        if a>0:
            lim = max(0,(1 - a)/6)
            v1 = np.linspace(np.sqrt(lim), np.sqrt(1/6), 100)
            v2 = np.sqrt(1/3 - v1**2)

            for s1, s2 in zip(sgn1, sgn2):
                ax.plot(s1*v1, s2*v2, color="green")
                ax.plot(s1*v2, s2*v1, color="green")
        
        sq = 1 / sqrt(2)
        ax.plot([sq, sq, -sq, -sq, sq], [sq, -sq, -sq, sq, sq], 'k--', alpha=.2, lw=3)

        ax.plot(phibar2, phibar1, 'ro')

        ax.set_xlabel("$v_1$")
        ax.set_ylabel("$v_2$")
        ax.set_xlim(-L, L)
        ax.set_ylim(-L, 0.1)
        ax.set_title("$a = %.1f$"%(a))


def plot_error(ax, phit, param):
    u, a, b, phibar1, phibar2, N, L, T, dt = param
    dx = L / N
    pt = np.einsum('txi->ti', phit) * dx / L
    dpt = (pt[1:] - pt[:-1])/dt
    frames = len(phit)
    t = np.linspace(0, frames*dt, frames-1)
    ax.plot(t, dpt[:,0], label="$\\frac{\\mathrm{d} \\bar \\varphi_1}{\\mathrm{d} t}$")
    ax.plot(t, dpt[:,1], label="$\\frac{\\mathrm{d} \\bar \\varphi_2}{\\mathrm{d} t}$")

    ax2 = ax.twinx()
    t = np.linspace(0, frames*dt, frames)
    ax2.plot(t, pt[:, 0]-pt[0,0], 'k--', label="$\\varphi_1(t) - \\varphi_1(0)$")
    ax2.plot(t, pt[:, 1]-pt[0,1], 'r--', label="$\\varphi_2(t) - \\varphi_2(0)$")

    ax.set_ylabel("$\\dot{\\bar\\varphi}$")
    ax.set_xlabel("$t$")
    ax2.set_ylabel("$\\Delta\\bar\\varphi$")

    ax.legend(loc=3)
    ax2.legend(loc=4)

def make_anim(folder, filename):
    filename = filename[:-4]

    phit, param = load_file(folder, filename)
    u, a, b, phibar1,  phibar2, N, L, T, dt = param
    dx = L / N
    x = np.linspace(0, L, N)

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
    ax3.set_xlabel("$\\varphi_2$")
    ax3.set_ylabel("$\\varphi_1$")

    ax = [ax1, ax2, ax3]
    fig.suptitle(", ".join(filename_from_param(param).split('_')))

    plot_error(ax4, phit, param)

    xx = np.linspace(0, L, 1000)
    ax[0].plot(xx, - np.tanh( sqrt(u / 2) * (xx - 8.2)) /np.sqrt(2) , label = "$\\tanh(\\sqrt{|r|/2} x)$", lw=4, alpha=.5)

    l1, = ax[0].plot([], [], 'r-', label='$\\varphi_1$')
    l2, = ax[0].plot([], [], 'k-', label='$\\varphi_2$')

    ax[0].plot([0, L], [phibar1, phibar1], 'r--')
    ax[0].plot([0, L], [phibar2, phibar2], 'k--')

    ax[0].set_xlim(0, L) 
    ax[0].set_ylim(-1.2, 1.2)
    ax[0].legend(loc=1)
    l5 = ax[0].text(L/10, 1, 'progress:')
    

    t = np.linspace(0, 2*pi)
    prange = 1.2
    m1, = ax[1].plot([], [], 'r--.')
    ax[1].plot(phibar2, phibar1, 'ro')
    ax[1].plot(np.cos(t), np.sin(t), 'k--') 
    ax[1].set_xlim(-prange, prange)
    ax[1].set_ylim(-prange, prange)

    add_phase(ax[2], phibar1, phibar2, a/u)

    frames = len(phit)

    n = 10
    def animate(m):
        m = m*n
        n2 = frames//10
        txt = 'progress:' + (m+1)//n2*'|'
        l5.set_text(txt)

        p = phit[m]
        l1.set_data(x, p[:, 0])
        l2.set_data(x, p[:, 1])

        m1.set_data([*p[:, 1], p[0, 1]], [*p[:, 0], p[0, 0]])

        n3 = frames//100
        if m//n3 - (m-n)//n3 == 1:
            txt = str((m+1)//n3) + "%"
            print(current_process().name, '\t', txt)

    anim = animation.FuncAnimation(fig, animate, cache_frame_data=False,   interval=1, frames=frames//n, repeat=False)
    plt.show()
    # anim.save(folder_vid+filename+".mp4", fps=30)

name = "sep3"

folder = "data/assym/" + name + "/"
folder_vid = "vid/assym/" + name + "/"

import os, shutil
from multiprocessing import Pool, current_process
if os.path.isdir(folder_vid):
    shutil.rmtree(folder_vid)
folders = folder_vid.split("/")
for i in range(len(folders)):
    fol = "/".join(folders[0:i+1]) + "/"
    if not os.path.isdir(fol):
        os.mkdir(fol)

fnames = get_all_filenames_in_folder(folder)


import time
startTime = time.time()

[make_anim(folder, fname) for fname in fnames[:]]

# folder_fname = [(folder, name) for name in fnames]
# with Pool(10) as pool:
#     pool.starmap(make_anim, folder_fname)

executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime))
