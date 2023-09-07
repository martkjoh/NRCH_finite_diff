import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.gridspec import GridSpec
import os, sys
sys.path.insert(0, os.path.abspath("./"))
from loadfiles import *

plt.rc("font", family="serif", size=16)
plt.rc("mathtext", fontset="cm")
plt.rc("lines", lw=2)


def add_phase(ax, phibar1, alpha):
    from numpy import pi, sin, cos
    from matplotlib import cm, ticker, colors as c

    rgba_to_hex = lambda rgba : '#'+''.join([f'{int(v*255):02x}' for v in rgba])
    color = rgba_to_hex(cm.viridis(.25))

    x0 = 1/np.sqrt(3)
    x1 = 1/np.sqrt(2)
    xx1 = np.linspace(x0, x1, 1000)
    xx2 = np.linspace(x1, 1, 1000)
    xx = np.linspace(x0, 1, 1000)
    xx3 = np.linspace(-x1, -x0, 1000)
    xx4 = np.linspace(-x0, x0, 1000)
    x = np.linspace(-1, 1)

    a = lambda x : np.sqrt(-3*x**4 + 4*x**2 - 1)

    ax.plot(-xx1, a(xx1), 'b-.')
    ax.plot(-xx2, a(xx2), 'k-')
    ax.plot([-x1, -x1], [a(x1), 1], 'g-')

    ax.fill_between(-xx, a(xx),np.zeros_like(xx), color=color, alpha=0.3, linewidth=0.0)
    ax.fill_between(xx3, a(xx3), 1, color=color, alpha=0.6, hatch='///', linewidth=0.0, edgecolor='#00000000')
    ax.fill_between(-xx3, a(xx3), 1, color=color, alpha=0.6, hatch='///', linewidth=0.0, edgecolor='#00000000')
    ax.fill_between(xx4, 0, 1, color=color, alpha=0.6, hatch='///', linewidth=0.0, edgecolor='#00000000')
    
    ax.plot(x, x**2, '--', color='purple', label='EL')
    ax.plot(-np.abs(phibar1), alpha, 'ro')

    ax.set_ylabel("$\\alpha/|r|$")
    ax.set_xlabel("$\\sqrt{u} \\bar \\varphi/|r|$")

    ax.set_ylim(0, .8)
    ax.set_xlim(-1, .1)


def plot_error(ax, phit, param):
    u, a, b, phibar1, phibar2, N, L, T, dt = param
    dx = L / N
    ax2 = ax.twinx()

    pt = np.einsum('txi->ti', phit) * dx / L
    dpt = (pt[1:] - pt[:-1])/dt
    frames = len(phit)
    t = np.linspace(0, frames*dt, frames-1)
    ax.plot(t, dpt[:,0], label="$\\frac{\\mathrm{d} \\bar \\varphi_1}{\\mathrm{d} t}$")
    ax.plot(t, dpt[:,1], label="$\\frac{\\mathrm{d} \\bar \\varphi_2}{\\mathrm{d} t}$")

    t = np.linspace(0, frames*dt, frames)
    ax2.plot(t, pt[:, 0]-pt[0,0], 'k--', label="$\\varphi_1(t) - \\varphi_1(0)$")
    ax2.plot(t, pt[:, 1]-pt[0,1], 'r--', label="$\\varphi_2(t) - \\varphi_2(0)$")

    ax.set_ylabel("$\\dot{\\bar\\varphi}$")
    ax.set_xlabel("$t$")
    ax2.set_ylabel("$\\Delta\\bar\\varphi$")

    legend_1 = ax.legend(loc=3)
    legend_1.remove()
    ax2.add_artist(legend_1)
    ax2.legend(loc=4)


def make_anim(folder, filename):
    filename = filename[:-4]

    phit, param = load_file(folder, filename)
    u, a, b, phibar1, phibar2, N, L, T, dt = param
    dx = L / N
    x = np.linspace(0, L, N)

    fig = plt.figure(layout="constrained", figsize=(18, 12))
    gs = GridSpec(3, 2, figure=fig)
    ax1 = fig.add_subplot(gs[0, :]) 
    ax2 = fig.add_subplot(gs[1:, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[2, 0])

    ax = [ax1, ax2, ax3]
    fig.suptitle(", ".join(filename_from_param(param).split('_')))

    plot_error(ax4, phit, param)

    l1, = ax[0].plot([], [], 'r-', label='$\\varphi_1$')
    l2, = ax[0].plot([], [], 'k-', label='$\\varphi_2$')
    ax[0].plot([0, L], [phibar1, phibar1], 'r--')
    ax[0].plot([0, L], [phibar2, phibar2], 'k--')

    ax[0].set_xlim(0, L) 
    ax[0].set_ylim(-1.2, 1.2)
    ax[0].legend(loc=1)
    l5 = ax[0].text(L/10, 1, 'progress:')
    
    t = np.linspace(0, 2*np.pi)
    prange = 1.2
    m1, = ax[1].plot([], [], 'r--.')
    ax[1].plot(0, phibar1, phibar2, 'ro')
    ax[1].plot(np.cos(t), np.sin(t), 'k--') 
    ax[1].set_xlim(-prange, prange)
    ax[1].set_ylim(-prange, prange)
    ax[1].set_xlabel("$\\varphi_2$")
    ax[1].set_ylabel("$\\varphi_1$")

    add_phase(ax[2], phibar1, a/u)

    frames = len(phit)

    n = 100
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
            print(txt)

    anim = animation.FuncAnimation(fig, animate, interval=1, frames=frames//n)
    # plt.show()
    anim.save(folder_vid+filename+".mp4", fps=30)

name = '1'
folder = "article/data/" + name + "/"
folder_vid = "article/vid/" + name + "/"
fnames = get_all_filenames_in_folder(folder)
[make_anim(folder, fname) for fname in fnames]
