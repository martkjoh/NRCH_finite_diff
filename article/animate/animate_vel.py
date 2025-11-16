import numpy as np
from numpy import pi, sqrt
from matplotlib import animation, cm
import matplotlib.pyplot as plt
import os, sys
sys.path.insert(0, os.path.abspath("./"))
from loadfiles import *

plt.rc("font", family="serif", size=16)
plt.rc("mathtext", fontset="cm")
plt.rc("lines", lw=2)

rgba_to_hex = lambda rgba : '#'+''.join([f'{int(v*255):02x}' for v in rgba])
color = rgba_to_hex(cm.viridis(.25))

bn = 0
def get_w(f, N, x, L):
    zerp = []
    zerm = []
    for i in range(N): # 1:N
        j = (i+1)%(N)
        if f[i]*f[j]<0:
            x1, x2 = x[i], x[j]
            y1, y2 = f[i], f[j]
            if x2<x1: x2 = x2 + L 
            x0 = x1 - (x2 - x1) * y1/(y2 - y1)
            x0 = x0 % L
            if f[i]<0: zerm.append(x0) # push!(zerm, x0)
            else: zerp.append(x0) # push!(zerp, x0)

    return [*zerm, *zerp]


def get_wf(folder, filename):
    filename = filename[:-4]

    phit, param = load_file(folder, filename)
    u, a, b, phibar1, phibar2, N, L, T, dt = param
    dx = L / N
    x = np.linspace(0, L, N)

    ft = phit[:, :, 0]
    gt = phit[:, :, 1]
    wf = np.array([get_w(f, N, x, L) for f in ft]).T
    wg = np.array([get_w(g, N, x, L) for g in gt]).T

    return ft, gt, wf, wg, x, param


def get_vs(folder, filename):
    ft, gt, wf, wg, x, param = get_wf(folder, filename)
    M = len(wf[0])
    u, a, b, phibar1, phibar2, N, L, T, dt = param
    dt = T / M
    t = np.linspace(0, T, M)
    w = [wf, wg][bn]
    S = np.abs(w[1, 0] - w[0, 0])
    v = a / S * L / (L - S)
    dxdt = (w[1, 2:] - w[1, 0:-2]) / (2*dt)
    vav = np.mean(dxdt[len(dxdt)//2:]) # last half
    return a, v, vav, S

name = 'vel'

sizes = [50, 75, 100, 125, 150, 200, 250, 300,][::-1]
fig, ax = plt.subplots(1,2, sharey=True)
Ss = []
for i, size in enumerate(sizes):
    folder = "article/data/" + name + "/{size}/".format(size=size)
    fnames = get_all_filenames_in_folder(folder)
    a, v0, va, S = np.array([get_vs(folder, fname) for fname in fnames]).T
    indx = np.argsort(a)[:10]
    a = a[indx]
    v0 = v0[indx]
    va = va[indx]
    color = cm.viridis(i/(len(sizes)-1))
    ax[0].plot(a, va, '-o', color=color, label="$S = {S:.0f}$".format(S=S[0]))
    ax[1].plot(S[-1], va[-1], 'o', color=color)
    Ss.append(S[-1])

ax[0].set_xlabel('$\\alpha$')
ax[0].set_ylabel('$v$')

Ss = np.array(Ss)
L = 50
ax[1].plot(Ss, va[-1] * Ss[-1] / Ss , 'k--', zorder=0, label="$1 / S$")
ax[1].set_xlabel('$S$')

ax[0].legend(fontsize=12)
ax[1].legend()
plt.show()
