"""Does the argmax(|tau|) location in OUR fork's tau field move across
time, now that we have a real cs mask (job 5083674) instead of the
analytic reconstruction? Settles whether the previously-flagged
"stationary vortex" (diary.md 2026-08-19) was a genuine flow feature or
a mask-reconstruction artifact.
"""
import glob, math
import numpy as np

N = 1024
rho_w, mu_w = 1.0e3, 1.0e-3
mu_a = 1.81e-5
L_bio, ANGLE, RPM = 0.25, 7.0, 32.5
Th_max = math.radians(ANGLE)
T_per = 60.0 / RPM
Ly = 0.03575 / L_bio
H_bio = L_bio * Ly
V_bio = L_bio/4*(H_bio + 0.5*L_bio*math.tan(Th_max))
U_bio = V_bio/(H_bio*0.5)/T_per
Re_w = rho_w*U_bio*L_bio/mu_w
mu1 = 1.0/Re_w
mu2 = (mu_a/mu_w)*mu1

def mu_of_f(f):
    fc = np.clip(f, 0, 1)
    return 1.0/(fc*(1.0/mu1 - 1.0/mu2) + 1.0/mu2)

def load(pattern):
    files = sorted(glob.glob(pattern))
    chunks = [np.loadtxt(f) for f in files]
    return np.vstack(chunks)

def to_grid(data):
    dx = 1.0/N
    ix = np.clip(np.round((data[:,0]+0.5-dx/2)/dx).astype(int), 0, N-1)
    iy = np.clip(np.round((data[:,1]+0.5-dx/2)/dx).astype(int), 0, N-1)
    ux = np.full((N,N), np.nan); uy = np.full((N,N), np.nan)
    f = np.full((N,N), np.nan); cs = np.full((N,N), np.nan)
    ux[iy,ix]=data[:,2]; uy[iy,ix]=data[:,3]; f[iy,ix]=data[:,4]; cs[iy,ix]=data[:,5]
    return ux,uy,f,cs

import re
files = glob.glob("/oscar/scratch/eaguerov/tmp/fork_l10_rampmatch/rundir/DataOurs_1024_*_*.txt")
pat = re.compile(r"DataOurs_1024_([^_]+)_\d+\.txt$")
all_times = sorted({float(pat.search(f).group(1)) for f in files if pat.search(f)})
targets = [1.19, 3.57, 5.96, 8.34, 10.72, 13.10]
TIMES = [min(all_times, key=lambda x: abs(x - g)) for g in targets]
for t in TIMES:
    data = load(f"/oscar/scratch/eaguerov/tmp/fork_l10_rampmatch/rundir/DataOurs_1024_{t:.6g}_*.txt")
    ux,uy,f,cs = to_grid(data)
    dx=1.0/N
    du_dy=np.full((N,N),np.nan); dv_dx=np.full((N,N),np.nan)
    du_dy[:,1:-1]=(ux[:,2:]-ux[:,:-2])/(2*dx)
    dv_dx[1:-1,:]=(uy[2:,:]-uy[:-2,:])/(2*dx)
    tau = mu_of_f(f)*(du_dy+dv_dx)
    mask = (f>0.5)&(cs>0.5)&~np.isnan(tau)
    tau_m = np.where(mask, np.abs(tau), np.nan)
    iy, ix = np.unravel_index(np.nanargmax(tau_m), tau_m.shape)
    x = (ix+0.5)/N - 0.5
    y = (iy+0.5)/N - 0.5
    print(f"t={t:8.4f}  argmax|tau| at (x={x:+.4f}, y={y:+.4f})  |tau|={tau_m[iy,ix]:.4g}")
