import glob, math
import numpy as np

N = 1024
T = 12.1466  # a settled snapshot
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

def load(pattern, ncols):
    files = sorted(glob.glob(pattern))
    chunks = []
    for fp in files:
        try:
            arr = np.loadtxt(fp)
        except ValueError:
            arr = np.loadtxt(fp, skiprows=1)
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        chunks.append(arr[:, :ncols])
    return np.vstack(chunks)

def to_grid(data, cols):
    dx = 1.0/N
    ix = np.clip(np.round((data[:,0]+0.5-dx/2)/dx).astype(int), 0, N-1)
    iy = np.clip(np.round((data[:,1]+0.5-dx/2)/dx).astype(int), 0, N-1)
    g = {}
    for name,k in cols.items():
        a = np.full((N,N), np.nan)
        a[iy,ix] = data[:,k]
        g[name]=a
    return g

def fields(pattern, ncols, cols, is_up):
    data = load(pattern, ncols)
    g = to_grid(data, cols)
    dx = 1.0/N
    du_dy = np.full((N,N), np.nan); dv_dx = np.full((N,N), np.nan)
    du_dy[:,1:-1] = (g['ux'][:,2:]-g['ux'][:,:-2])/(2*dx)
    dv_dx[1:-1,:] = (g['uy'][2:,:]-g['uy'][:-2,:])/(2*dx)
    tau = mu_of_f(g['f'])*(du_dy+dv_dx)
    mask = (g['f']>0.5)&(g['cs']>0.5)&~np.isnan(tau)
    return tau, mask

# find nearest matching times
import re
def nearest(directory, prefix, ncols, ndecpat):
    files = glob.glob(f"{directory}/{prefix}_1024_*_*.txt")
    pat = re.compile(rf"{prefix}_1024_([^_]+)_\d+\.txt$")
    times = sorted({float(pat.search(f).group(1)) for f in files if pat.search(f)})
    j = min(times, key=lambda x: abs(x-T))
    return j

t_ours = nearest("/oscar/scratch/eaguerov/tmp/fork_l10_rampmatch/rundir", "DataOurs", 6, None)
t_up = nearest("/oscar/scratch/eaguerov/tmp/upstream_l10_video/rundir/Data_all", "Data_all", 12, None)
print("t_ours", t_ours, "t_up", t_up)

tau1, m1 = fields(f"/oscar/scratch/eaguerov/tmp/fork_l10_rampmatch/rundir/DataOurs_1024_{t_ours:g}_*.txt", 6, dict(ux=2,uy=3,f=4,cs=5), False)
tau2, m2 = fields(f"/oscar/scratch/eaguerov/tmp/upstream_l10_video/rundir/Data_all/Data_all_1024_{t_up:.12g}_*.txt", 7, dict(ux=2,uy=3,f=4,cs=6), True)
valid = m1 & m2
a, b = tau1[valid], tau2[valid]
mag = np.abs(b)
print("overall corr", np.corrcoef(a,b)[0,1], "sign_agree", np.mean(np.sign(a)==np.sign(b)))
for pct in [0, 25, 50, 75, 90, 95, 99]:
    thresh = np.percentile(mag, pct)
    sel = mag >= thresh
    print(f"|tau_up| >= p{pct} ({thresh:.3g}): n={sel.sum()}, sign_agree={np.mean(np.sign(a[sel])==np.sign(b[sel]))*100:.1f}%")
