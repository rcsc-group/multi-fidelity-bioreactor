"""Replica of Kim et al. (2024) Fig. 13, in Kim's own layout.

Kim's caption: "...shown for (a) different rocking frequencies at
theta_max=7 deg and (b) different rocking ANGLES at f_b=32.5 rpm. The left
axis represents the shear stress, and the right axis represents the energy
dissipation rate." Four series per panel: the absolute maximum over a period
(solid) and the maximum of the spatially averaged value (hollow), for both
shear stress and EDR.

[2026-09-15] Supersedes scripts/plot_fig13a_current.py, which plotted an
invented "EDR vs rpm" second panel rather than Kim's angle sweep, and whose
<tau> series used mean|tau| where Kim's tau_liq_mean is the amplitude of the
SIGNED spatial mean.

All of our values are post-mask-fix (diary.md 2026-09-15 (2)): the tau/EDR
spatial averages must be taken over liquid INSIDE the bag, not over every
cell with f>0.5 (which reaches outside the embedded geometry and diluted
every mean by ~3.5x).

Provenance of our series differs by level and is deliberate:
  L6, L8 -- per-timestep KPIs from the mask-fixed reruns (full time
            resolution). shear_stress.dat is NON-dimensional; scaled here by
            rho*U^2 and rho*U^3/L.
  L9     -- recomputed from saved field frames, because a34fc4d4 predates the
            mask fix and its own shear_stress.dat is therefore invalid.
            NOTE these frames are on the old, period-commensurate cadence, so
            L9's peaks are slight LOWER BOUNDS (diary.md 2026-09-15 (3)).

Kim's panel (b) data was digitised from Figures/Fig_tau_Ediss_rpm_deg.pdf
into csv_raw/shear_ediss_vs_angle.csv; the digitiser reproduces the existing
rpm CSV to 0.1% median error.

Usage:  uv run python scripts/plot_fig13.py
"""
import json, math
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).parent.parent
CSV_RPM = ROOT/'experiments/kimetal2024/csv_raw/shear_ediss_vs_frequency.csv'
CSV_DEG = ROOT/'experiments/kimetal2024/csv_raw/shear_ediss_vs_angle.csv'
OUT = ROOT/'experiments/kimetal2024/figure_replicas/replicated_Fig13.png'
RPMS = [17.5,20,22.5,25,27.5,30,32.5,35,37.5]
# L9 angle sweep at 32.5 rpm (frames present for all six)
DEG_RUNS = {2.0:'7e103866',3.0:'52e007f3',4.0:'81f191ed',5.0:'fcc0fc84',6.0:'fbe93e9a',7.0:'a34fc4d4'}
RPM_RUNS_L9 = {20.0:'20a14369',22.5:'255e0b87',25.0:'bce29aa5',27.5:'9ec56180',
               30.0:'60d09a80',32.5:'a34fc4d4',35.0:'a281a16f',37.5:'0f0ad3ea'}

def scales(run):
    p=json.load(open(ROOT/'runs'/run/'params.json'))
    L=p['geometry']['a']; b=p['geometry']['b']/L
    tm=p['theta_max']; th=math.radians(tm[0] if isinstance(tm,list) else tm)
    om=p['omega_b']; Tp=2*math.pi/om; H=2*L*b
    V=L/4*(H+0.5*L*math.tan(th)); U=V/(H*0.5)/Tp
    return b, 1000*U**2, 1000*U**3/L

def from_timeseries(prefix):
    """L6/L8: per-timestep KPIs from the mask-fixed reruns."""
    rows=[]
    for r in RPMS:
        run=f'{prefix}{r:g}'
        f=ROOT/'runs'/run/'shear_stress.dat'
        if not f.exists(): continue
        a=np.loadtxt(f,skiprows=1); _,ts,es=scales(run)
        t=a[:,1]; m=t>=t[0]+0.6*(t[-1]-t[0])          # settled tail
        # col 8 = tau_100_signed. Kim takes max() of the SIGNED field
        # (bio_stress.m:583-585, no abs()); col 4 is max|tau|, a different
        # statistic. Measured difference here: 0-3.5%.
        rows.append(dict(x=r, tau_max=a[m,8].max()*ts,
                         tau_mean=0.5*(a[m,10].max()-a[m,10].min())*ts,
                         ed_mean=a[m,9].max()*es))
    return pd.DataFrame(rows).sort_values('x')

def from_frames(runs):
    """L9: recomputed from fields with the correct bag mask."""
    import struct
    def load(p):
        with open(p,'rb') as fh:
            (n,)=struct.unpack('i',fh.read(4)); (t,)=struct.unpack('d',fh.read(8)); fh.read(16)
            f=np.frombuffer(fh.read(n*n*4),dtype=np.float32).reshape(n,n)
            tau=np.frombuffer(fh.read(n*n*4),dtype=np.float32).reshape(n,n)
            ed=np.frombuffer(fh.read(n*n*4),dtype=np.float32).reshape(n,n)
        return f,tau,ed
    rows=[]
    for x,run in sorted(runs.items()):
        b,ts,es=scales(run)
        fl=sorted((ROOT/'runs'/run/'frames_tau').glob('frame_*.bin'))
        if not fl: continue
        fl=fl[len(fl)//2:]
        sm=[];em=[];tmax=0.0
        for fp in fl:
            f,tau,ed=load(fp); n=f.shape[0]
            yy=(-0.5+(np.arange(n)+0.5)/n)[:,None]*np.ones((1,n))
            bag=(f>0.5)&(np.abs(yy)<b)
            if not bag.any(): continue
            sm.append(tau[bag].mean()); em.append(ed[bag].mean())
            tmax=max(tmax,float(tau[bag].max()))   # signed, as Kim
        rows.append(dict(x=x, tau_max=tmax*ts,
                         tau_mean=0.5*(max(sm)-min(sm))*ts, ed_mean=max(em)*es))
    return pd.DataFrame(rows).sort_values('x')

kim_rpm=pd.read_csv(CSV_RPM,skiprows=[1]); kim_rpm['RPM']=pd.to_numeric(kim_rpm['RPM'])
kim_rpm=kim_rpm.sort_values('RPM')
kim_deg=pd.read_csv(CSV_DEG).sort_values('theta_deg')
l6=from_timeseries('fig13a_l6_mf_rpm'); l8=from_timeseries('fig13a_rm_mf_rpm')
l9=from_frames(RPM_RUNS_L9); l9d=from_frames(DEG_RUNS)
# L10 -- Kim's OWN mesh (n_L=2^10). Cross-level sweep, each point warm-started
# from the converged L9 state at the same rpm, on BioReactor-mpi-xlevel-chain
# (built 2026-09-15 10:41, i.e. WITH the bag-mask fix).
#
# These supersede the two old l10_fig13a_*_32rank_calib points, which must NOT
# be plotted: they ran on BioReactor-mpi-video-fixed (2026-09-08), predating
# the mask fix, and it shows -- their tau amplitude is 0.32x Kim, the ~3x
# signature of the broken mask, against 1.07x for the cross-level run at the
# same condition. The diary also records them as never reaching quasi-steady.
#
# Settling verified rather than assumed: per-cycle tau amplitude over the 5
# usable cycles varies 0.8% (32.5 rpm), 5.7% (25), 2.5% (37.5) with no drift,
# so the warm start lands quasi-steady from cycle 0.
l10=from_timeseries('l10_fig13a_xlevel_rpm')

# ── encoding ──────────────────────────────────────────────────────────────
# Each visual channel carries exactly ONE meaning, and the fill convention is
# Kim's own (caption of his Fig. 13: "maximum shear stress (solid blue
# circles) and energy dissipation rate (solid red squares) ... along with the
# maximum of the spatially averaged shear stress (hollow blue circles) and
# energy dissipation rate (hollow red squares)"):
#
#   COLOUR  = dataset        Kim / L6 / L8 / L9
#   MARKER  = quantity       circle = tau      square = epsilon
#   FILL    = statistic      filled = absolute max      hollow = max of the
#                            SPATIAL MEAN
#
# Previously colour encoded dataset AND quantity at once, which put our L8
# (darkred) next to Kim's EDR (firebrick) -- two different datasets in
# near-identical red -- and our <eps> was a FILLED triangle while every other
# spatially-averaged series was hollow. Axis labels are neutral now: with
# colour meaning dataset, colouring them blue/red would point at nothing.
#
# Okabe-Ito, distinguishable in common colour-vision deficiencies.
C_KIM, C_L6, C_L8, C_L9, C_L10 = ("black", "#E69F00", "#0072B2",
                                  "#009E73", "#CC79A7")
MK_TAU, MK_EPS = "o", "s"


def _series(axis, x, y, colour, marker, filled, lw=1.2, ms=5.5):
    axis.plot(x, y, color=colour, marker=marker, ms=ms, lw=lw,
              ls="-" if filled else "--",
              mfc=colour if filled else "w", mec=colour, mew=1.1)


fig,(ax,ax2)=plt.subplots(1,2,figsize=(13.5,4.8))


def panel(a, kx, kdf, ours, xlabel, xticks):
    a2 = a.twinx()
    # Kim: all four series
    _series(a,  kdf[kx], kdf["tau_liq_max"],    C_KIM, MK_TAU, True)
    _series(a,  kdf[kx], kdf["tau_liq_mean"],   C_KIM, MK_TAU, False)
    _series(a2, kdf[kx], kdf["Ediss_liq_max"],  C_KIM, MK_EPS, True)
    _series(a2, kdf[kx], kdf["Ediss_liq_mean"], C_KIM, MK_EPS, False)
    # ours: tau max + tau mean, and EDR mean only -- these runs predate the
    # ediss_kim_max column, so the filled square is Kim's alone. Shown as the
    # absence of a series rather than by substituting a different statistic.
    for d, colour in ours:
        if d.empty:
            continue
        _series(a,  d["x"], d["tau_max"],  colour, MK_TAU, True)
        _series(a,  d["x"], d["tau_mean"], colour, MK_TAU, False)
        _series(a2, d["x"], d["ed_mean"],  colour, MK_EPS, False)
    a.set_yscale("log"); a2.set_yscale("log")
    a.set_xlabel(xlabel, fontsize=11)
    a.set_ylabel(r"Shear stress $\tau'_w$ (Pa)   $\bullet$ circles", fontsize=11)
    a2.set_ylabel(r"EDR $\epsilon'_w$ (W/m$^3$)   $\blacksquare$ squares", fontsize=11)
    a.tick_params(which="both", direction="in")
    a2.tick_params(which="both", direction="in")
    a.grid(True, which="major", ls=":", alpha=0.35)
    a.set_xticks(xticks)
    return a2


panel(ax, "RPM", kim_rpm, [(l6, C_L6), (l8, C_L8), (l9, C_L9), (l10, C_L10)],
      r"Rocking frequency $f_b$ (rpm)", RPMS)
ax.set_title(r"(a)  $\theta_{b,max}=7°$", loc="left", fontsize=10)
panel(ax2, "theta_deg", kim_deg, [(l9d, C_L9)],
      r"Rocking angle $\theta_{b,max}$ (deg)", [2, 3, 4, 5, 6, 7])
ax2.set_title(r"(b)  $f_b=32.5$ rpm", loc="left", fontsize=10)

# Two small legends instead of one list of eleven: one decodes colour
# (dataset), the other decodes marker+fill (quantity, statistic). Both sit
# outside the axes.
from matplotlib.lines import Line2D
ds = [Line2D([], [], color=c, lw=1.4, label=l) for c, l in
      [(C_KIM, "Kim et al."), (C_L6, "L6"), (C_L8, "L8"), (C_L9, "L9"),
       (C_L10, "L10  (Kim's mesh)")]]
enc = [Line2D([], [], color="0.3", marker=MK_TAU, ls="-",  mfc="0.3", ms=6,
              label=r"$\tau'_{w,max}$   (absolute max)"),
       Line2D([], [], color="0.3", marker=MK_TAU, ls="--", mfc="w",   ms=6,
              label=r"$\langle\tau'_w\rangle$   (max of spatial mean)"),
       Line2D([], [], color="0.3", marker=MK_EPS, ls="-",  mfc="0.3", ms=6,
              label=r"$\epsilon'_{w,max}$   (absolute max)"),
       Line2D([], [], color="0.3", marker=MK_EPS, ls="--", mfc="w",   ms=6,
              label=r"$\langle\epsilon'_w\rangle$   (max of spatial mean)")]
leg1 = fig.legend(handles=ds, fontsize=8.5, loc="upper left",
                  bbox_to_anchor=(1.0, 0.93), frameon=False, title="dataset")
leg2 = fig.legend(handles=enc, fontsize=8.5, loc="upper left",
                  bbox_to_anchor=(1.0, 0.63), frameon=False, title="marker / fill")
for lg in (leg1, leg2):
    lg.get_title().set_fontsize(8.5)
fig.add_artist(leg1)

fig.tight_layout()
fig.savefig(OUT, dpi=150, bbox_inches="tight")
print("saved", OUT)
