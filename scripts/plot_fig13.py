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
        rows.append(dict(x=r, tau_max=a[m,4].max()*ts,
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
            tmax=max(tmax,float(np.abs(tau[bag]).max()))
        rows.append(dict(x=x, tau_max=tmax*ts,
                         tau_mean=0.5*(max(sm)-min(sm))*ts, ed_mean=max(em)*es))
    return pd.DataFrame(rows).sort_values('x')

kim_rpm=pd.read_csv(CSV_RPM,skiprows=[1]); kim_rpm['RPM']=pd.to_numeric(kim_rpm['RPM'])
kim_rpm=kim_rpm.sort_values('RPM')
kim_deg=pd.read_csv(CSV_DEG).sort_values('theta_deg')
l6=from_timeseries('fig13a_l6_mf_rpm'); l8=from_timeseries('fig13a_rm_mf_rpm')
l9=from_frames(RPM_RUNS_L9); l9d=from_frames(DEG_RUNS)

fig,(ax,ax2)=plt.subplots(1,2,figsize=(14,4.8))
def panel(a,kx,kdf,ours,xlabel,xticks):
    a2=a.twinx()
    a.plot(kdf[kx],kdf['tau_liq_max'],color='royalblue',marker='o',ms=6,lw=1.2,label=r"Kim $\tau'_{w,max}$")
    a.plot(kdf[kx],kdf['tau_liq_mean'],color='royalblue',marker='o',ms=6,lw=1.2,ls='--',mfc='w',label=r"Kim $\langle\tau'_w\rangle$")
    a2.plot(kdf[kx],kdf['Ediss_liq_max'],color='firebrick',marker='s',ms=5,lw=1.2,label=r"Kim $\epsilon'_{w,max}$")
    a2.plot(kdf[kx],kdf['Ediss_liq_mean'],color='firebrick',marker='s',ms=5,lw=1.2,ls='--',mfc='w',label=r"Kim $\langle\epsilon'_w\rangle$")
    for d,c,l in ours:
        if d.empty: continue
        a.plot(d['x'],d['tau_max'],color=c,marker='D',ms=5,lw=1.1,label=rf"$\tau'_{{w,max}}$ ({l})")
        a.plot(d['x'],d['tau_mean'],color=c,marker='D',ms=5,lw=1.1,ls='--',mfc='w',label=rf"$\langle\tau'_w\rangle$ ({l})")
        a2.plot(d['x'],d['ed_mean'],color=c,marker='^',ms=5,lw=1.1,ls=':',label=rf"$\langle\epsilon'_w\rangle$ ({l})")
    a.set_yscale('log'); a2.set_yscale('log')
    a.set_xlabel(xlabel); a.set_ylabel('Shear stress (Pa)',color='royalblue')
    a2.set_ylabel(r'EDR (W/m$^3$)',color='firebrick')
    a.tick_params(axis='y',colors='royalblue'); a2.tick_params(axis='y',colors='firebrick')
    a.tick_params(which='both',direction='in'); a2.tick_params(which='both',direction='in')
    a.grid(True,which='major',ls=':',alpha=0.4); a.set_xticks(xticks)
    return a2

a2a=panel(ax,'RPM',kim_rpm,[(l6,'darkorange','L6'),(l8,'darkred','L8'),(l9,'seagreen','L9')],
          r'Rocking frequency $f_b$ (rpm)',RPMS)
ax.set_title(r'(a)  $\theta_{b,max}=7°$',loc='left',fontsize=10)
a2b=panel(ax2,'theta_deg',kim_deg,[(l9d,'seagreen','L9')],
          r'Rocking angle $\theta_{b,max}$ (deg)',[2,3,4,5,6,7])
ax2.set_title(r'(b)  $f_b=32.5$ rpm',loc='left',fontsize=10)
h1,l1=ax.get_legend_handles_labels(); h2,l2=a2a.get_legend_handles_labels()
fig.legend(h1+h2,l1+l2,fontsize=7,loc='center left',bbox_to_anchor=(1.0,0.5),frameon=False)
fig.tight_layout(); fig.savefig(OUT,dpi=150,bbox_inches='tight')
print('saved',OUT)
