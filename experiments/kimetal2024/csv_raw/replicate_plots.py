import os
import pandas as pd
import matplotlib.pyplot as plt

# Reconstruct the exact Matlab style aesthetics
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 10,
    'axes.labelsize': 11,
    'xtick.labelsize': 9.5,
    'ytick.labelsize': 9.5,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.bbox': 'tight',
    
    # Matlab style tick markers
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': False, # Will handle manually for twinx
    'axes.grid': False # Grid handles manually
})

# Matlab's exact default color palette (RGB values from [0,1] range)
MATLAB_BLUE = '#0072BD'
MATLAB_ORANGE = '#D95319'
MATLAB_GREEN = '#77AC30'
MATLAB_PURPLE = '#7E2F8E'

script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = script_dir
output_dir = os.path.abspath(os.path.join(script_dir, "..", "figure_replicas"))
os.makedirs(output_dir, exist_ok=True)

# Setup the standard axes box layout in pixels (margins)
def create_axes_layout(fig, is_double=False):
    # Standard margins and sizes mimicking Matlab layout
    if not is_double:
        # left, bottom, width, height (proportions of figure)
        return fig.add_axes([0.13, 0.15, 0.74, 0.76])
    else:
        # For double panels
        ax_left = fig.add_axes([0.07, 0.15, 0.39, 0.76])
        ax_right = fig.add_axes([0.54, 0.15, 0.39, 0.76])
        return ax_left, ax_right

def style_axes(ax1, ax2, x_label, y1_label, y2_label, is_angle=False):
    # Set labels
    ax1.set_xlabel(x_label, labelpad=6)
    ax1.set_ylabel(y1_label, labelpad=6)
    ax2.set_ylabel(y2_label, labelpad=8)
    
    # Configure grid lines (dotted, light gray)
    ax1.grid(True, which='both', linestyle=':', color='#c8c8c8', linewidth=0.5)
    
    # Configure tick lines pointing inwards on all sides
    ax1.tick_params(direction='in', top=True, right=False)
    ax2.tick_params(direction='in', top=True, left=False, right=True)
    
    # X axis ticks and limits
    if is_angle:
        ax1.set_xlim([1.7, 7.3])
        ax1.set_xticks([2, 3, 4, 5, 6, 7])
    else:
        ax1.set_xlim([14.0, 38.5])
        ax1.set_xticks([15, 20, 25, 30, 35])

# ==============================================================================
# Figure 9: Mixing Time vs Frequency (7 deg)
# ==============================================================================
def plot_figure_9():
    csv_path = os.path.join(data_dir, "mixing_kla_vs_frequency.csv")
    df = pd.read_csv(csv_path)
    
    # Exactly 4.5" x 3.15" figure (676x472 px at 150 dpi)
    fig = plt.figure(figsize=(4.5, 3.15))
    ax1 = create_axes_layout(fig)
    ax2 = ax1.twinx()
    
    # Primary axis: Mixing time (s)
    h1 = ax1.plot(df["RPM"], df["dtmix_strict_0.95"], 'o-', color=MATLAB_BLUE, 
                  markerfacecolor=MATLAB_BLUE, markeredgecolor=MATLAB_BLUE, 
                  markersize=5, linewidth=1.2, label=r'$\chi = 0.95$')
    h2 = ax1.plot(df["RPM"], df["dtmix_strict_0.75"], '^-', color=MATLAB_ORANGE, 
                  markerfacecolor=MATLAB_ORANGE, markeredgecolor=MATLAB_ORANGE, 
                  markersize=5, linewidth=1.2, label=r'$\chi = 0.75$')
    h3 = ax1.plot(df["RPM"], df["dtmix_strict_0.5"], 'v-', color=MATLAB_GREEN, 
                  markerfacecolor=MATLAB_GREEN, markeredgecolor=MATLAB_GREEN, 
                  markersize=5, linewidth=1.2, label=r'$\chi = 0.50$')
    
    # Secondary axis: vorticity magnitude (1/s)
    h4 = ax2.plot(df["RPM"], df["vor_meanabs_steady_streaming"], 's--', color=MATLAB_PURPLE, 
                  markerfacecolor='none', markeredgecolor=MATLAB_PURPLE, 
                  markersize=5, linewidth=1.2, label=r'$\langle |\overline{\xi}_b^\prime| \rangle$')
    
    style_axes(ax1, ax2, 
               x_label=r"$f_b$ (rpm)", 
               y1_label=r"$\Delta t_{\mathrm{mix}}$ (s)", 
               y2_label=r"$\langle |\overline{\xi}_b^\prime| \rangle$ ($\mathrm{s}^{-1}$)")
    
    ax1.set_ylim([0, 750])
    ax2.set_ylim([0, 2.0])
    
    # Combine legend with neat borders
    lns = h1 + h2 + h3 + h4
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc='upper right', frameon=True, edgecolor='black', 
               facecolor='white', framealpha=1.0, fancybox=False)
    
    plt.savefig(os.path.join(output_dir, "replicated_Fig9.png"), dpi=150)
    plt.close()
    print("Replicated Figure 9 saved.")

# ==============================================================================
# Figure 10: Mixing Time vs Angle (32.5 rpm)
# ==============================================================================
def plot_figure_10():
    csv_path = os.path.join(data_dir, "mixing_kla_vs_angle.csv")
    df = pd.read_csv(csv_path)
    
    fig = plt.figure(figsize=(4.5, 3.15))
    ax1 = create_axes_layout(fig)
    ax2 = ax1.twinx()
    
    h1 = ax1.plot(df["Angle_deg"], df["dtmix_strict_0.5"], 'v-', color=MATLAB_GREEN, 
                  markerfacecolor=MATLAB_GREEN, markeredgecolor=MATLAB_GREEN, 
                  markersize=5, linewidth=1.2, label=r'$\chi = 0.50$')
    
    h2 = ax2.plot(df["Angle_deg"], df["vor_meanabs_steady_streaming"], 's--', color=MATLAB_PURPLE, 
                  markerfacecolor='none', markeredgecolor=MATLAB_PURPLE, 
                  markersize=5, linewidth=1.2, label=r'$\langle |\overline{\xi}_b^\prime| \rangle$')
    
    style_axes(ax1, ax2, 
               x_label=r"$\theta_{b,\mathrm{max}}$ (deg)", 
               y1_label=r"$\Delta t_{\mathrm{mix}}$ (s)", 
               y2_label=r"$\langle |\overline{\xi}_b^\prime| \rangle$ ($\mathrm{s}^{-1}$)",
               is_angle=True)
    
    ax1.set_ylim([0, 400])
    ax2.set_ylim([0, 1.5])
    
    lns = h1 + h2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc='upper right', frameon=True, edgecolor='black', 
               facecolor='white', framealpha=1.0, fancybox=False)
    
    plt.savefig(os.path.join(output_dir, "replicated_Fig10.png"), dpi=150)
    plt.close()
    print("Replicated Figure 10 saved.")

# ==============================================================================
# Figure 11: kLa vs Frequency (7 deg)
# ==============================================================================
def plot_figure_11():
    csv_path = os.path.join(data_dir, "mixing_kla_vs_frequency.csv")
    df = pd.read_csv(csv_path)
    
    fig = plt.figure(figsize=(4.5, 3.15))
    ax1 = create_axes_layout(fig)
    ax2 = ax1.twinx()
    
    h1 = ax1.plot(df["RPM"], df["kLa_exp5pts_50"], 'o-', color=MATLAB_BLUE, 
                  markerfacecolor=MATLAB_BLUE, markeredgecolor=MATLAB_BLUE, 
                  markersize=5, linewidth=1.2, label=r'$C_{w,\mathrm{oxy}}^* = 0.50$')
    h2 = ax1.plot(df["RPM"], df["kLa_exp5pts_25"], '^-', color=MATLAB_ORANGE, 
                  markerfacecolor=MATLAB_ORANGE, markeredgecolor=MATLAB_ORANGE, 
                  markersize=5, linewidth=1.2, label=r'$C_{w,\mathrm{oxy}}^* = 0.25$')
    h3 = ax1.plot(df["RPM"], df["kLa_exp5pts_10"], 'v-', color=MATLAB_GREEN, 
                  markerfacecolor=MATLAB_GREEN, markeredgecolor=MATLAB_GREEN, 
                  markersize=5, linewidth=1.2, label=r'$C_{w,\mathrm{oxy}}^* = 0.10$')
    
    h4 = ax2.plot(df["RPM"], df["vor_meanabs_steady_streaming"], 's--', color=MATLAB_PURPLE, 
                  markerfacecolor='none', markeredgecolor=MATLAB_PURPLE, 
                  markersize=5, linewidth=1.2, label=r'$\langle |\overline{\xi}_b^\prime| \rangle$')
    
    style_axes(ax1, ax2, 
               x_label=r"$f_b$ (rpm)", 
               y1_label=r"$k_La$ ($\mathrm{h}^{-1}$)", 
               y2_label=r"$\langle |\overline{\xi}_b^\prime| \rangle$ ($\mathrm{s}^{-1}$)")
    
    ax1.set_ylim([0, 60])
    ax2.set_ylim([0, 2.0])
    
    lns = h1 + h2 + h3 + h4
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc='upper left', frameon=True, edgecolor='black', 
               facecolor='white', framealpha=1.0, fancybox=False)
    
    plt.savefig(os.path.join(output_dir, "replicated_Fig11.png"), dpi=150)
    plt.close()
    print("Replicated Figure 11 saved.")

# ==============================================================================
# Figure 12: kLa vs Angle (32.5 rpm)
# ==============================================================================
def plot_figure_12():
    csv_path = os.path.join(data_dir, "mixing_kla_vs_angle.csv")
    df = pd.read_csv(csv_path)
    
    fig = plt.figure(figsize=(4.5, 3.15))
    ax1 = create_axes_layout(fig)
    ax2 = ax1.twinx()
    
    h1 = ax1.plot(df["Angle_deg"], df["kLa_exp5pts_25"], '^-', color=MATLAB_ORANGE, 
                  markerfacecolor=MATLAB_ORANGE, markeredgecolor=MATLAB_ORANGE, 
                  markersize=5, linewidth=1.2, label=r'$C_{w,\mathrm{oxy}}^* = 0.25$')
    h2 = ax1.plot(df["Angle_deg"], df["kLa_exp5pts_10"], 'v-', color=MATLAB_GREEN, 
                  markerfacecolor=MATLAB_GREEN, markeredgecolor=MATLAB_GREEN, 
                  markersize=5, linewidth=1.2, label=r'$C_{w,\mathrm{oxy}}^* = 0.10$')
    
    h3 = ax2.plot(df["Angle_deg"], df["vor_meanabs_steady_streaming"], 's--', color=MATLAB_PURPLE, 
                  markerfacecolor='none', markeredgecolor=MATLAB_PURPLE, 
                  markersize=5, linewidth=1.2, label=r'$\langle |\overline{\xi}_b^\prime| \rangle$')
    
    style_axes(ax1, ax2, 
               x_label=r"$\theta_{b,\mathrm{max}}$ (deg)", 
               y1_label=r"$k_La$ ($\mathrm{h}^{-1}$)", 
               y2_label=r"$\langle |\overline{\xi}_b^\prime| \rangle$ ($\mathrm{s}^{-1}$)",
               is_angle=True)
    
    ax1.set_ylim([0, 25])
    ax2.set_ylim([0, 1.5])
    
    lns = h1 + h2 + h3
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc='upper left', frameon=True, edgecolor='black', 
               facecolor='white', framealpha=1.0, fancybox=False)
    
    plt.savefig(os.path.join(output_dir, "replicated_Fig12.png"), dpi=150)
    plt.close()
    print("Replicated Figure 12 saved.")

# ==============================================================================
# Figure 13: Shear Stress & EDR (Double panel)
# ==============================================================================
def plot_figure_13():
    csv_a = os.path.join(data_dir, "shear_ediss_vs_frequency.csv")
    df_a = pd.read_csv(csv_a)
    
    csv_b = os.path.join(data_dir, "shear_ediss_vs_angle.csv")
    df_b = pd.read_csv(csv_b)
    
    # 9.6" x 3.15" double panel figure (1440x472 px at 150 dpi)
    fig = plt.figure(figsize=(9.6, 3.15))
    ax_a, ax_b = create_axes_layout(fig, is_double=True)
    
    # PANEL A: Frequency sweep
    ax_a_r = ax_a.twinx()
    
    # Left Axis: Shear stress (blue)
    h1 = ax_a.plot(df_a["RPM"], df_a["tau_liq_max"], 'o-', color=MATLAB_BLUE, 
                  markerfacecolor=MATLAB_BLUE, markeredgecolor=MATLAB_BLUE, 
                  markersize=5, linewidth=1.2, label=r'$\tau_{w,\mathrm{max}}^\prime$ (abs max)')
    h2 = ax_a.plot(df_a["RPM"], df_a["tau_liq_mean"], 'o--', color=MATLAB_BLUE, 
                  markerfacecolor='none', markeredgecolor=MATLAB_BLUE, 
                  markersize=5, linewidth=1.2, label=r'$\langle \tau_w^\prime \rangle$ (sp averaged max)')
    
    # Right Axis: EDR (red)
    h3 = ax_a_r.plot(df_a["RPM"], df_a["Ediss_liq_max"], 's-', color=MATLAB_ORANGE, 
                  markerfacecolor=MATLAB_ORANGE, markeredgecolor=MATLAB_ORANGE, 
                  markersize=5, linewidth=1.2, label=r'$\epsilon_{w,\mathrm{max}}^\prime$ (abs max)')
    h4 = ax_a_r.plot(df_a["RPM"], df_a["Ediss_liq_mean"], 's--', color=MATLAB_ORANGE, 
                  markerfacecolor='none', markeredgecolor=MATLAB_ORANGE, 
                  markersize=5, linewidth=1.2, label=r'$\langle \epsilon_w^\prime \rangle$ (sp averaged max)')
    
    style_axes(ax_a, ax_a_r, 
               x_label=r"$f_b$ (rpm)", 
               y1_label=r"$\tau_{w}^\prime$ (Pa)", 
               y2_label=r"$\epsilon_{w}^\prime$ ($\mathrm{W/m}^3$)")
    
    ax_a.set_ylim([0, 1.4])
    ax_a_r.set_ylim([0, 2000])
    
    lns_a = h1 + h2 + h3 + h4
    labs_a = [l.get_label() for l in lns_a]
    ax_a.legend(lns_a, labs_a, loc='upper left', frameon=True, edgecolor='black', 
               facecolor='white', framealpha=1.0, fancybox=False)
    ax_a.text(0.02, 0.93, '(a)', transform=ax_a.transAxes, fontweight='bold')
    
    # PANEL B: Angle sweep
    ax_b_r = ax_b.twinx()
    
    h5 = ax_b.plot(df_b["Angle_deg"], df_b["tau_liq_max"], 'o-', color=MATLAB_BLUE, 
                  markerfacecolor=MATLAB_BLUE, markeredgecolor=MATLAB_BLUE, 
                  markersize=5, linewidth=1.2, label=r'$\tau_{w,\mathrm{max}}^\prime$ (abs max)')
    h6 = ax_b.plot(df_b["Angle_deg"], df_b["tau_liq_mean"], 'o--', color=MATLAB_BLUE, 
                  markerfacecolor='none', markeredgecolor=MATLAB_BLUE, 
                  markersize=5, linewidth=1.2, label=r'$\langle \tau_w^\prime \rangle$ (sp averaged max)')
    
    h7 = ax_b_r.plot(df_b["Angle_deg"], df_b["Ediss_liq_max"], 's-', color=MATLAB_ORANGE, 
                  markerfacecolor=MATLAB_ORANGE, markeredgecolor=MATLAB_ORANGE, 
                  markersize=5, linewidth=1.2, label=r'$\epsilon_{w,\mathrm{max}}^\prime$ (abs max)')
    h8 = ax_b_r.plot(df_b["Angle_deg"], df_b["Ediss_liq_mean"], 's--', color=MATLAB_ORANGE, 
                  markerfacecolor='none', markeredgecolor=MATLAB_ORANGE, 
                  markersize=5, linewidth=1.2, label=r'$\langle \epsilon_w^\prime \rangle$ (sp averaged max)')
    
    style_axes(ax_b, ax_b_r, 
               x_label=r"$\theta_{b,\mathrm{max}}$ (deg)", 
               y1_label=r"$\tau_{w}^\prime$ (Pa)", 
               y2_label=r"$\epsilon_{w}^\prime$ ($\mathrm{W/m}^3$)",
               is_angle=True)
    
    ax_b.set_ylim([0, 0.3])
    ax_b_r.set_ylim([0, 60])
    
    lns_b = h5 + h6 + h7 + h8
    labs_b = [l.get_label() for l in lns_b]
    ax_b.legend(lns_b, labs_b, loc='upper left', frameon=True, edgecolor='black', 
               facecolor='white', framealpha=1.0, fancybox=False)
    ax_b.text(0.02, 0.93, '(b)', transform=ax_b.transAxes, fontweight='bold')
    
    plt.savefig(os.path.join(output_dir, "replicated_Fig13.png"), dpi=150)
    plt.close()
    print("Replicated Figure 13 saved.")

if __name__ == "__main__":
    plot_figure_9()
    plot_figure_10()
    plot_figure_11()
    plot_figure_12()
    plot_figure_13()
    print("All replicated plots generated successfully!")
