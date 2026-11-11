import os
import sys

# Ensure UTF-8 output on Windows consoles
if sys.platform == 'win32' and hasattr(sys.stdout, 'reconfigure'):
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def michaelis_menten(s, v_max, km):
    """
    Standard hyperbolic Michaelis-Menten rate equation:
    v = (V_max * [S]) / (K_m + [S])
    """
    return (v_max * s) / (km + s)

def load_kinetics_data(filepath='data/alkaline_phosphatase_kinetics.csv'):
    """
    Loads raw experimental velocity replicates across substrate concentrations.
    Computes mean initial velocity and standard deviation.
    """
    df = pd.read_csv(filepath)
    rep_cols = [c for c in df.columns if 'velocity' in c]
    df['mean_velocity'] = df[rep_cols].mean(axis=1)
    df['std_velocity'] = df[rep_cols].std(axis=1)
    return df

def fit_nonlinear_model(s_data, v_data):
    """
    Performs non-linear least-squares fitting of the Michaelis-Menten model.
    Returns estimated V_max, K_m, standard errors, and R-squared.
    """
    # Initial parameter guess: max velocity and median substrate concentration
    p0 = [np.max(v_data), np.median(s_data)]
    popt, pcov = curve_fit(michaelis_menten, s_data, v_data, p0=p0, bounds=(0, [np.inf, np.inf]))
    v_max, km = popt[0], popt[1]
    v_max_err, km_err = np.sqrt(np.diag(pcov))

    # Calculate R-squared
    residuals = v_data - michaelis_menten(s_data, v_max, km)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((v_data - np.mean(v_data))**2)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    return v_max, km, v_max_err, km_err, r_squared

def fit_lineweaver_burk(s_data, v_data):
    """
    Lineweaver-Burk double reciprocal transformation:
    1 / v = (K_m / V_max) * (1 / [S]) + (1 / V_max)
    """
    inv_s = 1.0 / s_data
    inv_v = 1.0 / v_data
    slope, intercept = np.polyfit(inv_s, inv_v, 1)

    v_max = 1.0 / intercept
    km = slope * v_max

    # R-squared in reciprocal space
    fit_inv_v = slope * inv_s + intercept
    ss_res = np.sum((inv_v - fit_inv_v)**2)
    ss_tot = np.sum((inv_v - np.mean(inv_v))**2)
    r_squared = 1.0 - (ss_res / ss_tot)

    return v_max, km, slope, intercept, r_squared

def fit_hanes_woolf(s_data, v_data):
    """
    Hanes-Woolf transformation:
    [S] / v = (1 / V_max) * [S] + (K_m / V_max)
    """
    s_over_v = s_data / v_data
    slope, intercept = np.polyfit(s_data, s_over_v, 1)

    v_max = 1.0 / slope
    km = intercept * v_max

    fit_vals = slope * s_data + intercept
    ss_res = np.sum((s_over_v - fit_vals)**2)
    ss_tot = np.sum((s_over_v - np.mean(s_over_v))**2)
    r_squared = 1.0 - (ss_res / ss_tot)

    return v_max, km, slope, intercept, r_squared

def plot_michaelis_menten(s_data, v_data, v_std, v_max, km, r2, output_path='figures/michaelis_menten_curve.png'):
    """
    Generates publication-quality hyperbolic saturation plot with V_max and K_m annotations.
    """
    s_dense = np.linspace(0, np.max(s_data) * 1.1, 300)
    v_dense = michaelis_menten(s_dense, v_max, km)

    plt.figure(figsize=(7.5, 4.8))
    plt.errorbar(s_data, v_data, yerr=v_std, fmt='o', color='#1E3A5F', ecolor='#9E472A',
                 elinewidth=1.5, capsize=4, label='Observed Data (mean ± s.d.)', markersize=6)
    plt.plot(s_dense, v_dense, color='#0284C7', linewidth=2, label=f'Nonlinear Fit (R² = {r2:.4f})')

    # Reference asymptotes
    plt.axhline(y=v_max, color='#DC2626', linestyle='--', alpha=0.7, label=f'V_max = {v_max:.2f} µM/min')
    plt.axhline(y=v_max / 2, color='#16A34A', linestyle=':', alpha=0.7, label=f'½ V_max = {v_max/2:.2f} µM/min')
    plt.axvline(x=km, color='#9333EA', linestyle=':', alpha=0.7, label=f'K_m = {km:.3f} mM')

    plt.title('Alkaline Phosphatase Kinetics: Michaelis-Menten Fit', fontsize=12, pad=12, fontweight='bold')
    plt.xlabel('Substrate Concentration [pNPP] (mM)', fontsize=10)
    plt.ylabel('Initial Velocity v₀ (µM · min⁻¹)', fontsize=10)
    plt.xlim(left=0)
    plt.ylim(bottom=0, top=v_max * 1.18)
    plt.legend(frameon=True, fontsize=8.5, loc='lower right')
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_lineweaver_burk(s_data, v_data, v_max, km, slope, intercept, r2, output_path='figures/lineweaver_burk_plot.png'):
    """
    Generates Lineweaver-Burk double reciprocal transformation plot.
    """
    inv_s = 1.0 / s_data
    inv_v = 1.0 / v_data

    # Extrapolate to negative x-intercept (-1/Km)
    x_intercept = -1.0 / km
    x_vals = np.linspace(x_intercept * 1.15, np.max(inv_s) * 1.05, 200)
    y_vals = slope * x_vals + intercept

    plt.figure(figsize=(7.5, 4.8))
    plt.scatter(inv_s, inv_v, color='#1E3A5F', s=45, zorder=5, label='Experimental Reciprocals (1/v vs 1/[S])')
    plt.plot(x_vals, y_vals, color='#D97706', linewidth=1.8, label=f'Linear Regression (R² = {r2:.4f})')

    # Axes crossings
    plt.axvline(x=0, color='black', linewidth=0.8)
    plt.axhline(y=0, color='black', linewidth=0.8)
    plt.scatter([x_intercept], [0], color='#DC2626', s=55, zorder=6, label=f'-1/K_m = {x_intercept:.2f} mM⁻¹')
    plt.scatter([0], [intercept], color='#2563EB', s=55, zorder=6, label=f'1/V_max = {intercept:.4f} (µM/min)⁻¹')

    plt.title('Lineweaver-Burk Double Reciprocal Diagnostic Plot', fontsize=12, pad=12, fontweight='bold')
    plt.xlabel('1 / [S] (mM⁻¹)', fontsize=10)
    plt.ylabel('1 / v₀ (µM · min⁻¹)⁻¹', fontsize=10)
    plt.legend(frameon=True, fontsize=8.5, loc='upper left')
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def analyze_inhibition(filepath='data/inhibition_assay_data.csv'):
    """
    Performs kinetic comparative analysis across Control, Competitive Inhibitor,
    and Non-Competitive Inhibitor treatment groups.
    """
    df = pd.read_csv(filepath)
    s = df['substrate_conc_mM'].values

    conditions = {
        'Control (Uninhibited)': df['control_velocity'].values,
        'Competitive Inhibitor (+ Pi)': df['competitive_inhibitor_velocity'].values,
        'Non-Competitive Inhibitor (+ EDTA)': df['noncompetitive_inhibitor_velocity'].values
    }

    summary = []
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    s_dense = np.linspace(0.01, 11.0, 300)
    colors = {'Control (Uninhibited)': '#0284C7',
              'Competitive Inhibitor (+ Pi)': '#EA580C',
              'Non-Competitive Inhibitor (+ EDTA)': '#16A34A'}

    for name, v in conditions.items():
        v_max, km, _, _, r2 = fit_nonlinear_model(s, v)
        summary.append({
            'Condition': name,
            'Apparent_Vmax_uM_min': round(v_max, 2),
            'Apparent_Km_mM': round(km, 3),
            'Fit_R2': round(r2, 4)
        })

        # Panel 1: Hyperbolic saturation
        axes[0].scatter(s, v, color=colors[name], s=35, alpha=0.85)
        axes[0].plot(s_dense, michaelis_menten(s_dense, v_max, km), color=colors[name], lw=1.8,
                     label=f'{name} (Km={km:.2f}, Vmax={v_max:.1f})')

        # Panel 2: Lineweaver-Burk lines
        inv_s = 1.0 / s
        inv_v = 1.0 / v
        slope, intercept = np.polyfit(inv_s, inv_v, 1)
        x_ext = np.linspace(-1.0 / km * 1.1, np.max(inv_s), 100)
        axes[1].scatter(inv_s, inv_v, color=colors[name], s=30)
        axes[1].plot(x_ext, slope * x_ext + intercept, color=colors[name], lw=1.8, label=name)

    axes[0].set_title('Substrate Saturation Curves', fontsize=11, fontweight='bold')
    axes[0].set_xlabel('[S] (mM)', fontsize=10)
    axes[0].set_ylabel('Velocity (µM/min)', fontsize=10)
    axes[0].legend(fontsize=8, frameon=True)
    axes[0].grid(True, linestyle='--', alpha=0.3)

    axes[1].axvline(0, color='black', lw=0.7)
    axes[1].axhline(0, color='black', lw=0.7)
    axes[1].set_title('Lineweaver-Burk Inhibition Diagnostic', fontsize=11, fontweight='bold')
    axes[1].set_xlabel('1 / [S] (mM⁻¹)', fontsize=10)
    axes[1].set_ylabel('1 / v (µM/min)⁻¹', fontsize=10)
    axes[1].legend(fontsize=8, frameon=True)
    axes[1].grid(True, linestyle='--', alpha=0.3)

    plt.suptitle('Enzyme Inhibition Kinetics: Mode of Inhibition Analysis', fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig('figures/inhibition_diagnostic_plot.png', dpi=300)
    plt.close()

    res_df = pd.DataFrame(summary)
    res_df.to_csv('results/inhibition_kinetics_comparison.csv', index=False)
    return res_df

def main():
    print("============================================================")
    print("ENZYME KINETICS & MICHAELIS-MENTEN ANALYSIS ENGINE")
    print("Model System: Alkaline Phosphatase (EC 3.1.3.1) + pNPP Substrate")
    print("============================================================")

    # 1. Load data
    df = load_kinetics_data('data/alkaline_phosphatase_kinetics.csv')
    s = df['substrate_conc_mM'].values
    v = df['mean_velocity'].values
    v_std = df['std_velocity'].values

    # 2. Fit models
    nl_vmax, nl_km, nl_vmax_err, nl_km_err, nl_r2 = fit_nonlinear_model(s, v)
    lb_vmax, lb_km, lb_slope, lb_intercept, lb_r2 = fit_lineweaver_burk(s, v)
    hw_vmax, hw_km, hw_slope, hw_intercept, hw_r2 = fit_hanes_woolf(s, v)

    # Enzyme concentration: 10 nM = 0.010 uM
    enzyme_conc_uM = 0.010
    kcat = (nl_vmax / enzyme_conc_uM) / 60.0  # s^-1
    kcat_km = (kcat / (nl_km * 1e-3))        # M^-1 s^-1

    print("\nNonlinear Least-Squares Michaelis-Menten Fit:")
    print(f"  V_max: {nl_vmax:.2f} ± {nl_vmax_err:.2f} µM/min")
    print(f"  K_m:   {nl_km:.3f} ± {nl_km_err:.3f} mM")
    print(f"  R²:    {nl_r2:.4f}")
    print(f"  Turnover number (k_cat):        {kcat:.2f} s⁻¹")
    print(f"  Catalytic efficiency (k_cat/Km): {kcat_km:.2e} M⁻¹ s⁻¹")

    print("\nLinear Transformation Estimates:")
    print(f"  Lineweaver-Burk : V_max = {lb_vmax:.2f} µM/min, K_m = {lb_km:.3f} mM (R² = {lb_r2:.4f})")
    print(f"  Hanes-Woolf     : V_max = {hw_vmax:.2f} µM/min, K_m = {hw_km:.3f} mM (R² = {hw_r2:.4f})")

    # Save summary table
    summary_data = [
        {'Method': 'Nonlinear Regression (curve_fit)', 'Vmax_uM_min': nl_vmax, 'Km_mM': nl_km, 'R_squared': nl_r2},
        {'Method': 'Lineweaver-Burk Double Reciprocal', 'Vmax_uM_min': lb_vmax, 'Km_mM': lb_km, 'R_squared': lb_r2},
        {'Method': 'Hanes-Woolf Linearization', 'Vmax_uM_min': hw_vmax, 'Km_mM': hw_km, 'R_squared': hw_r2}
    ]
    pd.DataFrame(summary_data).to_csv('results/kinetic_constants_summary.csv', index=False)

    # 3. Generate diagnostic plots
    plot_michaelis_menten(s, v, v_std, nl_vmax, nl_km, nl_r2, 'figures/michaelis_menten_curve.png')
    plot_lineweaver_burk(s, v, lb_vmax, lb_km, lb_slope, lb_intercept, lb_r2, 'figures/lineweaver_burk_plot.png')

    # 4. Inhibition analysis
    print("\nAnalyzing Inhibition Assay Data...")
    inhibition_results = analyze_inhibition('data/inhibition_assay_data.csv')
    print(inhibition_results.to_string(index=False))

    print("\nAll figures generated in 'figures/' and summary CSVs saved in 'results/'.")

if __name__ == '__main__':
    main()
