# Enzyme Kinetics & Michaelis-Menten Analysis Engine

A computational biochemistry toolkit for analyzing enzyme catalytic velocity, fitting nonlinear Michaelis-Menten kinetics, comparing double-reciprocal linear transformations, and classifying enzyme inhibition mechanisms.

## Overview

Enzymes are biological catalysts that dramatically accelerate biochemical reactions without being consumed. Understanding their reaction kinetics, substrate affinity ($K_m$), maximum turnover rates ($V_{\max}$, $k_{\text{cat}}$), and regulatory inhibition profiles is fundamental across biochemistry, molecular pharmacology, and bioengineering.

This project investigates:
1. **Steady-State Saturation Kinetics**: Fitting nonlinear hyperbolic rate equations to experimental initial velocity ($v_0$) data across varying substrate concentrations ($[S]$).
2. **Linear Transformation Diagnostics**: Evaluating Lineweaver-Burk and Hanes-Woolf transformations to benchmark parameter estimations against nonlinear least-squares regression.
3. **Turnover & Catalytic Efficiency**: Deriving the catalytic rate constant ($k_{\text{cat}}$) and second-order specificity constant ($k_{\text{cat}} / K_m$).
4. **Inhibition Mechanism Classification**: Quantifying the distinct kinetic shifts induced by Competitive and Non-Competitive inhibitors.

### Model System: Alkaline Phosphatase (ALP)

* **Enzyme**: Calf Intestinal / Bacterial Alkaline Phosphatase (EC 3.1.3.1)
* **Substrate**: $p$-Nitrophenyl Phosphate (pNPP)
* **Product**: $p$-Nitrophenol (chromogenic yellow product, measured spectrophotometrically at $\lambda = 405\text{ nm}$) + Inorganic Phosphate ($P_i$)
* **Inhibitors**:
  * *Competitive Inhibitor*: Inorganic Phosphate ($P_i$, product analog binding catalytic active site)
  * *Non-Competitive Inhibitor*: EDTA / Chelators (disrupting essential active site $\text{Zn}^{2+}/\text{Mg}^{2+}$ cofactors)

---

## Mathematical Formulation

### 1. Michaelis-Menten Rate Equation
Under the Briggs-Haldane steady-state approximation ($d[ES]/dt = 0$):

$$v_0 = \frac{V_{\max} [S]}{K_m + [S]}$$

Where:
* $v_0$ = Initial reaction velocity ($\mu\text{M} \cdot \text{min}^{-1}$)
* $V_{\max}$ = Maximum theoretical velocity at saturating substrate concentration
* $[S]$ = Substrate concentration ($\text{mM}$)
* $K_m$ = Michaelis constant ($\text{mM}$), the substrate concentration at which $v_0 = \frac{1}{2} V_{\max}$

### 2. Catalytic Rate & Efficiency
$$k_{\text{cat}} = \frac{V_{\max}}{[E]_0} \quad (\text{s}^{-1})$$

$$\text{Catalytic Efficiency} = \frac{k_{\text{cat}}}{K_m} \quad (\text{M}^{-1}\text{s}^{-1})$$

### 3. Lineweaver-Burk Double Reciprocal Transformation
$$\frac{1}{v_0} = \frac{K_m}{V_{\max}} \cdot \frac{1}{[S]} + \frac{1}{V_{\max}}$$

* **y-intercept**: $\frac{1}{V_{\max}}$
* **x-intercept**: $-\frac{1}{K_m}$
* **Slope**: $\frac{K_m}{V_{\max}}$

### 4. Diagnostic Criteria for Enzyme Inhibition

| Inhibition Mode | Apparent $K_m$ ($K_m^{\text{app}}$) | Apparent $V_{\max}$ ($V_{\max}^{\text{app}}$) | Lineweaver-Burk Visual Pattern |
| :--- | :--- | :--- | :--- |
| **Competitive** | Increases ($\alpha K_m$) | Unchanged ($V_{\max}$) | Common intersection at y-axis ($1/V_{\max}$) |
| **Non-Competitive** | Unchanged ($K_m$) | Decreases ($V_{\max}/\alpha'$) | Common intersection at negative x-axis ($-1/K_m$) |
| **Uncompetitive** | Decreases | Decreases | Parallel lines (identical slope $K_m/V_{\max}$) |

---

## Project Structure

```text
enzyme-kinetics-analysis/
|-- data/
|   |-- alkaline_phosphatase_kinetics.csv
|   |-- inhibition_assay_data.csv
|-- figures/
|   |-- michaelis_menten_curve.png
|   |-- lineweaver_burk_plot.png
|   |-- inhibition_diagnostic_plot.png
|-- results/
|   |-- kinetic_constants_summary.csv
|   |-- inhibition_kinetics_comparison.csv
|-- analyze_kinetics.py
|-- explain.txt
|-- requirements.txt
|-- .gitignore
|-- README.md
```

---

## Getting Started

### Prerequisites
* Python 3.8 or higher
* pip package manager

### Installation
```bash
cd enzyme-kinetics-analysis
pip install -r requirements.txt
```

### Running the Analysis
```bash
python analyze_kinetics.py
```

---

## Experimental Results Summary

### 1. Kinetic Parameter Estimates (Alkaline Phosphatase)
* **Substrate**: pNPP ($0.05 - 10.0\text{ mM}$)
* **Enzyme Concentration $[E]_0$**: $10\text{ nM}$

| Method | Estimated $V_{\max}$ ($\mu\text{M/min}$) | Estimated $K_m$ ($\text{mM}$) | Goodness of Fit ($R^2$) |
| :--- | :--- | :--- | :--- |
| **Nonlinear Least-Squares (`curve_fit`)** | **$52.88 \pm 0.44$** | **$0.468 \pm 0.011$** | **$0.9998$** |
| **Lineweaver-Burk Transformation** | $51.72$ | $0.442$ | $0.9989$ |
| **Hanes-Woolf Linearization** | $52.79$ | $0.465$ | $0.9999$ |

* **Turnover Number ($k_{\text{cat}}$)**: $88.13\text{ s}^{-1}$
* **Catalytic Efficiency ($k_{\text{cat}} / K_m$)**: $1.88 \times 10^5\text{ M}^{-1}\text{s}^{-1}$

### 2. Inhibition Assay Profiling

| Condition | Apparent $V_{\max}$ ($\mu\text{M/min}$) | Apparent $K_m$ ($\text{mM}$) | Kinetic Diagnosis |
| :--- | :--- | :--- | :--- |
| **Control (Uninhibited)** | $52.88$ | $0.468$ | Baseline catalytic activity |
| **Competitive Inhibitor ($+ P_i$)** | $52.48$ | $1.264$ | $K_m$ shifted $\approx 2.7\times$, $V_{\max}$ preserved |
| **Non-Competitive Inhibitor ($+ \text{EDTA}$)** | $26.43$ | $0.465$ | $V_{\max}$ reduced by $50\%$, $K_m$ invariant |

---

