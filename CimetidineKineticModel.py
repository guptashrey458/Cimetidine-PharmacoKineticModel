# cimetidine.py

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class CimetidineKineticModel:
    def __init__(self, logFC_dict=None):
        """Initialize PBPK model with optional logFC adjustments."""
        
        # Dose & absorption
        self.Dose = 400                # mg dose
        self.F = 0.506                 # Bioavailability (fraction absorbed)
        self.ka = 1.54                 # Absorption rate constant (1/hr)
        
        # Volumes (L)
        self.V_gut = 1.0               # Gut volume
        self.V_plasma = 3.0            # Plasma volume
        self.V_liver = 1.8             # Liver volume

        # Hepatic blood flow
        self.Q_hepatic = 90.0          # L/hr
        
        # Default clearances (L/hr)
        self.CL_OCT1 = 5.0             # Hepatic uptake
        self.CL_OCT3 = 2.0             # Additional hepatic uptake
        self.CL_MATE1 = 2.0            # Hepatic efflux
        self.CL_Pgp = 1.5              # Efflux transporter
        
        self.CL_CYP1A2 = 3.0           # Metabolic enzyme
        self.CL_CYP2C9 = 2.0           # Metabolic enzyme
        self.CL_CYP2D6 = 1.0           # Metabolic enzyme
        
        self.CL_OCT2 = 1.0             # Renal uptake
        self.CL_renal_passive = 3.0    # Passive renal clearance
        self.CL_biliary = 0.5          # Biliary clearance

        # Monte Carlo variability parameters (standard deviations)
        self.ka_sigma = 0.2
        self.CL_sigma = 0.3
        self.F_sigma = 0.05
        self.Q_hepatic_sigma = 5.0

        # Apply logFC adjustments using 2^logFC formula
        if logFC_dict:
            self.CL_OCT1 *= 2**logFC_dict.get('OCT1', 0.0)
            self.CL_OCT3 *= 2**logFC_dict.get('OCT3', 0.0)
            self.CL_MATE1 *= 2**logFC_dict.get('MATE1', 0.0)
            self.CL_Pgp *= 2**logFC_dict.get('Pgp', 0.0)
            self.CL_CYP1A2 *= 2**logFC_dict.get('CYP1A2', 0.0)
            self.CL_CYP2C9 *= 2**logFC_dict.get('CYP2C9', 0.0)
            self.CL_CYP2D6 *= 2**logFC_dict.get('CYP2D6', 0.0)

    def odes(self, y, t):
        """Final corrected PBPK ODE system"""
        A_gut, A_plasma, A_liver, A_eliminated = y

        # Concentrations
        C_gut = A_gut / self.V_gut
        C_plasma = A_plasma / self.V_plasma
        C_liver = A_liver / self.V_liver

        # Total metabolic clearance in liver
        total_metabolic_CL = self.CL_CYP1A2 + self.CL_CYP2C9 + self.CL_CYP2D6

        # Gut: simple first-order absorption only
        dA_gut = -self.ka * A_gut

        # Plasma: gain from gut, loss to liver, renal, etc.
        dA_plasma = (
            self.F * self.ka * A_gut
            - self.Q_hepatic * (C_plasma - C_liver)
            - (self.CL_OCT1 + self.CL_OCT3) * C_plasma  # hepatic uptake
            + self.CL_MATE1 * C_liver                  # hepatic efflux back
            - (self.CL_OCT2 + self.CL_renal_passive) * C_plasma  # renal
        )

        # Liver: uptake, metabolism, biliary/P-gp efflux
        dA_liver = (
            self.Q_hepatic * (C_plasma - C_liver)
            + (self.CL_OCT1 + self.CL_OCT3) * C_plasma
            - self.CL_MATE1 * C_liver
            - total_metabolic_CL * C_liver
            - (self.CL_biliary + self.CL_Pgp) * C_liver  # both biliary
        )

        # Eliminated: all pathways out of the system
        dA_eliminated = (
            (self.CL_OCT2 + self.CL_renal_passive) * C_plasma  # renal
            + total_metabolic_CL * C_liver                     # metabolic
            + (self.CL_biliary + self.CL_Pgp) * C_liver        # biliary
        )

        return [dA_gut, dA_plasma, dA_liver, dA_eliminated]

    def simulate(self, 
                 ka=None, F=None,
                 CL_OCT1=None, CL_CYP1A2=None, CL_CYP2C9=None, CL_Pgp=None,
                 duration=12, points=1000):
        """
        Simulate PBPK model with optional parameter overrides.
        """
        # Store original values
        original_params = {
            'ka': self.ka,
            'F': self.F,
            'CL_OCT1': self.CL_OCT1,
            'CL_CYP1A2': self.CL_CYP1A2,
            'CL_CYP2C9': self.CL_CYP2C9,
            'CL_Pgp': self.CL_Pgp
        }
        
        # Update with passed values if provided
        if ka is not None:
            self.ka = ka
        if F is not None:
            self.F = F
        if CL_OCT1 is not None:
            self.CL_OCT1 = CL_OCT1
        if CL_CYP1A2 is not None:
            self.CL_CYP1A2 = CL_CYP1A2
        if CL_CYP2C9 is not None:
            self.CL_CYP2C9 = CL_CYP2C9
        if CL_Pgp is not None:
            self.CL_Pgp = CL_Pgp

        # Run simulation
        t = np.linspace(0, duration, points)
        y0 = [self.Dose, 0, 0, 0]  # All drug in gut initially

        sol = odeint(self.odes, y0, t)
        A_plasma = sol[:, 1]
        C_plasma = A_plasma / self.V_plasma

        # Restore original values
        for param, value in original_params.items():
            setattr(self, param, value)

        return t, C_plasma

    def plot_model(self, duration=24, points=1000, label='Baseline'):
        """Plot plasma concentration-time profile"""
        t, C_plasma = self.simulate(duration=duration, points=points)
        plt.plot(t, C_plasma, label=label)
        plt.xlabel("Time (hr)")
        plt.ylabel("Plasma Concentration (mg/L)")
        plt.title(f"Cimetidine PBPK Model — Dose: {self.Dose} mg")
        plt.legend()
        plt.grid(True)
        plt.show()

    def knockout(self, transporters=['OCT1']):
        """Perform transporter/enzyme knockout(s) and plot results"""
        t, C_baseline = self.simulate()
        baseline_auc = np.trapz(C_baseline, t)
        baseline_peak = np.max(C_baseline)

        knockout_logFC = {t: -10 for t in transporters}  # Use -10 for complete knockout
        model_ko = CimetidineKineticModel(logFC_dict=knockout_logFC)
        t, C_ko = model_ko.simulate()

        auc_ko = np.trapz(C_ko, t)
        peak_ko = np.max(C_ko)

        auc_ratio = float(auc_ko) / float(baseline_auc)
        peak_ratio = float(peak_ko) / float(baseline_peak)

        plt.plot(t, C_baseline, label="Baseline")
        plt.plot(t, C_ko, label=f"Knockout: {','.join(transporters)}")
        plt.xlabel("Time (hr)")
        plt.ylabel("Plasma Concentration (mg/L)")
        plt.title(f"Knockout — AUC Ratio: {auc_ratio:.2f}, Cmax Ratio: {peak_ratio:.2f}")
        plt.legend()
        plt.grid(True)
        plt.show()

        return {'auc_ratio': auc_ratio, 'cmax_ratio': peak_ratio}
    def calculate_pk_metrics(self):
        t, C_plasma = self.simulate()
        auc = np.trapz(C_plasma, t)
        cmax = np.max(C_plasma)
        tmax = t[np.argmax(C_plasma)]
        half_life = self.calculate_half_life_robust(t, C_plasma)
        clearance = float(self.Dose * self.F) / float(auc)

        return {
            'AUC': auc,
            'Cmax': cmax,
            'Tmax': tmax,
            'Half_life': half_life,
            'CL/F': clearance
        }

    def calculate_half_life_robust(self, t, C_plasma):
        try:
            # Use more data points for better fit (up to 100 points)
            n_points = min(100, len(C_plasma) // 2)
            
            # Get terminal phase data
            terminal_t = t[-n_points:]
            terminal_C = C_plasma[-n_points:]
            
            # Filter out non-positive concentrations and very small values
            positive_mask = terminal_C > 1e-10  # Avoid log of very small numbers
            
            if np.sum(positive_mask) < 10:  # Need at least 10 points
                return np.nan
            
            # Use only valid concentrations
            valid_t = terminal_t[positive_mask]
            valid_C = terminal_C[positive_mask]
            
            # Check if concentrations are generally decreasing
            if len(valid_C) < 10 or valid_C[-1] >= valid_C[0]:
                return np.nan
            
            # Calculate slope of log-linear terminal phase
            log_C = np.log(valid_C)
            slope, intercept = np.polyfit(valid_t, log_C, 1)
            
            # Calculate half-life only if slope is negative (elimination)
            if slope < 0:
                half_life = np.log(2) / abs(slope)
                
                # Sanity check: half-life should be reasonable (0.1 to 100 hours)
                if 0.1 <= half_life <= 100:
                    return half_life
            
            return np.nan
            
        except Exception as e:
            print(f"Warning: Could not calculate half-life: {e}")
            return np.nan

    def simulate_monte_carlo(self, n_sim=1000, duration=24, points=1000, logFC_dict=None, seed=None):
        """
        Perform Monte Carlo simulation with parameter variability
        and return mean & 95% CI for plasma concentration.
        
        Parameters:
        -----------
        n_sim : int
            Number of Monte Carlo simulations
        duration : float
            Simulation duration in hours
        points : int
            Number of time points
        logFC_dict : dict, optional
            logFC adjustments for transporters/enzymes
        seed : int, optional
            Random seed for reproducibility
        """
        if seed is not None:
            np.random.seed(seed)
            
        t = np.linspace(0, duration, points)
        all_C_plasma = []

        for i in range(n_sim):
            # Sample log-normal variability for PK parameters using exposed sigmas
            ka = np.random.lognormal(mean=np.log(self.ka), sigma=self.ka_sigma)
            F = np.random.normal(loc=self.F, scale=self.F_sigma)
            Q_hepatic = np.random.normal(loc=self.Q_hepatic, scale=self.Q_hepatic_sigma)
            
            # Sample variability for clearances (log-normal for positive values)
            CL_OCT1 = np.random.lognormal(mean=np.log(self.CL_OCT1), sigma=self.CL_sigma)
            CL_OCT3 = np.random.lognormal(mean=np.log(self.CL_OCT3), sigma=self.CL_sigma)
            CL_MATE1 = np.random.lognormal(mean=np.log(self.CL_MATE1), sigma=self.CL_sigma)
            CL_Pgp = np.random.lognormal(mean=np.log(self.CL_Pgp), sigma=self.CL_sigma)
            CL_CYP1A2 = np.random.lognormal(mean=np.log(self.CL_CYP1A2), sigma=self.CL_sigma)
            CL_CYP2C9 = np.random.lognormal(mean=np.log(self.CL_CYP2C9), sigma=self.CL_sigma)
            CL_CYP2D6 = np.random.lognormal(mean=np.log(self.CL_CYP2D6), sigma=self.CL_sigma)

            # Store original parameters
            original_params = {
                'ka': self.ka, 'F': self.F, 'Q_hepatic': self.Q_hepatic,
                'CL_OCT1': self.CL_OCT1, 'CL_OCT3': self.CL_OCT3,
                'CL_MATE1': self.CL_MATE1, 'CL_Pgp': self.CL_Pgp,
                'CL_CYP1A2': self.CL_CYP1A2, 'CL_CYP2C9': self.CL_CYP2C9,
                'CL_CYP2D6': self.CL_CYP2D6
            }

            # Update with sampled parameters
            self.ka = ka
            self.F = np.clip(F, 0, 1)  # keep F in [0, 1]
            self.Q_hepatic = max(0, Q_hepatic)
            self.CL_OCT1 = CL_OCT1
            self.CL_OCT3 = CL_OCT3
            self.CL_MATE1 = CL_MATE1
            self.CL_Pgp = CL_Pgp
            self.CL_CYP1A2 = CL_CYP1A2
            self.CL_CYP2C9 = CL_CYP2C9
            self.CL_CYP2D6 = CL_CYP2D6

            # Apply logFC adjustments if provided
            if logFC_dict:
                self.CL_OCT1 *= 2**logFC_dict.get('OCT1', 0.0)
                self.CL_OCT3 *= 2**logFC_dict.get('OCT3', 0.0)
                self.CL_MATE1 *= 2**logFC_dict.get('MATE1', 0.0)
                self.CL_Pgp *= 2**logFC_dict.get('Pgp', 0.0)
                self.CL_CYP1A2 *= 2**logFC_dict.get('CYP1A2', 0.0)
                self.CL_CYP2C9 *= 2**logFC_dict.get('CYP2C9', 0.0)
                self.CL_CYP2D6 *= 2**logFC_dict.get('CYP2D6', 0.0)

            # Run simulation
            sol = odeint(self.odes, [self.Dose, 0, 0, 0], t)
            C_plasma = sol[:, 1].astype(float) / float(self.V_plasma)
            all_C_plasma.append(C_plasma)

            # Reset parameters for next iteration
            for param, value in original_params.items():
                setattr(self, param, value)

        all_C_plasma = np.array(all_C_plasma)
        mean_C = np.mean(all_C_plasma, axis=0)
        lower_CI = np.percentile(all_C_plasma, 2.5, axis=0)
        upper_CI = np.percentile(all_C_plasma, 97.5, axis=0)

        return t, mean_C, lower_CI, upper_CI

    def plot_monte_carlo_comparison(self, baseline_logFC=None, knockout_logFC=None, 
                                  n_sim=500, duration=24, seed=None):
        """
        Plot Monte Carlo comparison between baseline and knockout conditions
        
        Parameters:
        -----------
        baseline_logFC : dict, optional
            logFC adjustments for baseline condition
        knockout_logFC : dict, optional
            logFC adjustments for knockout condition
        n_sim : int
            Number of Monte Carlo simulations
        duration : float
            Simulation duration in hours
        seed : int, optional
            Random seed for reproducibility
        """
        # Run baseline Monte Carlo
        t_baseline, mean_baseline, lower_baseline, upper_baseline = self.simulate_monte_carlo(
            n_sim=n_sim, duration=duration, logFC_dict=baseline_logFC, seed=seed
        )
        
        # Run knockout Monte Carlo (use different seed if same seed provided)
        knockout_seed = seed + 1 if seed is not None else None
        t_ko, mean_ko, lower_ko, upper_ko = self.simulate_monte_carlo(
            n_sim=n_sim, duration=duration, logFC_dict=knockout_logFC, seed=knockout_seed
        )
        
        # Plot comparison
        plt.figure(figsize=(12, 8))
        
        # Baseline
        plt.plot(t_baseline, mean_baseline, 'b-', linewidth=2, label='Baseline (Mean)')
        plt.fill_between(t_baseline, lower_baseline, upper_baseline, 
                        alpha=0.2, color='blue', label='Baseline (95% CI)')
        
        # Knockout
        plt.plot(t_ko, mean_ko, 'r-', linewidth=2, label='Knockout (Mean)')
        plt.fill_between(t_ko, lower_ko, upper_ko, 
                        alpha=0.2, color='red', label='Knockout (95% CI)')
        
        plt.xlabel("Time (hours)", fontsize=12)
        plt.ylabel("Plasma Concentration (mg/L)", fontsize=12)
        plt.title("Monte Carlo Comparison: Baseline vs Knockout", fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
        
        # Calculate and print summary statistics
        auc_baseline = np.trapz(mean_baseline, t_baseline)
        auc_ko = np.trapz(mean_ko, t_ko)
        cmax_baseline = np.max(mean_baseline)
        cmax_ko = np.max(mean_ko)
        
        print(f"📊 Monte Carlo Comparison Summary:")
        print(f"   Baseline AUC: {auc_baseline:.2f} mg·hr/L")
        print(f"   Knockout AUC: {auc_ko:.2f} mg·hr/L")
        print(f"   AUC Ratio: {float(auc_ko)/float(auc_baseline):.2f}")
        print(f"   Baseline Cmax: {cmax_baseline:.2f} mg/L")
        print(f"   Knockout Cmax: {cmax_ko:.2f} mg/L")
        print(f"   Cmax Ratio: {float(cmax_ko)/float(cmax_baseline):.2f}")
        
        return {
            'baseline': (t_baseline, mean_baseline, lower_baseline, upper_baseline),
            'knockout': (t_ko, mean_ko, lower_ko, upper_ko)
        }

    def plot_monte_carlo(self, n_sim=1000, duration=24, logFC_dict=None, label='Monte Carlo', seed=None):
        """Plot Monte Carlo simulation with confidence intervals"""
        t, mean_C, lower_CI, upper_CI = self.simulate_monte_carlo(
            n_sim=n_sim, duration=duration, logFC_dict=logFC_dict, seed=seed
        )

        plt.figure(figsize=(10, 6))
        plt.plot(t, mean_C, label=f'{label} (Mean)', color='blue', linewidth=2)
        plt.fill_between(t, lower_CI, upper_CI, color='blue', alpha=0.2, 
                        label=f'{label} (95% CI)')
        plt.xlabel("Time (hours)", fontsize=12)
        plt.ylabel("Plasma Concentration (mg/L)", fontsize=12)
        plt.title("Cimetidine PBPK Monte Carlo Simulation", fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

        return t, mean_C, lower_CI, upper_CI

    def plot_gene_regulation_impact(self):
        """
        Plot the impact of gene regulation (logFC × fitted CL) for each gene/transporter.
        """
        import matplotlib.pyplot as plt

        # ✅ Averaged logFCs you provided
        logFC_dict = {
            'Pgp': 0.08207506625,
            'CYP1A2': -0.2293257425,
            'CYP2C9': -0.1893353167,
            'CYP2D6': 0.016903843,
            'OCT1': -0.011315867,
            'OCT3': 0.0860044,
            'MATE1': -0.0164414245
        }

        # ✅ Fitted clearances from self
        fitted_CLs = {
            'Pgp': self.CL_Pgp,
            'CYP1A2': self.CL_CYP1A2,
            'CYP2C9': self.CL_CYP2C9,
            'CYP2D6': self.CL_CYP2D6,
            'OCT1': self.CL_OCT1,
            'OCT3': self.CL_OCT3,
            'MATE1': self.CL_MATE1
        }

        # ✅ Calculate impact: logFC × fitted CL
        impact_values = {gene: logFC_dict[gene] * fitted_CLs[gene] for gene in logFC_dict}

        print("✅ Gene regulation impact (logFC × fitted CL):")
        for gene, impact in impact_values.items():
            print(f"   {gene}: {impact:.3f}")

        # ✅ Bar plot with color coding
        genes = list(impact_values.keys())
        values = [impact_values[g] for g in genes]
        colors = ['green' if v >= 0 else 'red' for v in values]

        plt.figure(figsize=(10, 6))
        bars = plt.bar(genes, values, color=colors)

        for bar, val in zip(bars, values):
            height = bar.get_height()
            plt.text(
                bar.get_x() + bar.get_width() / 2,
                height,
                f'{val:.3f}',
                ha='center',
                va='bottom' if height >= 0 else 'top'
            )

        plt.axhline(0, color='black', linewidth=0.8)
        plt.title("Gene Regulation Impact (logFC × CL) in Microgravity", fontsize=14)
        plt.ylabel("Impact Value (logFC × CL)")
        plt.xlabel("Gene / Transporter")
        plt.tight_layout()
        plt.show()

    def print_parameter_comparison(self, logFC_dict):
        """
        Print a neat table comparing Normal vs Microgravity parameter values
        and log2FC for each gene/transporter.
        """
        # Define fitted normal clearances (baseline)
        normal_CLs = {
            'Pgp': 1.5,
            'CYP1A2': 3.0,
            'CYP2C9': 2.0,
            'CYP2D6': 1.0,
            'OCT1': 5.0,
            'OCT3': 2.0,
            'MATE1': 2.0
        }

        # Microgravity-adjusted: Normal × 2^logFC
        microgravity_CLs = {}
        for gene, normal_val in normal_CLs.items():
            logFC = logFC_dict.get(gene, 0.0)
            microgravity_CLs[gene] = normal_val * (2 ** logFC)

        # Print header
        print(f"{'Gene':<10} {'Normal':>10} {'Microgravity':>15} {'log2FC':>10}")
        print("-" * 45)

        for gene in normal_CLs:
            normal_val = normal_CLs[gene]
            micro_val = microgravity_CLs[gene]
            logFC = logFC_dict.get(gene, 0.0)
            print(f"{gene:<10} {normal_val:>10.4f} {micro_val:>15.4f} {logFC:>10.4f}")

    def plot_parameter_comparison(self, logFC_dict):
        """
        Bar chart: Normal vs Microgravity parameter values for each gene/transporter.
        """
        import numpy as np
        import matplotlib.pyplot as plt

        # ✅ Baseline fitted clearances
        normal_CLs = {
            'Pgp': 1.5,
            'CYP1A2': 3.0,
            'CYP2C9': 2.0,
            'CYP2D6': 1.0,
            'OCT1': 5.0,
            'OCT3': 2.0,
            'MATE1': 2.0
        }

        # ✅ Microgravity-adjusted: Normal × 2^logFC
        microgravity_CLs = {
            gene: normal_CLs[gene] * (2 ** logFC_dict.get(gene, 0.0))
            for gene in normal_CLs
        }

        genes = list(normal_CLs.keys())
        normal_values = [normal_CLs[g] for g in genes]
        micro_values = [microgravity_CLs[g] for g in genes]

        x = np.arange(len(genes))
        width = 0.35

        plt.figure(figsize=(12, 6))
        plt.bar(x - width/2, normal_values, width, label='Normal (Baseline)', color='gray')
        plt.bar(x + width/2, micro_values, width, label='Microgravity', color='blue')

        plt.xticks(x, genes)
        plt.ylabel("Parameter Value (Clearance, L/hr)")
        plt.title("Gene Parameter Value: Normal vs Microgravity")
        plt.legend()
        plt.tight_layout()
        plt.show()
