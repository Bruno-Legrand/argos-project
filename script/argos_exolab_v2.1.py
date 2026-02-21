import matplotlib.pyplot as plt
import lightkurve as lk
import numpy as np
import pandas as pd
from astroquery.mast import Catalogs
import os
from datetime import datetime

# --- SYSTEM CONFIGURATION ---
script_dir = os.path.dirname(os.path.abspath(__file__))
base_path = os.path.abspath(os.path.join(script_dir, ".."))
export_dir = os.path.join(base_path, "exports")

# Modular namespacing to support future ARGOS branches (Flarelab, Novalab, etc.)
history_path = os.path.join(export_dir, "argos_exolab_targets.txt")
log_path = os.path.join(export_dir, "argos_exolab_observations.log")
csv_path = os.path.join(export_dir, "argos_exolab_database.csv")

if not os.path.exists(export_dir): 
    os.makedirs(export_dir)

def estimate_spectral_type(teff):
    """Classify host star based on Effective Temperature."""
    if teff < 3800: return "M"
    elif teff < 5300: return "K"
    elif teff < 6000: return "G"
    else: return "F"

# --- TARGET LIST & BENCHMARK SELECTION ---
# Replace with your 50-target list for full production runs
tic_list = [261155555, 261259521, 234520440]

print(f"🚀 ARGOS v2.1 - FAULT-TOLERANT DETECTION ENGINE")

processed = []
if os.path.exists(history_path):
    with open(history_path, "r") as f: 
        processed = f.read().splitlines()

for tic in tic_list:
    target = f"TIC {tic}"
    if str(tic) in processed: 
        continue

    try:
        # --- BLOC A: STELLAR SOURCING (MAST CATALOGS) ---
        star_info = Catalogs.query_object(target, radius=0.001, catalog="TIC")
        r_star = star_info['rad'][0] if not np.isnan(star_info['rad'][0]) else 1.0
        t_eff = star_info['Teff'][0] if not np.isnan(star_info['Teff'][0]) else 5778.0
        mag = star_info['Tmag'][0]
        ra, dec = star_info['ra'][0], star_info['dec'][0]
        type_spectral = estimate_spectral_type(t_eff)
        
        # --- BLOC B: PIPELINE LOGIC (SINGLE-SECTOR STABLE) ---
        search = lk.search_lightcurve(target, author="SPOC")
        if len(search) == 0: 
            search = lk.search_lightcurve(target)
        if len(search) == 0: 
            continue
        
        # Sourcing only the first available sector to ensure data stability
        lc = search[0].download(quality_bitmask='default')
        if lc is None: 
            continue
            
        # 401-point window flattening as defined in White Paper Section 2
        lc = lc.normalize().flatten(window_length=401).remove_outliers()
        
        flux_std = np.std(lc.flux.value) * 1e6
        max_dip = np.max(np.abs(lc.flux.value - 1.0))

        # --- BLOC C: BLS DETECTION ENGINE ---
        # Frequency factor set to 5.0 for optimized RAM usage on consumer hardware
        bls = lc.to_periodogram(method='bls', period=np.linspace(0.5, 20, 5000), frequency_factor=5.0)
        best_p = bls.period_at_max_power.value
        t0 = bls.transit_time_at_max_power.value
        sde = bls.max_power.value
        duration = bls.duration_at_max_power.value * 24 
        depth = bls.depth_at_max_power.value

        # --- BLOC D: ASTROPHYSICAL CHARACTERIZATION ---
        # 1. Mass & Orbital distance estimation
        m_star = 0.6 * r_star 
        dist_ua = (best_p / 365.25)**(2/3) * (m_star)**(1/3)
        
        # 2. Planetary Radius (Earth Radii)
        r_planet = np.sqrt(depth) * (r_star * 109.1)
        
        # 3. Equilibrium Temperatures (Stefan-Boltzmann derived)
        t_eq = t_eff * np.sqrt(r_star / (2 * (dist_ua * 215))) * (0.7**0.25)
        t_eq_c = t_eq - 273.15
        
        # 4. Habitable Zone Assessment (Insolation flux)
        insolation = (t_eff / 5778)**4 * (r_star / dist_ua)**2
        is_habitable = 0.5 <= insolation <= 2.0
        hz_status = "IN HABITABLE ZONE" if is_habitable else "OUTSIDE HZ"

        # --- BLOC E: AUTOMATED TRIAGE LOGIC ---
        # Criteria for filtering potential Single Transits vs Periodic Signals
        is_single_transit = sde < 10 and (depth > 0.0015) and (max_dip > 3 * np.std(lc.flux.value))

        if is_single_transit:
            confidence = "SINGLE TRANSIT"
            conf_color = "blue"
            icon = "🔵"
        elif sde >= 15:
            confidence = "HIGH"
            conf_color = "green"
            icon = "🟢"
        else:
            confidence = "LOW"
            conf_color = "red"
            icon = "🔴"

        # --- BLOC F: EXOPLANET IDENTITY CARD GENERATION ---
        plt.style.use('default')
        fig = plt.figure(figsize=(10, 12), facecolor='white')
        plt.figtext(0.5, 0.94, "EXOPLANET IDENTITY CARD", fontsize=22, fontweight='bold', ha='center')
        plt.figtext(0.5, 0.91, f"System: {target}", fontsize=16, color='gray', ha='center')
        
        if is_single_transit:
            plt.figtext(0.5, 0.88, "⚠️ POTENTIAL SINGLE TRANSIT DETECTED", fontsize=12, color='blue', fontweight='bold', ha='center', bbox=dict(facecolor='none', edgecolor='blue', pad=5.0))

        hemisphere = "North" if dec >= 0 else "South"
        text_info = (
            f"STAR:\nCoordinates: RA {ra:.4f} / DEC {dec:.4f}\n"
            f"Visibility: {hemisphere} | Mag: {mag:.2f} | Type: {type_spectral}\n"
            f"Stellar Radius: {r_star:.2f} R_Sun | Temperature: {t_eff:.0f} K\n\n"
            f"PLANET:\nRadius: {r_planet:.2f} Earth Radii | HZ Status: {hz_status}\n"
            f"Orbital Distance: {dist_ua:.4f} AU | Equilibrium Temp: {t_eq:.1f} K\n"
            f"Insolation: {insolation:.2f}x Earth flux\n\n"
            f"TRANSIT:\nPeriod: {best_p:.4f} days | Duration: {duration:.2f} hrs\n"
            f"Depth: {depth*100:.3f} % | Detection Score (SDE): {sde:.1f}\n"
        )
        plt.figtext(0.12, 0.86, text_info, fontsize=10, va='top', linespacing=1.6)
        plt.figtext(0.12, 0.58, f"Confidence Level: {confidence}", fontsize=10, fontweight='bold', color=conf_color, va='top')

        ax = fig.add_axes([0.15, 0.12, 0.70, 0.35]) 
        lc_folded = lc.fold(period=best_p, epoch_time=t0)
        lc_folded.scatter(ax=ax, color='#86241a', s=2, alpha=0.3, label="Raw Data")
        bin_size = best_p/100 if not is_single_transit else 0.01
        lc_folded.bin(time_bin_size=bin_size).plot(ax=ax, color='black', lw=2, label="Binned Transit")
        ax.set_title(f"Photometric Signature: {target}", fontsize=14, fontweight='bold', pad=15)
        ax.legend(loc='upper right'); ax.grid(True, linestyle=':', alpha=0.5)

        plt.savefig(os.path.join(export_dir, f"REPORT_{tic}_ARGOS_v2.1.webp"), dpi=120, bbox_inches='tight')
        plt.close()

        # --- BLOC G: DATABASE & LOGGING ---
        # 1. Operational Log
        with open(log_path, "a") as f:
            f.write(f"[{datetime.now().strftime('%H:%M')}] {icon} {target} | SDE: {sde:.1f} | Conf: {confidence}\n")
        
        # 2. Research Database (CSV)
        if not os.path.exists(csv_path):
            with open(csv_path, "w") as f:
                header = "Date;Target;RA;DEC;Hemi;Type;Radius_Re;Dist_AU;T_eq_K;T_eq_C;Insolation;HZ_Status;Dur_Hrs;Noise_PPM;SDE;Mag;Confidence\n"
                f.write(header)
        
        row = (f"{datetime.now().strftime('%Y-%m-%d')};{target};{ra:.4f};{dec:.4f};{hemisphere};{type_spectral};"
               f"{r_planet:.2f};{dist_ua:.4f};{t_eq:.1f};{t_eq_c:.1f};{insolation:.2f};{hz_status};"
               f"{duration:.2f};{flux_std:.0f};{sde:.1f};{mag:.2f};{confidence}\n")
        
        with open(csv_path, "a") as f:
            f.write(row)

        # 3. Processing History
        with open(history_path, "a") as f:
            f.write(f"{tic}\n")
            
        print(f"✅ {target} processed ({confidence})")

    except Exception as e:
        print(f"❌ {target} Critical Error: {e}")

print(f"\n🎯 Benchmark complete. Results available in: {export_dir}")