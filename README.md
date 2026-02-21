# ARGOS: Automated Reporting & General Observation System
**An Autonomous Framework for Multi-Domain Astronomical Detection**

## 🌌 Overview
ARGOS is an independent framework designed for high-precision "Triage" and automated reporting of astronomical data. It bridges the gap between raw public data (TESS, Gaia) and actionable scientific candidates.

### 🔬 Current Module: EXOLAB v2.1
The **EXOLAB** module is the first specialized branch of the ARGOS framework, optimized for detecting exoplanetary transits in TESS data.

**Key Features:**
* **Fault-Tolerant Architecture:** Optimized for stability over long processing runs.
* **Single-Sector Processing:** Ensures signal integrity by avoiding multi-sector artifacts.
* **Automated Triage Logic:** Integrated filtering to distinguish between High-Confidence candidates, Single Transits, and false positives.

---

## 🛠 Modular Architecture (The ARGOS Roadmap)
ARGOS is built as a modular "Core" capable of supporting various specialized laboratories:

| Module | Status | Detection Focus |
| :--- | :--- | :--- |
| **EXOLAB** | 🟢 **Active** | Exoplanetary transits & Habitability |
| **FLARELAB** | 🟡 *Planned* | Stellar flare energetics |
| **NOVALAB** | 🟡 *Planned* | Transient events & Supernovae |

---

## 🚀 Usage
```bash
# To run the exoplanet detection engine:
python3 scripts/argos_exolab_v2.1.py
