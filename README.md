# 🔋 Tier-1 Battery Digital Twin: Multi-Node Electro-Thermal Simulator

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_svg)](https://okechukwuikenna-battery-physics-sim.streamlit.app)

## 🚀 Overview
This repository hosts a high-fidelity **Electro-Thermal Digital Twin** for Lithium-ion battery cells. Developed as a physics-informed predictive tool, it allows engineers to simulate complex internal cell behavior—such as radial thermal gradients and local anode potential drops—under aggressive fast-charging protocols.



## 🧠 Key Scientific Features
* **Multi-Node Core-Shell Thermal Architecture:** Resolves internal radial temperature gradients by coupling core heat generation with surface-to-ambient convection.
* **Closed-Loop CC-CV Controller:** Simulates real-world Battery Management System (BMS) logic, including dynamic current tapering and voltage clamping.
* **Anode Potential Tracking ($\Phi_{anode}$):** A microscopic look at the anode/electrolyte interface to predict the onset of **Lithium Plating**—the primary driver of fast-charging failure.
* **Dual-Pathway Degradation:** Tracks competing aging mechanisms—**Chemical SEI Growth** (High-Temp) vs. **Metallic Plating** (High-Rate/Low-Temp).
* **Arrhenius Kinetics:** Dynamically adjusts internal resistance based on thermally activated ion transport.



## 🛠️ Tech Stack
* **Core:** Python 3.x
* **Math/Physics:** `NumPy` for kinetics, `SciPy` for non-linear ODE solving.
* **Visualization:** `Matplotlib` for high-resolution electrochemical diagnostic plots.
* **UI/Deployment:** `Streamlit` for a reactive, web-based engineering workspace.

## 📊 Methodology
The simulator couples a 1-RC Equivalent Circuit Model (ECM) with a 2-node lumped capacitance thermal model. By resolving the spatial thermal gradients and local anode electrochemistry, it identifies "hidden" safety risks like core-localized hotspots and anode potential drops that surface sensors cannot detect.

## 🚀 Getting Started
1. **Clone the repo:**
   ```bash
   git clone [https://github.com/okechukwuikenna/battery-physics-sim.git](https://github.com/okechukwuikenna/battery-physics-sim.git)
