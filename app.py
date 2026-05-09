import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# --- PAGE CONFIGURATION ---
st.set_page_config(page_title="Tier-1 Battery Electro-Thermal Simulator", layout="wide")
st.title("🔋 Tier-1 Multi-Node Electro-Thermal & Degradation Simulator")
st.markdown("""
This advanced workspace implements a **Core-Shell 2-Node Thermal Model** coupled with a 
**1-RC Equivalent Circuit Model (ECM)**, closed-loop **CC-CV charging control**, and 
**local anode potential tracking** to evaluate true physical lithium-plating thresholds.
""")

# --- SIDEBAR CONTROLS ---
st.sidebar.header("⚡ Charge Controller Setup")
I_cc = st.sidebar.slider("CC Stage Current [A]", 1.0, 40.0, 20.0, step=1.0)
V_max = st.sidebar.slider("Voltage Cut-off (CV Clamping) [V]", 4.0, 4.4, 4.2, step=0.05)
I_cutoff = 0.5  # Amps to terminate CV charge

st.sidebar.header("🌬️ Thermal Management")
h = st.sidebar.slider("Convective Cooling (h) [W/m²/K]", 2.0, 150.0, 15.0, step=5.0)
T_amb_c = st.sidebar.slider("Ambient Air Temp [°C]", -10.0, 50.0, 25.0, step=5.0)
T_amb = T_amb_c + 273.15

with st.sidebar.expander("🔬 Cell Micro-Structure & Chemistry", expanded=False):
    st.markdown("### Thermal Network")
    C_core = st.slider("Core Heat Cap. [J/K]", 100.0, 1000.0, 400.0, step=50.0)
    C_surf = st.slider("Surface Heat Cap. [J/K]", 50.0, 500.0, 150.0, step=50.0)
    R_cond = st.slider("Radial Conduction Res. [K/W]", 0.1, 5.0, 1.2, step=0.1)
    
    st.markdown("### Electrochemical Kinetics")
    R_ref = st.slider("Bulk Resistance (at 298K) [Ω]", 0.005, 0.050, 0.012, step=0.001, format="%.3f")
    Ea_R = st.slider("Kinetics Activation Energy [kJ/mol]", 10.0, 50.0, 24.0, step=2.0) * 1000.0
    R_anode_share = st.slider("Anode Fraction of Total Resistance", 0.1, 0.5, 0.3, step=0.05)
    
    st.markdown("### Degradation Mechanics")
    k_sei = st.number_input("Thermal SEI Growth Rate", value=2.0e-6, format="%.2e")
    k_plating = st.number_input("Metallic Plating Rate", value=8.0e-6, format="%.2e")

# --- CORE PHYSICAL ENGINE ---
# Normalized OCV curve fit for typical high-energy NMC cell
def get_ocv(soc):
    # Safe polynomial + exponential fit for clean derivative/smooth transition
    return 3.2 + 0.5 * soc + 0.4 * (soc ** 2) + 0.1 * (soc ** 5) - 0.15 * np.exp(-25 * soc)

def get_R_int(T_core):
    return R_ref * np.exp((Ea_R / 8.314) * (1.0 / T_core - 1.0 / 298.15))

def run_single_charge_simulation(I_cc_input, h_input):
    # Simulation Parameters
    dt = 1.0  # 1-second steps for high-fidelity convergence
    max_steps = 7200  # 2 hours max run
    
    # Initialize states
    soc = 0.1
    T_core = T_amb
    T_surf = T_amb
    V_RC = 0.0
    
    # Component Values
    A = 0.012  # Surface Area [m^2]
    C_nom = 5.0 * 3600.0  # 5 Ah capacity in Amp-seconds
    R_conv = 1.0 / (h_input * A)
    
    # Transient RC branch values
    R_RC = 0.008
    tau_RC = 45.0  # 45-second diffusion time constant
    
    # History logs
    history = {
        "time": [], "I": [], "V_cell": [], "SoC": [], 
        "T_core": [], "T_surf": [], "Phi_anode": [], "Q_gen": []
    }
    
    for step in range(max_steps):
        R_int = get_R_int(T_core)
        V_ocv = get_ocv(soc)
        
        # Check CC-CV status
        # Terminal voltage equation: V = OCV + I*R + V_RC
        I_trial = I_cc_input
        V_trial = V_ocv + I_trial * R_int + V_RC
        
        if V_trial >= V_max:
            # We are in CV Mode: Solve for I that keeps V_cell = V_max
            I = (V_max - V_ocv - V_RC) / R_int
            I = max(0.0, min(I, I_cc_input))  # Prevent non-physical current
        else:
            # We are in CC Mode
            I = I_cc_input
            
        # Stopping condition (fully charged, current tapered to cutoff)
        if soc >= 0.99 or (V_trial >= V_max and I < I_cutoff):
            break
            
        # Update Electrochemical states
        soc += (I / C_nom) * dt
        soc = min(1.0, soc)
        
        # RC Overpotential State Equation
        V_RC += ((I * R_RC - V_RC) / tau_RC) * dt
        V_cell = V_ocv + I * R_int + V_RC
        
        # Compute dynamic core-shell thermodynamics
        Q_gen = (I**2) * R_int + I * T_core * (-0.1e-3)  # Joule heating + entropic loss
        
        dT_core = ((Q_gen - (T_core - T_surf) / R_cond) / C_core) * dt
        dT_surf = (((T_core - T_surf) / R_cond - (T_surf - T_amb) / R_conv) / C_surf) * dt
        
        T_core += dT_core
        T_surf += dT_surf
        
        # Physics-based Local Anode Potential vs Li/Li+
        U_anode_eq = 0.35 - 0.3 * soc  # Anode drops towards 0.05V as lithium intercalates
        R_anode = R_int * R_anode_share
        Phi_anode = U_anode_eq - I * R_anode
        
        # Append stats
        history["time"].append(step * dt)
        history["I"].append(I)
        history["V_cell"].append(V_cell)
        history["SoC"].append(soc)
        history["T_core"].append(T_core)
        history["T_surf"].append(T_surf)
        history["Phi_anode"].append(Phi_anode)
        history["Q_gen"].append(Q_gen)
        
    return {k: np.array(v) for k, v in history.items()}

# --- RUN ENGINE ---
with st.spinner("Executing real-time physical solvers..."):
    # Run user configuration
    cust = run_single_charge_simulation(I_cc, h)
    
    # Run standard baseline (5A, moderate cooling h=10) for delta comparisons
    base = run_single_charge_simulation(5.0, 10.0)

# --- LONG-TERM CYCLING DEGRADATION LOOP ---
# Integrates chemistry metrics over the charging profiles
def run_lifetime_cycling(sim_data, k_s, k_p):
    Q_health = 1.0
    Q_history = []
    
    # Extract cycle metrics
    I_profile = sim_data["I"]
    T_profile = sim_data["T_core"]
    Phi_profile = sim_data["Phi_anode"]
    
    # 1. Thermal SEI growth rate (Arrhenius dependency on core temperature)
    Ea_sei = 38000.0  # J/mol
    sei_decay = k_s * np.mean(np.exp(-Ea_sei / (8.314 * T_profile)))
    
    # 2. Metallic plating rate (activates strictly when local anode potential < 0V)
    plating_active = np.where(Phi_profile < 0.0, -Phi_profile, 0.0)
    plating_decay = k_p * np.mean(plating_active) if np.any(plating_active > 0) else 0.0
    
    total_decay_rate = sei_decay + plating_decay
    
    for _ in range(300): # simulate 300 cycles
        Q_health *= (1.0 - total_decay_rate)
        Q_history.append(Q_health)
        
    return np.array(Q_history), sei_decay, plating_decay

Q_base_hist, base_sei, base_plate = run_lifetime_cycling(base, k_sei, k_plating)
Q_cust_hist, cust_sei, cust_plate = run_lifetime_cycling(cust, k_sei, k_plating)

# --- KEY PERFORMANCE METRICS ---
col1, col2, col3, col4 = st.columns(4)

total_charge_time_min = cust["time"][-1] / 60.0
base_time_min = base["time"][-1] / 60.0
time_delta = total_charge_time_min - base_time_min
col1.metric("Total Charge Time", f"{total_charge_time_min:.1f} min", f"{time_delta:.1f} min vs Baseline")

max_core_temp = np.max(cust["T_core"]) - 273.15
temp_delta = max_core_temp - (np.max(base["T_core"]) - 273.15)
col2.metric("Max Core Temperature", f"{max_core_temp:.1f} °C", f"{temp_delta:+.1f} °C vs Baseline")

min_anode_potential = np.min(cust["Phi_anode"]) * 1000.0
col3.metric(
    "Minimum Anode Potential", 
    f"{min_anode_potential:.1f} mV", 
    "🚨 PLATING RISK" if min_anode_potential < 0.0 else "✅ SAFE",
    delta_color="inverse" if min_anode_potential < 0.0 else "normal"
)

soh_300 = Q_cust_hist[-1] * 100.0
soh_delta = soh_300 - (Q_base_hist[-1] * 100.0)
col4.metric("State of Health (300 cycles)", f"{soh_300:.1f} %", f"{soh_delta:.1f}% vs Baseline")

st.divider()

# --- DIAGNOSTIC PLOTS ---
tab1, tab2, tab3 = st.tabs(["⚡ Single Cycle Charge Dynamics", "🔬 Local Anode Potential & Degradation", "📈 Long-Term Capacity Degradation"])

with tab1:
    col_t1, col_t2 = st.columns(2)
    with col_t1:
        # Dynamic Current / Voltage curves
        fig, ax1 = plt.subplots(figsize=(6, 3.5))
        ax2 = ax1.twinx()
        
        ax1.plot(cust["time"]/60.0, cust["I"], color="#0984e3", label="Current (I)", linewidth=2.5)
        ax2.plot(cust["time"]/60.0, cust["V_cell"], color="#d63031", label="Terminal Voltage (V)", linewidth=2.5)
        
        ax1.set_xlabel("Time [minutes]")
        ax1.set_ylabel("Current [A]", color="#0984e3")
        ax2.set_ylabel("Voltage [V]", color="#d63031")
        ax1.tick_params(axis='y', labelcolor="#0984e3")
        ax2.tick_params(axis='y', labelcolor="#d63031")
        
        plt.title("CC-CV Charging Profile")
        ax1.grid(True, alpha=0.3)
        st.pyplot(fig)
        
    with col_t2:
        # Core vs Surface thermal profiles
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.plot(cust["time"]/60.0, cust["T_core"] - 273.15, color="#d63031", label="Core Temp", linewidth=2.5)
        ax.plot(cust["time"]/60.0, cust["T_surf"] - 273.15, color="#e17055", label="Surface Temp", linestyle="--", linewidth=2)
        ax.axhline(T_amb_c, color="gray", linestyle=":", label="Ambient Boundary")
        
        ax.set_xlabel("Time [minutes]")
        ax.set_ylabel("Temperature [°C]")
        ax.set_title("Core-Shell Core vs Surface Temp Gradient")
        ax.legend()
        ax.grid(True, alpha=0.3)
        st.pyplot(fig)

with tab2:
    col_d1, col_d2 = st.columns(2)
    with col_d1:
        # Local Anode Potential
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.plot(cust["time"]/60.0, cust["Phi_anode"] * 1000.0, color="#2ecc71", label="Selected Charge Rate", linewidth=2.5)
        ax.plot(base["time"]/60.0, base["Phi_anode"] * 1000.0, color="gray", linestyle="--", label="Baseline (5A)")
        ax.axhline(0.0, color="#d63031", linestyle="-.", label="Lithium Plating Limit (0 mV)")
        
        ax.set_xlabel("Time [minutes]")
        ax.set_ylabel("Local Anode Potential [mV vs Li/Li+]")
        ax.set_title("Lithium Plating Threshold Warning Tool")
        ax.legend()
        ax.grid(True, alpha=0.3)
        st.pyplot(fig)
        
    with col_d2:
        # Bar Chart of Competing Degradation Pathways
        fig, ax = plt.subplots(figsize=(6, 3.5))
        mechanisms = ["Thermal SEI (Cust)", "Metallic Plating (Cust)", "Thermal SEI (Base)", "Metallic Plating (Base)"]
        values = [cust_sei * 1e6, cust_plate * 1e6, base_sei * 1e6, base_plate * 1e6]
        colors = ["#fdcb6e", "#e17055", "#ffeaa7", "#ff7675"]
        
        ax.bar(mechanisms, values, color=colors)
        ax.set_ylabel("Degradation Speed Indicator [ppm/cycle]")
        ax.set_title("Competing Lifetime Wear Mechanisms")
        plt.xticks(rotation=20, ha="right")
        ax.grid(True, axis='y', alpha=0.3)
        st.pyplot(fig)

with tab3:
    fig, ax = plt.subplots(figsize=(10, 3.5))
    ax.plot(range(1, 301), Q_base_hist * 100.0, label="Baseline Setup (5A)", color="gray", linestyle="--")
    ax.plot(range(1, 301), Q_cust_hist * 100.0, label=f"Custom Setup ({I_cc}A)", color="#0984e3", linewidth=3)
    ax.set_ylabel("State of Health (Capacity) [%]")
    ax.set_xlabel("Cycle")
    ax.set_title("300-Cycle Predictive Aging Curve")
    ax.legend()
    ax.grid(True, alpha=0.3)
    st.pyplot(fig)