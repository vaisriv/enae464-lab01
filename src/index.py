import os
import sys
import csv
import numpy as np
import matplotlib.pyplot as plt

def main():
    # ─── Physical constants & ambient conditions ─────────────────────────────────
    P_AMB     = 100900.0 # Ambient pressure [Pa]
    T_AMB     = 298.15   # Ambient temperature [K]
    RHO_AIR   = 1.18     # Air density [kg/m^3]
    RHO_WATER = 998.0    # Water density [kg/m^3]
    G         = 9.81     # Gravitational acceleration [m/s^2]
    RADIUS    = 0.635e-2 # Cylinder Radius [m]

    # ─── File paths ──────────────────────────────────────────────────────────────
    INPUT_CSV  = "./data/pressure_vs_theta.csv"
    OUTPUT_CSV = "./outputs/text/pressure_vs_theta.csv"
    OUTPUT_FILE = "./outputs/text/summary.txt"

    # ─── Read CSV ────────────────────────────────────────────────────────────────
    theta_deg_list = []
    P_inf_raw_list = []
    P_0_raw_list   = []
    P_raw_list     = []

    with open(INPUT_CSV, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            theta_deg_list.append(float(row["theta_deg"]))
            P_inf_raw_list.append(float(row["P_inf"]))
            P_0_raw_list.append(float(row["P_0"]))
            P_raw_list.append(float(row["P"]))

    N = len(theta_deg_list)
    print(f"Read {N} data points from '{INPUT_CSV}'.")

    # ─── Convert water column heights to differential pressures ──────────────────
    # The manometer readings are heights (in cm of water):
    #     ΔP = ρ_water · g · Δh
    #     P_actual = P_ref - ρ_water · g · h_reading
    #     P_surface - P_freestream  =  ρ_water · g · (h_surface - h_inf)
    #     q_inf = P_0 - P_inf = ρ_water · g · (h_0 - h_inf)
    SCALE = 1e-2  # Measurements taken in cm (convert to meters)

    # ─── Compute differential pressures ─────────────────────────────────────────
    # For each data point we have a corresponding P_inf and P_0 reading,
    # so we compute per-row to account for any drift.
    delta_P_list = [] # P_surface - P_freestream  [Pa]
    q_inf_list   = [] # dynamic pressure  [Pa]
    U_inf_list   = [] # freestream velocity [m/s]
    Cp_list      = [] # pressure coefficient

    for i in range(N):
        h_inf = P_inf_raw_list[i] * SCALE # freestream static tap height [m]
        h_0   = P_0_raw_list[i]   * SCALE # stagnation tap height [m]
        h_p   = P_raw_list[i]     * SCALE # surface tap height [m]

        # (P_surface - P_inf) in Pa
        # Higher reading = higher pressure, so:
        delta_P = RHO_WATER * G * (h_p - h_inf)

        # Dynamic pressure  q = P_total - P_static = rho_w * g * (h_0 - h_inf)
        q_inf = RHO_WATER * G * (h_0 - h_inf)

        if q_inf <= 0:
            print(f"WARNING row {i}: q_inf = {q_inf:.2f} Pa <= 0  (h_inf={h_inf:.4f}, h_0={h_0:.4f})")
            # Use absolute value as fallback but flag it
            q_inf = abs(q_inf) if abs(q_inf) > 1e-6 else 1e-6

        U_inf = (2.0 * q_inf / RHO_AIR) ** 0.5 # Bernoulli: q = 0.5 * rho * U^2

        Cp = delta_P / q_inf # Cp = (P - P_inf) / q_inf  = (P - P_inf) / (0.5 * rho * U^2)

        delta_P_list.append(delta_P)
        q_inf_list.append(q_inf)
        U_inf_list.append(U_inf)
        Cp_list.append(Cp)

    # TODO: fig: pressure coefficient vs angle of attack
    # A plot of \(C_p\) vs \(\theta\) from the experiment, plotted simultaneous (on the same graph) with the pressure coefficient that is given by inviscid theory.


    # ─── Summary statistics ─────────────────────────────────────────────────────
    q_avg   = sum(q_inf_list) / N
    U_avg   = sum(U_inf_list) / N

    with open(OUTPUT_FILE, "w", newline="") as f:
        f.write(f"Ambient conditions:")
        f.write(f"  P_amb   = {P_AMB} Pa")
        f.write(f"  T_amb   = {T_AMB} K")
        f.write(f"  rho_air = {RHO_AIR} kg/m^3")
        f.write(f"\nDerived quantities (averages over all rows):")
        f.write(f"  q_inf   = {q_avg:.2f} Pa")
        f.write(f"  U_inf   = {U_avg:.2f} m/s")

    print(f"\nSummary written to '{OUTPUT_FILE}'.")

    # ─── Write output table ─────────────────────────────────────────────────────
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["theta_deg", "P-P_inf_Pa", "Cp", "P_Pa"])

        for i in range(N):
            # P_actual = P_amb + (P_surface - P_inf) i.e. ambient + differential
            P_actual = P_AMB + delta_P_list[i]
            writer.writerow([
                f"{theta_deg_list[i]:.1f}",
                f"{delta_P_list[i]:.3f}",
                f"{Cp_list[i]:.4f}",
                f"{P_actual:.3f}",
            ])

    print(f"\nResults written to '{OUTPUT_CSV}'.")


if __name__ == "__main__":
    main()
