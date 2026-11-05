"""Phase A0: 1D freely propagating NH3/H2/N2/air flames and the flame properties
needed for Phases A/B (proposal.md, sections 4.3 and 4.4).

Definitions (proposal.md section 4.3):
  fuel(gamma)  : X_NH3 = (1-g)/(1+g), X_H2 = 1.5 g/(1+g), X_N2 = 0.5 g/(1+g)
  delta_T      : (T_b - T_u) / max|dT/dx|
  l_D          : lambda_u / (rho_u cp_u s_L)
  tau_F        : delta_T / s_L
  Ze           : E/R (T_b - T_u)/T_b^2, E/R = -2 dln(rho_u s_L)/d(1/T_b)
                 (a) Lehmann 2024: inert N2 mass fraction perturbed by +/- dY_N2
                 (b) check     : T_u perturbed by +/- dT_u
  Le_eff       : Joulin-Mitani / Bechtold-Matalon, lean: D = fuel, E = O2,
                 A = 1 + Ze (1/phi - 1); Le_fuel with three weightings
                 (molar = Chabot, diffusivity = Dinkelacker/Lehmann,
                 heat release = Law); species Le evaluated in the unburnt
                 and in the burnt gas (Schneider uses the burnt gas).
  Pe_cd        : Rieth 2022 / Schneider II eq. (7): max|convection| / max|diffusion|
                 of a fuel species; diffusion rate from the steady 1D balance
                 d(j_k)/dx = omega_k W_k - mdot dY_k/dx.
"""
import json
import time

import cantera as ct
import numpy as np

AIR = "O2:0.21, N2:0.79"


def fuel_composition(gamma):
    """Mole fractions of the partially cracked fuel (stoichiometric cracking)."""
    return {"NH3": (1.0 - gamma) / (1.0 + gamma),
            "H2": 1.5 * gamma / (1.0 + gamma),
            "N2": 0.5 * gamma / (1.0 + gamma)}


def x_h2_fuel(gamma):
    """H2 fraction of the N2-free fuel (Lehmann's X_H2,F)."""
    return 1.5 * gamma / (1.0 + 0.5 * gamma)


def set_mixture(gas, gamma, phi, T, P_atm, dY_N2=0.0):
    """Set the unburnt state. dY_N2 shifts the N2 mass fraction at fixed phi
    (the other species are rescaled), as used for Ze."""
    fuel = {k: v for k, v in fuel_composition(gamma).items() if v > 0}
    gas.set_equivalence_ratio(phi, fuel, AIR)
    gas.TP = T, P_atm * ct.one_atm
    if dY_N2:
        Y = gas.Y.copy()
        iN2 = gas.species_index("N2")
        rest = 1.0 - Y[iN2]
        Y_new_N2 = Y[iN2] + dY_N2
        Y *= (1.0 - Y_new_N2) / rest
        Y[iN2] = Y_new_N2
        gas.TPY = T, P_atm * ct.one_atm, Y
    return gas


def initial_width(P_atm, T_u):
    """Domain width [m]: thick lean NH3 flames at low P and low T_u."""
    return 0.06 * (1.0 / P_atm) ** 0.6 * (300.0 / T_u) ** 1.0 + 0.004


def solve_free_flame(gas, width, soret=True, flame=None, loglevel=0):
    """Solve a freely propagating flame. If `flame` is given, it is re-solved
    from its current solution with the new inlet state (continuation)."""
    if flame is None:
        flame = ct.FreeFlame(gas, width=width)
        flame.transport_model = "mixture-averaged"
        flame.soret_enabled = soret
        flame.max_time_step_count = 20000
        # coarse solution first, then the production grid
        flame.set_refine_criteria(ratio=3, slope=0.15, curve=0.3, prune=0.05)
        flame.solve(loglevel=loglevel, auto=True)
        flame.set_refine_criteria(ratio=3, slope=0.04, curve=0.08, prune=0.01)
        flame.solve(loglevel=loglevel, auto=False)
    else:
        flame.inlet.T = gas.T
        flame.inlet.Y = gas.Y
        flame.P = gas.P
        flame.solve(loglevel=loglevel, auto=False)
    return flame


def lewis_numbers(gas, species):
    """Le_k = lambda / (rho cp D_k,mix) at the current gas state."""
    gas.transport_model = "mixture-averaged"
    D = gas.mix_diff_coeffs
    alpha = gas.thermal_conductivity / (gas.density * gas.cp_mass)
    return {k: alpha / D[gas.species_index(k)] for k in species}


def le_fuel(Le, gamma, gas):
    """Fuel Lewis number with three weightings (two-fuel blend NH3/H2)."""
    xh = x_h2_fuel(gamma)
    out = {}
    if gamma >= 1.0:
        return {"molar": Le["H2"], "diff": Le["H2"], "hrr": Le["H2"]}
    if gamma <= 0.0:
        return {"molar": Le["NH3"], "diff": Le["NH3"], "hrr": Le["NH3"]}
    out["molar"] = xh * Le["H2"] + (1 - xh) * Le["NH3"]
    out["diff"] = 1.0 / (xh / Le["H2"] + (1 - xh) / Le["NH3"])
    q = heat_release_per_mole(gas)
    wh, wn = xh * q["H2"], (1 - xh) * q["NH3"]
    out["hrr"] = (wh * Le["H2"] + wn * Le["NH3"]) / (wh + wn)
    return out


_LHV = {}


def heat_release_per_mole(gas):
    """Lower heating values [J/kmol] of H2 and NH3 at 298 K
    (H2 + 1/2 O2 -> H2O ; NH3 + 3/4 O2 -> 3/2 H2O + 1/2 N2)."""
    if not _LHV:
        h = {}
        for k in ("H2", "NH3", "O2", "H2O", "N2"):
            gas.TPX = 298.15, ct.one_atm, f"{k}:1"
            h[k] = gas.enthalpy_mole
        _LHV["H2"] = h["H2"] + 0.5 * h["O2"] - h["H2O"]
        _LHV["NH3"] = h["NH3"] + 0.75 * h["O2"] - 1.5 * h["H2O"] - 0.5 * h["N2"]
    return _LHV


def le_eff(Le_fuel, Le_O2, Ze, phi):
    """Joulin-Mitani effective Lewis number (Lehmann eq. 9-10)."""
    if phi <= 1.0:
        Le_D, Le_E, phi_star = Le_fuel, Le_O2, 1.0 / phi
    else:
        Le_D, Le_E, phi_star = Le_O2, Le_fuel, phi
    A = 1.0 + Ze * (phi_star - 1.0)
    return 1.0 + ((Le_E - 1.0) + (Le_D - 1.0) * A) / (1.0 + A)


def adiabatic_T(mech_gas, gas):
    """Equilibrium (HP) adiabatic flame temperature of the state in `gas`."""
    mech_gas.TPX = gas.T, gas.P, gas.X
    mech_gas.equilibrate("HP")
    return mech_gas.T, mech_gas.density


def flame_props(f, Tb, rho_b):
    """Basic properties of a converged free flame; Tb, rho_b at equilibrium
    (NH3 flames relax slowly, the end of the domain is not at T_ad)."""
    x, T = f.grid, f.T
    rho_u = f.density[0]
    sL = f.velocity[0]
    dTdx = np.gradient(T, x)
    dT = (Tb - T[0]) / np.max(np.abs(dTdx))
    return {"sL": sL, "rho_u": rho_u, "rho_b": rho_b, "Tu": T[0], "Tb": Tb,
            "T_end": T[-1],
            "delta_T": dT, "x_Tmaxgrad": x[np.argmax(np.abs(dTdx))],
            "hrr_max": float(np.max(f.heat_release_rate)), "npts": len(x)}


def pe_cd(f, gas, species):
    """Max |convection rate| / max |diffusion rate| of a species (Rieth/Schneider)."""
    k = gas.species_index(species)
    x = f.grid
    mdot = f.density[0] * f.velocity[0]
    Y = f.Y[k]
    conv = mdot * np.gradient(Y, x)
    src = f.net_production_rates[k] * gas.molecular_weights[k]
    diff = src - conv          # = d(j_k)/dx in a steady 1D flame
    return float(np.max(np.abs(conv)) / np.max(np.abs(diff)))


def outlet_mole_fractions(f, gas, species):
    return {k: float(f.X[gas.species_index(k), -1]) for k in species}


def run_case(mech, gamma, phi, P_atm, T_u, soret=True, dY_N2=0.003, dT_u=10.0,
             profile_path=None, loglevel=0, phi_start=None, dphi=0.02):
    """All A0 quantities for one (gamma, phi, P, T_u). Returns a dict.
    phi_start: if given, the flame is first solved at phi_start and then
    continued to phi in steps of dphi (cases whose cold start fails)."""
    t0 = time.time()
    gas = ct.Solution(mech)
    res = {"mech": mech, "gamma": gamma, "phi": phi, "P_atm": P_atm, "T_u": T_u,
           "soret": soret, "X_H2_F": x_h2_fuel(gamma), "status": "ok"}

    set_mixture(gas, gamma, phi, T_u, P_atm)
    gas_u = {"X": dict(zip(gas.species_names, gas.X))}
    res["X_u"] = {k: v for k, v in gas_u["X"].items() if v > 0}
    lam_u, cp_u, rho_u, mu_u = (gas.thermal_conductivity, gas.cp_mass,
                                gas.density, gas.viscosity)
    res["Pr_u"] = mu_u * cp_u / lam_u
    Le_u = lewis_numbers(gas, ["H2", "NH3", "O2"])
    g_eq = ct.Solution(mech)
    Tb, rho_b = adiabatic_T(g_eq, gas)

    if phi_start is None:
        f = solve_free_flame(gas, initial_width(P_atm, T_u), soret=soret, loglevel=loglevel)
    else:
        g_c = ct.Solution(mech)
        set_mixture(g_c, gamma, phi_start, T_u, P_atm)
        f = solve_free_flame(g_c, initial_width(P_atm, T_u), soret=soret, loglevel=loglevel)
        step = -dphi if phi < phi_start else dphi
        for ph in np.append(np.arange(phi_start + step, phi, step), phi):
            set_mixture(g_c, gamma, ph, T_u, P_atm)
            solve_free_flame(g_c, None, soret=soret, flame=f, loglevel=loglevel)
        res["phi_start"] = phi_start
    base = flame_props(f, Tb, rho_b)
    res.update(base)
    res["sigma"] = base["rho_u"] / base["rho_b"]
    res["l_D"] = lam_u / (rho_u * cp_u * base["sL"])
    res["tau_F"] = base["delta_T"] / base["sL"]
    res["delta_ratio"] = base["delta_T"] / res["l_D"]   # Lehmann's delta = l_F / l_D
    if profile_path:
        f.save(profile_path, basis="mole", overwrite=True)

    # burnt-gas Lewis numbers
    gas.TPY = f.T[-1], f.P, f.Y[:, -1]
    Le_b = lewis_numbers(gas, ["H2", "NH3", "O2"])
    res["Le_u"], res["Le_b"] = Le_u, Le_b

    # Pe_cd for both fuels (H2 = Schneider convention)
    res["Pe_cd_H2"] = pe_cd(f, gas, "H2")
    res["Pe_cd_NH3"] = pe_cd(f, gas, "NH3") if gamma < 1.0 else float("nan")

    # outlet (end of domain) composition: equilibrium-like, for reference only
    res["X_out"] = outlet_mole_fractions(f, gas, ["NO", "N2O", "NO2", "NH3", "H2", "O2", "H2O"])

    # Ze: continuation from the current solution (small perturbations)
    def burning_flux(dY=0.0, dT=0.0):
        set_mixture(gas, gamma, phi, T_u + dT, P_atm, dY_N2=dY)
        Tb_p, _ = adiabatic_T(g_eq, gas)
        solve_free_flame(gas, None, soret=soret, flame=f, loglevel=loglevel)
        return f.density[0] * f.velocity[0], Tb_p

    ze = {}
    for tag, kw in (("N2", [{"dY": +dY_N2}, {"dY": -dY_N2}]),
                    ("Tu", [{"dT": +dT_u}, {"dT": -dT_u}])):
        try:
            (m1, Tb1), (m2, Tb2) = burning_flux(**kw[0]), burning_flux(**kw[1])
            ER = -2.0 * (np.log(m1) - np.log(m2)) / (1.0 / Tb1 - 1.0 / Tb2)
            ze[tag] = {"E_R": ER, "Ze": ER * (base["Tb"] - T_u) / base["Tb"] ** 2}
        except Exception as e:  # noqa: BLE001
            ze[tag] = {"E_R": float("nan"), "Ze": float("nan"), "error": str(e)[:200]}
    res["Ze_method"] = ze
    Ze = ze["N2"]["Ze"]
    res["Ze"] = Ze

    # Le_eff: 3 fuel weightings x (unburnt, burnt)
    res["Le_fuel"], res["Le_eff"] = {}, {}
    for side, Le in (("u", Le_u), ("b", Le_b)):
        lf = le_fuel(Le, gamma, gas)
        res["Le_fuel"][side] = lf
        res["Le_eff"][side] = {w: le_eff(v, Le["O2"], Ze, phi) for w, v in lf.items()}
    res["Ze_Le_eff_m1"] = Ze * (res["Le_eff"]["u"]["diff"] - 1.0)
    res["Ze_over_Pe_cd"] = Ze / res["Pe_cd_H2"]
    res["I0_schneider"] = 0.978 + 0.0115 * res["Ze_over_Pe_cd"]
    res["wall_s"] = time.time() - t0
    return res


def to_json(res):
    def conv(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        raise TypeError(type(o))
    return json.dumps(res, default=conv)
