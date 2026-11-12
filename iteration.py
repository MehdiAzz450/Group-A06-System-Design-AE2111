"""
AE2111-I  WP3.4  -  Class I <-> Class II iteration to a converged design
Group A06 - reasonably sized passenger aircraft

Loop (one iteration):
    1. design point: W/S from the landing constraint and T/W = 0.316 from the
       updated wing/thrust loading (matching) diagram, both fixed in the loop
    2. resize for the current MTOM: wing area and span, tail areas (volume
       coefficients), fuel tank capacity (REQ-RAN-02), required thrust
    3. Class II drag for that geometry   -> class2_drag.class2_drag()  (C_D0, e, L/D)
    4. Class II OEM for that geometry    -> class2_weight.class2()
    5. Class I MTOM and fuel with the Class II OEM and L/D -> class1_weight.class1()
    repeat until the MTOM changes less than `tol` (relative).

Assumptions (stated in the report, Table "iteration assumptions"):
    - The landing constraint sets W/S; it does not depend on MTOM, C_D0 or e.
    - T/W = 0.316 is read from the matching diagram redrawn with the Class II polar;
      it only sets the required thrust (the engine and its mass are fixed).
    - Planform shape (AR, taper, sweep, t/c), tail arm and fuselage are fixed;
      the same engine (PW2037) is kept, so engine mass is constant.

Usage
-----
    python iteration.py
    from iteration import iterate;  res = iterate()
"""

import math

import class1_weight as c1
import class2_drag as cd
import class2_weight as c2

G = 9.81

LOOP = dict(
    WS=8303.0,               # N/m2 design wing loading (landing limit 8344)   [WP2 5.1.2, compliance matrix]
    TW=0.316,                # -    design T/W, updated matching diagram       [WP3.4 matching diagram]
    eta_WP1=0.357,           # -    WP1 engine efficiency, ONLY for the WP1 start point [WP1 Excel D13]
    use_class2_drag=True,    # -    False -> keep the Class I drag polar
    V_h=0.70,                # -    horizontal tail volume coefficient         [WP3 3.4.2]
    V_v=0.105,               # -    vertical tail volume coefficient           [WP3 3.4.3]
    V_wing_fuel_ref=21.15,   # m3   wing fuel volume at S = 119.33 m2          [WP2 5.1.2]
    S_ref=119.33,            # m2   wing area belonging to V_wing_fuel_ref     [WP2 Table 5.2]
    T_available=2 * 165.5e3, # N    2x PW2037 take-off thrust                  [WP3 2.4]
    payload_B=8_531.0,       # kg   payload at MTOM with full fuel             [REQ-RAN-02]
    R_B=11_716e3,            # m    range at MTOM with full fuel               [REQ-RAN-02]
    R_ferry=12_697e3,        # m    ferry range                                [REQ-RAN-03]
    LD_req=18.0,             # -    cruise L/D                                 [REQ-SYS-WING-06]
    D_wing_req=103_159.0,    # N    wing cruise drag                           [REQ-SYS-WING-05]
)


def mac(S, b, taper):
    c_r = 2 * S / (b * (1 + taper))
    return 2 / 3 * c_r * (1 + taper + taper ** 2) / (1 + taper)


def resize(mtom, V_fuel, p2=c2.INPUTS, loop=LOOP):
    """Class II inputs that change with MTOM (SI)."""
    S = mtom * G / loop["WS"]
    b = math.sqrt(p2["AR"] * S)
    c_bar = mac(S, b, p2["taper"])
    return dict(
        MTOM=mtom,
        S_w=S,
        b_w=b,
        S_ht=loop["V_h"] * S * c_bar / p2["L_t"],
        S_vt=loop["V_v"] * S * b / p2["L_t"],
        V_fuel=V_fuel,                    # total fuel tank capacity [m3] (Raymer 15.34 V_t, Torenbeek V_ft)
    )


def evaluate(mtom, V_cap, p1, p2, loop):
    """One pass: geometry, Class II drag and weight, Class I with the Class II OEM and L/D."""
    geo = resize(mtom, V_cap, p2, loop)
    p2i = {**p2, **geo}
    dr = cd.class2_drag(p2i)
    p1i = {**p1, "LD_cruise": dr["LD_max"]} if loop["use_class2_drag"] else p1
    r2 = c2.class2(p2i)
    r1 = c1.class1(p1i, oem=r2["OEM"])
    pr = payload_range(r1, p1i, loop)
    return geo, p1i, p2i, dr, r2, r1, pr


def payload_range(r1, p1i, loop=LOOP):
    """Fuel capacity from REQ-RAN-02 and the corner points of the payload-range diagram."""
    mtom, oem = r1["MTOM"], r1["OEM"]
    ff_B, _ = c1.fuel_fraction({**p1i, "R_nom": loop["R_B"]})
    fuel_cap = max(ff_B * mtom, r1["fuel"])                     # tank capacity needed [kg]
    pl_B = mtom - oem - fuel_cap
    tom_C = oem + fuel_cap
    return dict(
        fuel_cap=fuel_cap, V_cap=fuel_cap / p1i["rho_fuel"],
        A=dict(payload=r1["payload"], R=c1.nominal_range(p1i, r1["fuel"] / mtom)),
        B=dict(payload=pl_B, R=c1.nominal_range(p1i, fuel_cap / mtom)),
        C=dict(payload=0.0, R=c1.nominal_range(p1i, fuel_cap / tom_C)),
    )


def iterate(p1=c1.INPUTS, p2=c2.INPUTS, loop=LOOP, tol=1e-4, max_iter=50):
    """Converged Class I / II result and the iteration history."""
    p1_wp1 = {**p1, "eta_eng": loop["eta_WP1"]}
    start = c1.class1(p1_wp1)                                    # WP1 statistical estimate
    mtom = start["MTOM"]
    V_cap = payload_range(start, p1_wp1, loop)["V_cap"]
    history = [dict(it=0, OEM=start["OEM"], fuel=start["fuel"], change=None,
                    CD0=start["derived"]["CD0"], LD=start["derived"]["LD"], TW=loop["TW"], V_cap=V_cap,
                    **resize(mtom, V_cap, p2, loop))]
    for i in range(1, max_iter + 1):
        geo, p1i, p2i, dr, r2, r1, pr = evaluate(mtom, V_cap, p1, p2, loop)
        change = (r1["MTOM"] - mtom) / mtom
        change_V = (pr["V_cap"] - V_cap) / V_cap
        mtom, V_cap = r1["MTOM"], pr["V_cap"]
        history.append(dict(it=i, OEM=r2["OEM"], fuel=r1["fuel"], change=change, CD0=dr["CD0"], LD=dr["LD_max"],
                            TW=loop["TW"], V_cap=V_cap, **resize(mtom, V_cap, p2, loop)))
        if abs(change) < tol and abs(change_V) < tol:
            break
    else:
        raise RuntimeError("Class I / II iteration did not converge")
    geo, p1i, p2i, dr, r2, r1, pr = evaluate(mtom, V_cap, p1, p2, loop)
    TW = loop["TW"]
    m_cruise = r1["MTOM"] * (1 - 0.5 * r1["fuel_fraction"])     # mid-cruise mass
    checks = dict(
        TW=TW,
        T_required=TW * r1["MTOM"] * G,
        T_available=loop["T_available"],
        V_fuel=r1["V_fuel"],
        V_wing_fuel=loop["V_wing_fuel_ref"] * (geo["S_w"] / loop["S_ref"]) ** 1.5,
        MAC=mac(geo["S_w"], geo["b_w"], p2["taper"]),
        cruise=cd.cruise_point(dr, p2i, m_cruise),
    )
    return dict(class1=r1, class2=r2, drag=dr, inputs1=p1i, inputs2=p2i, geometry=geo,
                history=history, start=start, checks=checks, payload_range=pr)


def print_iteration(res):
    print("=" * 104)
    print(" CLASS I <-> CLASS II ITERATION")
    print("=" * 104)
    print(f"{'it':>3}{'MTOM [kg]':>12}{'OEM [kg]':>11}{'fuel [kg]':>11}{'S_w [m2]':>10}{'b [m]':>8}"
          f"{'S_h [m2]':>10}{'S_v [m2]':>10}{'C_D0':>9}{'L/D':>7}{'T/W':>8}{'change':>10}")
    for h in res["history"]:
        ch = "-" if h["change"] is None else f"{100 * h['change']:.3f}%"
        print(f"{h['it']:>3}{h['MTOM']:>12,.0f}{h['OEM']:>11,.0f}{h['fuel']:>11,.0f}{h['S_w']:>10.2f}"
              f"{h['b_w']:>8.2f}{h['S_ht']:>10.2f}{h['S_vt']:>10.2f}{h['CD0']:>9.5f}{h['LD']:>7.2f}"
              f"{h['TW']:>8.4f}{ch:>10}")
    s, r1, ck, pr = res["start"], res["class1"], res["checks"], res["payload_range"]
    print("-" * 104)
    print(f"{'':<28}{'Class I (WP1)':>14}{'converged':>14}{'change':>10}")
    for name, a, b in [("MTOM [kg]", s["MTOM"], r1["MTOM"]), ("OEM [kg]", s["OEM"], r1["OEM"]),
                       ("Fuel [kg]", s["fuel"], r1["fuel"]),
                       ("OEM / MTOM [-]", s["OEM_fraction"], r1["OEM_fraction"]),
                       ("Fuel fraction [-]", s["fuel_fraction"], r1["fuel_fraction"]),
                       ("Fuel volume [m3]", s["V_fuel"], r1["V_fuel"])]:
        print(f"{name:<28}{a:>14,.3f}{b:>14,.3f}{100 * (b - a) / a:>9.1f}%")
    cr = ck["cruise"]
    print(f"\nThrust:  T/W {ck['TW']:.4f}, required {ck['T_required'] / 1e3:.1f} kN, available "
          f"{ck['T_available'] / 1e3:.1f} kN (margin {100 * (ck['T_available'] / ck['T_required'] - 1):.1f} %)")
    print(f"Fuel:    harmonic {ck['V_fuel']:.2f} m3, capacity needed for REQ-RAN-02 {pr['V_cap']:.2f} m3, "
          f"wing tanks {ck['V_wing_fuel']:.2f} m3 (scaled with S^1.5)")
    print(f"Cruise:  C_L {cr['CL']:.3f}, L/D {cr['LD']:.2f} (incl. wave drag), wing drag {cr['D_wing'] / 1e3:.1f} kN")
    for k, name in [("A", "harmonic"), ("B", "MTOM, full fuel"), ("C", "ferry")]:
        print(f"Payload-range {k} ({name}): payload {pr[k]['payload']:,.0f} kg, range {pr[k]['R'] / 1e3:,.0f} km")
    print(f"Wing MAC {ck['MAC']:.3f} m")


def check_calibration(p2=c2.INPUTS, tol=0.005):
    """Verification: the factors in class2_weight.INPUTS must equal the ones that
    calibrate_raymer.py derives with the current equations and shared inputs."""
    import calibrate_raymer as cr
    mean, _ = cr.calibrate(verbose=False)
    bad = {g: (p2["fudge"].get(g, 1.0), round(f, 3)) for g, f in mean.items()
           if abs(p2["fudge"].get(g, 1.0) - f) > tol}
    if bad:
        print(f"WARNING: calibration factors out of date (used, recalculated): {bad}")
    else:
        print("Calibration factors match calibrate_raymer.py.")


if __name__ == "__main__":
    check_calibration()
    res = iterate()
    print_iteration(res)
