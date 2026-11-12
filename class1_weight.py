"""
AE2111-I  WP3.4  -  Class I weight estimation (ADSEE fuel-fraction method)
Group A06 - reasonably sized passenger aircraft

Method: TU Delft ADSEE Class I estimation.
    MTOM = (m_OE + m_pl) / (1 - m_f/MTOM),   m_OE = f_OE * MTOM   (statistics)
    m_f/MTOM = 1 - exp(-R_eq / (eta_eng * (e_f / g) * (L/D)))
    R_eq   = (R_nom + R_lost) * (1 + f_con) + 1.2 * R_div + t_E * V_cr
    R_lost = (1 / 0.7) * (L/D) * (h_cr + V_cr^2 / (2 g))
    eta_eng: fixed value (WP1), or V_cr / (TSFC * e_f) with TSFC = 22 * BPR^-0.19 [g/(kN s)]
    L/D = (L/D)_max = 0.5 * sqrt(pi * AR * e / C_D0),   C_D0 = C_f,eq * S_wet/S

The defaults reproduce the WP1 Excel ("Aircraft Design and Calculations - Y2.xlsx",
sheet "class I weight estimation", column "Harmonic (MTOM)").

Usage
-----
    python class1_weight.py                  -> prints the Class I result
    from class1_weight import class1, INPUTS
    res = class1(INPUTS, oem=41_007)          -> fixed OEM from Class II (WP3.4 loop)

Every input carries a comment with its source:
    [WP1]      -> value from the WP1 Excel (cell reference given)
    [WP2]      -> value from the WP2 report
"""

import math

G = 9.81                     # m/s2, as in the WP1 Excel
R_AIR = 287.05
GAMMA = 1.4

# ================================================================== INPUTS (SI)
INPUTS = dict(
    # ---- mission (harmonic: max payload over the harmonic range)
    payload=18_960.0,        # kg   max payload                                [WP1 Ref!E20]
    R_nom=9_545e3,           # m    harmonic range                             [WP1 Ref!H20]
    h_cr=9_448.8,            # m    cruise altitude FL310                      [WP1 D8]
    V_cr=257.0,              # m/s  cruise speed (M 0.85); None -> ISA at h_cr [WP1 D8]
    M_cr=0.85,               # -    only used if V_cr is None                  [WP1]
    f_con=0.05,              # -    contingency fuel fraction                  [WP1 D11]
    R_div=400e3,             # m    diversion range                            [WP1 D9]
    t_E=1800.0,              # s    loiter endurance                           [WP1 D9]

    # ---- statistics
    f_OE=0.50,               # -    OEM / MTOM                                 [WP1 D16]

    # ---- aerodynamics
    LD_cruise=None,          # -    set a value to bypass (L/D)_max below
    AR=10.5,                 # -                                               [WP1 Ref!R20]
    Cf_eq=0.0027,            # -    equivalent skin friction coefficient       [WP1 D5]
    Swet_S=6.5,              # -    wetted area ratio                          [WP1 D4]
    e=0.76,                  # -    Oswald efficiency factor                   [WP1 D6]

    # ---- propulsion
    eta_eng=None,            # -    None -> eta = V_cr / (TSFC e_f)            [reader Eq. 5.21]
                             #      (WP1 used 0.357; replaced for the selected PW2037)
    TSFC=None,               # kg/(N s)  None -> 22 B^-0.19 g/(kN s)          [reader Eq. 5.23]
    BPR=6.0,                 # -    PW2037 bypass ratio                        [WP3.1]
    e_f=44e6,                # J/kg specific energy of kerosene                [WP1 D13]
    rho_fuel=800.0,          # kg/m3                                           [WP1 D25]

    # ---- for comparison only
    MTOM_ref=100_999.0,      # kg   WP1 Excel / WP2 final value                [WP1 D20]
)


# ================================================================== helpers
def isa_speed_of_sound(h):
    """ISA speed of sound [m/s], troposphere."""
    return math.sqrt(GAMMA * R_AIR * (288.15 - 0.0065 * h))


def fuel_fraction(p):
    """Mission fuel mass fraction m_f / MTOM and the intermediate quantities."""
    V = p["V_cr"] if p["V_cr"] is not None else p["M_cr"] * isa_speed_of_sound(p["h_cr"])
    CD0 = p["Cf_eq"] * p["Swet_S"]
    LD = 0.5 * math.sqrt(math.pi * p["AR"] * p["e"] / CD0)
    if p["LD_cruise"] is not None:
        LD = p["LD_cruise"]
    if p["eta_eng"] is not None:
        eta = p["eta_eng"]
        TSFC = V / (eta * p["e_f"])
    else:
        TSFC = p["TSFC"] if p["TSFC"] is not None else 22e-6 * p["BPR"] ** -0.19
        eta = V / (TSFC * p["e_f"])
    R_lost = (1 / 0.7) * LD * (p["h_cr"] + V ** 2 / (2 * G))
    R_eq = (p["R_nom"] + R_lost) * (1 + p["f_con"]) + 1.2 * p["R_div"] + p["t_E"] * V
    ff = 1 - math.exp(-R_eq / (eta * (p["e_f"] / G) * LD))
    return ff, dict(V_cr=V, LD=LD, CD0=CD0, TSFC=TSFC, eta_eng=eta,
                    R_lost=R_lost, R_eq=R_eq)


def nominal_range(p, ff):
    """Nominal range [m] flown with fuel fraction ff (inverse of fuel_fraction)."""
    _, d = fuel_fraction(p)
    R_eq = -math.log(1 - ff) * d["eta_eng"] * (p["e_f"] / G) * d["LD"]
    return (R_eq - 1.2 * p["R_div"] - p["t_E"] * d["V_cr"]) / (1 + p["f_con"]) - d["R_lost"]


def class1(p=INPUTS, oem=None):
    """Class I masses in kg.

    oem=None  -> m_OE = f_OE * MTOM (statistics, the WP1 estimate)
    oem=value -> m_OE fixed, e.g. the Class II OEM inside the WP3.4 iteration."""
    ff, d = fuel_fraction(p)
    if oem is None:
        if 1 - p["f_OE"] - ff <= 0:
            raise ValueError(f"no solution: f_OE + fuel fraction = {p['f_OE'] + ff:.3f} >= 1")
        mtom = p["payload"] / (1 - p["f_OE"] - ff)
        oem = p["f_OE"] * mtom
    else:
        mtom = (oem + p["payload"]) / (1 - ff)
    m_f = ff * mtom
    return dict(MTOM=mtom, OEM=oem, fuel=m_f, payload=p["payload"], fuel_fraction=ff,
                OEM_fraction=oem / mtom, MZFM=oem + p["payload"],
                V_fuel=m_f / p["rho_fuel"], derived=d)


# ================================================================== reporting
def print_report(p=INPUTS, oem=None):
    r = class1(p, oem)
    d = r["derived"]
    print("=" * 52)
    print(" CLASS I WEIGHT ESTIMATION  (ADSEE fuel fraction)")
    print("=" * 52)
    print(f"{'Item':<28}{'Mass [kg]':>12}{'% MTOM':>10}")
    print("-" * 52)
    for name, key in [("OEM", "OEM"), ("Payload", "payload"), ("Fuel", "fuel"), ("MTOM", "MTOM")]:
        print(f"{name:<28}{r[key]:>12,.0f}{100 * r[key] / r['MTOM']:>10.2f}")
    print(f"{'MZFM (OEM + payload)':<28}{r['MZFM']:>12,.0f}{100 * r['MZFM'] / r['MTOM']:>10.2f}")
    print(f"{'Fuel volume [m3]':<28}{r['V_fuel']:>12.2f}")
    diff = r["MTOM"] - p["MTOM_ref"]
    print(f"{'MTOM (WP1 Excel)':<28}{p['MTOM_ref']:>12,.0f}")
    print(f"{'Difference':<28}{diff:>12,.0f}{100 * diff / p['MTOM_ref']:>9.1f}%")

    print("\nDerived quantities")
    for k, unit, scale in [("V_cr", "m/s", 1), ("CD0", "-", 1), ("LD", "-", 1), ("TSFC", "g/(kN s)", 1e6), ("eta_eng", "-", 1),
                           ("R_lost", "km", 1e-3), ("R_eq", "km", 1e-3)]:
        print(f"  {k:<10}{d[k] * scale:>10.4f} {unit}")
    return r


if __name__ == "__main__":
    print_report()
