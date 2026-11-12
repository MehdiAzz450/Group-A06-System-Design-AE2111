"""
AE2111-I  WP3.3a  -  Class II drag estimation (component drag build-up)
Group A06 - reasonably sized passenger aircraft

Method: component drag build-up (ADSEE-II; Raymer Ch. 12):
    C_D0 = (1/S_ref) * sum(C_f,c * FF_c * IF_c * S_wet,c) + C_D,misc + C_D,L&P
    C_f   = k_lam * 1.328/sqrt(Re) + (1 - k_lam) * 0.455 / ((log10 Re)^2.58 (1 + 0.144 M^2)^0.65)
    Re    = min(rho V l / mu, 44.62 (l/k)^1.053 M^1.16)        (cut-off Reynolds number, transonic)
    FF_wing/tail = [1 + 0.6/(x/c)_m (t/c) + 100 (t/c)^4] [1.34 M^0.18 (cos Lambda_m)^0.28]
    FF_fus  = 1 + 60/f^3 + f/400,     FF_nac = 1 + 0.35/f,    f = l / d
    S_wet,wing/tail = S_exposed (1.977 + 0.52 t/c)
    C_D,upsweep = 3.83 u^2.5 A_max / S_ref   (u in rad)
    e = 1 / (pi A psi + 1/phi)
    C_D,wave = 20 (M - M_cr)^4,  M_cr = M_DD - (0.1/80)^(1/3)   (for M > M_cr)

Geometry is read from the Class II weight inputs (class2_weight.INPUTS), so both
estimates always use the same aircraft.
"""

import math

import class2_weight as c2

DRAG = dict(
    # ---- flight condition (cruise)
    h=9_448.8,               # m    FL310                                      [WP1]
    V=257.0,                 # m/s  cruise speed, M 0.85                        [WP1]
    M=0.85,                  # -                                               [WP1]
    M_DD=0.861,              # -    wing drag divergence Mach number           [WP2 3.6]
    k_skin=0.634e-5,         # m    skin roughness, smooth paint               [Raymer Table 12.5]
    # ---- laminar flow fractions
    lam_wing=0.10,           # -                                               [ASSUMED, check ADSEE-II]
    lam_tail=0.10,           # -                                               [ASSUMED, check ADSEE-II]
    lam_fus=0.05,            # -                                               [ASSUMED, check ADSEE-II]
    lam_nac=0.05,            # -                                               [ASSUMED, check ADSEE-II]
    # ---- airfoil position of maximum thickness
    xc_max_wing=0.37,        # -    SC(2)-0612                                 [ASSUMED, check airfoil data]
    xc_max_tail=0.30,        # -    NACA 0012-type tail section                [ASSUMED]
    tc_ht=0.12,              # -    horizontal tail thickness                  [ASSUMED]
    # ---- interference factors
    IF_wing=1.0,             # -    low wing with fillet                       [Raymer 12.5.5]
    IF_fus=1.0,              # -                                               [Raymer 12.5.5]
    IF_tail=1.04,            # -    conventional tail                          [Raymer 12.5.5]
    IF_nac=1.3,              # -    pylon nacelle within one diameter of wing  [Raymer 12.5.5]
    # ---- miscellaneous
    upsweep_deg=4.0,         # deg  tail-cone upsweep: angle between the fuselage centreline and the
                             #      line from the start of the tail cone to the tail endpoint
                             #      [reader Sec. 6.1.5 / Fig. 6.8 definition, Sec. 6.5 Step 6;
                             #       value as Example 6.8; Raymer Eq. 12.36 u]
    f_LP=0.03,               # -    leakage & protuberance, fraction of C_D0   [ASSUMED, Raymer 2-5 %]
    # ---- Oswald factor
    psi=0.0075,              # -    parasite drag due to lift                  [ADSEE]
    phi=0.97,                # -    span efficiency                            [ADSEE]
)


def isa(h):
    T = 288.15 - 0.0065 * h
    rho = 1.225 * (T / 288.15) ** 4.2559
    mu = 1.458e-6 * T ** 1.5 / (T + 110.4)                      # Sutherland
    return T, rho, mu


def chord_at(y, c_root, b, taper):
    return c_root * (1 - (1 - taper) * 2 * y / b)


def planform(S, AR, taper, sweep_c4_deg, x_c=0.25):
    """Span, root chord, MAC and sweep at chord fraction x_c of a trapezoidal planform."""
    b = math.sqrt(AR * S)
    c_r = 2 * S / (b * (1 + taper))
    mac = 2 / 3 * c_r * (1 + taper + taper ** 2) / (1 + taper)
    tan_x = math.tan(math.radians(sweep_c4_deg)) - 4 / AR * (x_c - 0.25) * (1 - taper) / (1 + taper)
    return b, c_r, mac, math.atan(tan_x)


def exposed_area(S, b, c_r, taper, w_body):
    """Planform area outside a body of width w_body (both halves)."""
    return S - w_body / 2 * (c_r + chord_at(w_body / 2, c_r, b, taper))


def skin_friction(Re_actual, l, lam, d):
    Re_cut = 44.62 * (l / d["k_skin"]) ** 1.053 * d["M"] ** 1.16
    Re = min(Re_actual, Re_cut)
    cf_lam = 1.328 / math.sqrt(Re)
    cf_turb = 0.455 / (math.log10(Re) ** 2.58 * (1 + 0.144 * d["M"] ** 2) ** 0.65)
    return lam * cf_lam + (1 - lam) * cf_turb, Re


def ff_lifting(tc, xc_max, sweep_m, M):
    return (1 + 0.6 / xc_max * tc + 100 * tc ** 4) * (1.34 * M ** 0.18 * math.cos(sweep_m) ** 0.28)


def class2_drag(p=c2.INPUTS, d=DRAG):
    """Zero-lift drag build-up, Oswald factor and cruise polar for the aircraft in p."""
    _, rho, mu = isa(d["h"])
    Re_m = rho * d["V"] / mu
    S = p["S_w"]
    comps = {}

    # wing
    b, c_r, mac, sw_m = planform(S, p["AR"], p["taper"], p["sweep_c4_deg"], d["xc_max_wing"])
    S_exp = exposed_area(S, b, c_r, p["taper"], p["d_fus"])
    comps["Wing"] = dict(l=mac, S_wet=S_exp * (1.977 + 0.52 * p["tc_root"]), lam=d["lam_wing"],
                         FF=ff_lifting(p["tc_root"], d["xc_max_wing"], sw_m, d["M"]), IF=d["IF_wing"])
    # horizontal tail (taper 0.4, from WP3 3.4.2)
    b_h, c_rh, mac_h, sw_mh = planform(p["S_ht"], p["AR_ht"], 0.4, p["sweep_ht_deg"], d["xc_max_tail"])
    w_body = p["F_w"] if p.get("F_w") is not None else c2.derived(p)["F_w"]
    S_exp_h = exposed_area(p["S_ht"], b_h, c_rh, 0.4, w_body)
    comps["Horizontal tail"] = dict(l=mac_h, S_wet=S_exp_h * (1.977 + 0.52 * d["tc_ht"]), lam=d["lam_tail"],
                                    FF=ff_lifting(d["tc_ht"], d["xc_max_tail"], sw_mh, d["M"]), IF=d["IF_tail"])
    # vertical tail (taper 0.35; one half-wing of a planform with area 2 S_v and aspect ratio 2 A_v)
    _, _, mac_v, sw_mv = planform(2 * p["S_vt"], 2 * p["AR_vt"], 0.35, p["sweep_vt_deg"], d["xc_max_tail"])
    comps["Vertical tail"] = dict(l=mac_v, S_wet=p["S_vt"] * (1.977 + 0.52 * p["tc_root_vt"]), lam=d["lam_tail"],
                                  FF=ff_lifting(p["tc_root_vt"], d["xc_max_tail"], sw_mv, d["M"]), IF=d["IF_tail"])
    # fuselage
    f_fus = p["l_fus"] / p["d_fus"]
    comps["Fuselage"] = dict(l=p["l_fus"], lam=d["lam_fus"], IF=d["IF_fus"], FF=1 + 60 / f_fus ** 3 + f_fus / 400,
                             S_wet=c2.fuselage_wetted_area(p["d_fus"], p["l_fus"], p["nose_slender"], p["tail_slender"]))
    # nacelles
    f_nac = p["l_nac"] / p["d_nac"]
    comps["Nacelles"] = dict(l=p["l_nac"], lam=d["lam_nac"], IF=d["IF_nac"], FF=1 + 0.35 / f_nac,
                             S_wet=p["N_en"] * math.pi * (p["d_nac"] + p["d_nac_exit"]) / 2 * p["l_nac"])

    for c in comps.values():
        c["Cf"], c["Re"] = skin_friction(Re_m * c["l"], c["l"], c["lam"], d)
        c["CD0"] = c["Cf"] * c["FF"] * c["IF"] * c["S_wet"] / S

    CD0_comp = sum(c["CD0"] for c in comps.values())
    CD_upsweep = 3.83 * math.radians(d["upsweep_deg"]) ** 2.5 * (math.pi / 4 * p["d_fus"] ** 2) / S
    CD0_noLP = CD0_comp + CD_upsweep
    CD0 = CD0_noLP / (1 - d["f_LP"])                            # L&P as a fraction of the total
    CD_LP = CD0 - CD0_noLP
    e = 1 / (math.pi * p["AR"] * d["psi"] + 1 / d["phi"])
    M_cr = d["M_DD"] - (0.1 / 80) ** (1 / 3)
    CD_wave = 20 * (d["M"] - M_cr) ** 4 if d["M"] > M_cr else 0.0
    S_wet_total = sum(c["S_wet"] for c in comps.values())
    return dict(components=comps, CD0=CD0, CD0_components=CD0_comp, CD_upsweep=CD_upsweep, CD_LP=CD_LP,
                e=e, CD_wave=CD_wave, M_cr=M_cr, S_wet=S_wet_total, Swet_S=S_wet_total / S,
                LD_max=0.5 * math.sqrt(math.pi * p["AR"] * e / CD0), rho=rho, q=0.5 * rho * d["V"] ** 2,
                Re_per_m=Re_m)


def cruise_point(r, p, m_cruise):
    """Cruise C_L, C_D, L/D and drag [N] of the complete aircraft and of the wing alone."""
    CL = m_cruise * 9.81 / (r["q"] * p["S_w"])
    CDi = CL ** 2 / (math.pi * p["AR"] * r["e"])
    CD = r["CD0"] + CDi + r["CD_wave"]
    CD_wing = r["components"]["Wing"]["CD0"] + CDi + r["CD_wave"]
    return dict(CL=CL, CD=CD, LD=CL / CD, D=CD * r["q"] * p["S_w"], D_wing=CD_wing * r["q"] * p["S_w"])


def print_report(p=c2.INPUTS, d=DRAG):
    r = class2_drag(p, d)
    print("=" * 78)
    print(" CLASS II DRAG ESTIMATION  (component build-up, cruise)")
    print("=" * 78)
    print(f"{'Component':<18}{'l [m]':>8}{'Re [-]':>11}{'Cf [-]':>10}{'FF':>7}{'IF':>6}"
          f"{'S_wet [m2]':>12}{'C_D0 [cts]':>11}")
    for k, c in r["components"].items():
        print(f"{k:<18}{c['l']:>8.2f}{c['Re']:>11.3e}{c['Cf']:>10.5f}{c['FF']:>7.3f}{c['IF']:>6.2f}"
              f"{c['S_wet']:>12.1f}{1e4 * c['CD0']:>11.1f}")
    print(f"{'Fuselage upsweep':<72}{1e4 * r['CD_upsweep']:>6.1f}")
    print(f"{'Leakage & protuberance':<72}{1e4 * r['CD_LP']:>6.1f}")
    print(f"{'C_D0 total':<72}{1e4 * r['CD0']:>6.1f}")
    print(f"S_wet/S = {r['Swet_S']:.2f},  e = {r['e']:.3f},  (L/D)_max = {r['LD_max']:.2f},  "
          f"C_D,wave = {1e4 * r['CD_wave']:.1f} cts (M_cr = {r['M_cr']:.3f})")
    return r


if __name__ == "__main__":
    print_report()
