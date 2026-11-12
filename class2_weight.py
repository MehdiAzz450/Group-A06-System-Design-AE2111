"""
AE2111-I  WP3.3b  -  Class II weight estimation (Raymer, Cargo/Transport equations)
Group A06 - reasonably sized passenger aircraft

Method: D.P. Raymer, "Aircraft Design: A Conceptual Approach", Ch. 15,
        Cargo/Transport weights, Eqs. 15.25 - 15.44.
Supplemented by E. Torenbeek, "Synthesis of Subsonic Airplane Design", Ch. 8:
        Eq. 8-44 (furnishings) and Table 8-13 (operational items), because
        Raymer's cargo/transport set does not cover a passenger cabin or OEM items.
        Raymer 15.44 (handling gear) is left out: it is inside Torenbeek's furnishing group.
        Raymer Sec. 15.4 calibration factors per group (see calibrate_raymer.py).
        Torenbeek Eq. 8-12 (wing) is only printed as a cross-check, never summed.
Raymer's equations are in imperial units (lb, ft, ft^2, in, gal, ft^3).
All inputs below are given in SI and converted inside `raymer_weights()`.

Usage
-----
    python class2_weight.py                 -> prints the breakdown for the inputs below
    from class2_weight import class2, INPUTS
    res = class2({**INPUTS, "MTOM": 99000})  -> reuse in the WP3.4 Class I / II iteration

Every input carries a comment with its source:
    [WP1] / [WP2] / [WP3.1] / [WP3.2]  -> value from your reports / Excel
    [WP3 x.y]                          -> section / drawing of the WP3 report (latest)
    [ASSUMED]                          -> typical value, to be justified in the report
    [PLACEHOLDER]                      -> still being calculated by the team, replace!
"""

import math

# ------------------------------------------------------------------ unit conversions
KG2LB = 2.20462262
M2FT = 3.28083990
M2IN = 39.3700787
M22FT2 = M2FT ** 2
M32FT3 = M2FT ** 3
M32GAL = 264.172052          # US gallons
MS2FTS = M2FT
G = 9.80665

# ================================================================== INPUTS (SI)
INPUTS = dict(
    # ---- masses & loads
    MTOM=100_999.0,          # kg   design gross mass W_dg                     [WP2 Table 5.1]
    OEM_class1=50_500.0,     # kg   Class I OEM, for comparison                [WP2 Table 5.1]
    f_MLM=0.776,             # -    landing mass fraction (W_l = f_MLM*MTOM)   [WP1/WP2]
    n_limit=2.5,             # -    limit manoeuvring load factor: 2.1+24000/(W+10000)=2.21 < 2.5 -> 2.5 [CS 25.337(b)]; n_ult = 1.5 n [CS 25.303]
    N_gear=None,             # -    None -> derived from the shock absorber energy balance below
    V_sink=3.05,             # m/s  limit descent velocity at design landing weight [CS 25.473(a)(2)(i), 10 ft/s]
    stroke=0.40,             # m    shock absorber stroke, reader range 0.30-0.50 m  [reader Sec. 9.2; WP3 gear drawing]
    eta_oleo=0.80,           # -    oleo-pneumatic shock absorber efficiency   [ASSUMED, typical value, Raymer Ch. 11]

    # ---- wing
    S_w=119.33,              # m2   wing area                                  [WP2 Table 5.2]
    b_w=35.398,              # m    span                                       [WP2]
    AR=10.5,                 # -                                               [WP2]
    taper=0.29,              # -                                               [WP2]
    sweep_c4_deg=30.77,      # deg  quarter-chord sweep                        [WP2]
    tc_root=0.12,            # -    SC(2)-0612                                 [WP2]
    wing_method="raymer",    # -    "raymer" (15.25) or "torenbeek" (8-12) in the OEM sum
    wing_spoilers=True,      # -    spoilers / speed brakes: +2 % (Torenbeek 8-12) [ASSUMED]
    N_en_wing=2,             # -    wing-mounted engines: -5 % (2), -10 % (4)   [WP3.1]
    # wing-mounted control surfaces (for S_csw)
    flap_eta=(0.114, 0.621), # -    Fowler flap spanwise extent                [WP2 Table 4.7]
    flap_cf_c=0.25,          # -    flap chord ratio                           [WP2]
    ail_eta=(0.77, 0.95),    # -    aileron spanwise extent                    [WP2 Table 5.2]
    ail_ca_c=0.275,          # -    aileron chord ratio                        [WP2 Table 5.2]

    # ---- fuselage
    l_fus=44.525,            # m    total fuselage length                      [WP3 drawing]
    d_fus=4.107,             # m    outer diameter (= width = depth)           [WP3.2]
    l_cabin=32.4,            # m                                               [WP3.2]
    d_fus_inner=3.85,        # m                                               [WP3.2]
    nose_slender=2.0,        # -    l_nc / d_fus  (l_nc = 8.00 m)              [WP3 drawing]
    tail_slender=3.0,        # -    l_tc / d_fus  (l_tc = 12.00 m)             [WP3 drawing]
    l_radome_tailcap=1.5,    # m    excluded from structural length            [ASSUMED]
    K_door=1.12,             # -    two side cargo doors (fwd + aft hold)      [Raymer 15.28 K_door definition]
    l_cockpit=3.0,           # m    cockpit length added to pressurised length [ASSUMED]

    # ---- empennage
    S_ht=15.648,             # m2   V_h = 0.70                                 [WP3 3.4.2, drawing]
    AR_ht=4.0,               # -                                               [WP3.2]
    sweep_ht_deg=35.0,       # deg                                             [WP3.2]
    Se_Sht=0.30,             # -    elevator / HT area                         [ASSUMED]
    K_uht=1.0,               # -    1.143 if all-moving HT                     [ASSUMED]
    S_vt=22.455,             # m2   V_v = 0.105                                [WP3 3.4.3, drawing]
    AR_vt=1.6,               # -                                               [WP3.2]
    sweep_vt_deg=40.0,       # deg                                             [WP3.2]
    tc_root_vt=0.12,         # -                                               [ASSUMED]
    Sr_Svt=0.30,             # -    rudder / VT area                           [ASSUMED]
    Ht_Hv=0.0,               # -    0 = conventional tail, 1 = T-tail          [WP3.2c conventional]
    x_ac_tail=40.073,        # m    tail aerodynamic centre, 0.9 l_fus         [WP3 drawing]
    L_t=19.751,              # m    tail arm (aft CG = 25 % MAC to tail a.c.)  [WP3 drawing]
    F_w=2.33,                # m    fuselage width at the HT a.c.; None -> cone [WP3 drawing]

    # ---- landing gear
    L_m=2.326,               # m    main gear strut length                     [WP3 drawing, sheet 5; PLACEHOLDER until gear design final]
    L_n=2.478,               # m    nose gear strut length                     [WP3 drawing, sheet 5; PLACEHOLDER until gear design final]
    N_mw=4,                  # -    main wheels                                [WP3.2]
    N_mss=2,                 # -    main gear shock struts                     [WP3.2]
    N_nw=2,                  # -    nose wheels                                [WP3.2]
    K_mp=1.0, K_np=1.0,      # -    1.126 / 1.15 if kneeling gear              [ASSUMED]
    K_Lg=1.0,                # -    1.12 if fuselage-mounted main gear         [ASSUMED wing-mounted]
    V_stall=66.2,            # m/s  V_SR0 at max landing mass                  [WP2 compliance matrix]

    # ---- propulsion (PW2037)
    N_en=2,                  # -                                               [WP3.1]
    m_en=7300 / KG2LB,       # kg   dry engine mass, PW2037 (7300 lb)          [WP3.1]
    T_max_en=165.5e3,        # N    take-off thrust per engine                 [WP3 2.4]
    l_en=146.8 / M2IN,       # m    engine length (146.8 in)                   [WP3.1]
    d_en=84.8 / M2IN,        # m    nominal engine diameter (84.8 in)          [WP3.1]
    l_nac=7.0988,            # m    nacelle length excl. cone, l_n (Type B)    [WP3.1]
    d_nac=2.4780,            # m    maximum nacelle diameter D_n               [WP3.1]
    d_nac_exit=1.6521,       # m    fan cowling exit diameter D_ef             [WP3.1]
    K_ng=1.017,              # -    pylon-mounted nacelle                      [Raymer]
    K_tr=1.18,               # -    thrust reversers                           [WP3.1]
    x_engine_front=13.48,    # m    engine inlet x from nose                   [WP3 drawing, sheet 6: front spar - 0.85 l_n, reader 8.3.2]

    # ---- fuel system
    V_fuel=39.43,            # m3   total fuel volume                          [WP2 Table 5.1]
    integral_fraction=1.0,   # -    fraction of volume in integral tanks       [ASSUMED]
    N_t=3,                   # -    tanks: 2 wing + centre                     [ASSUMED]

    # ---- systems & payload
    N_pax=180,               # -                                               [WP3 3.2.3]
    N_flightcrew=2,          # -                                               [ASSUMED]
    N_cabincrew=4,           # -    1 per 50 passenger seats, 180/50 -> 4      [EU Air Ops ORO.CC.100]
    m_cargo_max=3660.0,      # kg   hold luggage 3600 + cargo 60               [WP3 3.2.4]
    m_uld=6 * 72.0,          # kg   six empty AKH containers                   [WP3 3.2.4]
    N_f=6,                   # -    functions performed by controls            [ASSUMED, Raymer 4-7]
    N_m=1,                   # -    mechanical functions                       [ASSUMED, Raymer 0-2]
    R_zbar=0.39,             # -    non-dim. yaw radius of gyration            [ASSUMED, Roskam]
    m_apu_uninst=150.0,      # kg                                              [ASSUMED]
    R_kva=50.0,              # kVA                                             [ASSUMED, Raymer 40-60]
    N_gen=2,                 # -    = N_en                                     [ASSUMED]
    m_uav=500.0,             # kg   uninstalled avionics (~1100 lb)            [ASSUMED, Raymer 800-1400 lb]
    payload_max=18_960.0,    # kg   max payload (for MZFM)                     [WP1]

    # ---- furnishings & operational items (Torenbeek, Ch. 8)
    # Raymer's cargo/transport furnishings eq. (15.41) is fitted on military cargo
    # aircraft and omits passenger seats, galleys and lavatories. For a passenger
    # aircraft Torenbeek Eq. 8-44 (furnishings = 0.196 * MZFM^0.91, kg) is used.
    furnishings_method="torenbeek",   # "torenbeek" or "raymer"
    long_range=True,         # -    Torenbeek Table 8-13 long-range values     [WP1: 9545 km]
    overwater=True,          # -    extended overwater (life jackets)          [WP2: transatlantic]

    # ---- calibration factors per group (Raymer Sec. 15.4, "fudge factors").
    # 1.0 = pure statistical result.
    # Factors = actual / Raymer-calculated group mass, mean of B737-200 and A300B2
    # (actual masses: ADSEE-II weight lecture slide 41, source Torenbeek; geometry: Jenkinson et al.).
    # Computed with calibrate_raymer.py on THIS file. 'surface_controls' multiplies Raymer 15.35 only.
    # Rerun calibrate_raymer.py after ANY change to the equations or to inputs shared with
    # the reference aircraft (e.g. K_door, N_gear), and copy the new means here.
    fudge=dict(wing=1.405, tail=1.163, fuselage=1.303, landing_gear=1.064, nacelle=1.119,
               surface_controls=2.257),
)

# Validation data: group mass as % of MTOW for the twin-engine jet transports on
# ADSEE-II weight lecture slide 41 ("Some reference data to compare your results with",
# source Torenbeek). Order: wing, tail, fuselage, landing gear, surface controls, nacelle.
# None = not given on the slide.
SLIDE41_TWIN = {
    "VFW-Fokker 614":   (14.1, 2.74, 12.8, 3.45, 1.82, 2.37),
    "Fokker F-28/1000": (11.3, 2.46, 10.8, 4.24, 2.13, 1.28),
    "Fokker F-28/2000": (11.3, 2.46, 11.8, 4.24, 2.15, 1.28),
    "Fokker F-28/5000": (11.6, 2.31, 9.95, 3.90, 2.35, 1.20),
    "Fokker F-28/6000": (11.6, 2.31, 10.8, 3.94, 2.36, 1.20),
    "BAC 1-11/300":     (11.1, 2.72, 11.2, 3.29, 1.76, None),
    "BAC 1-11/400":     (11.1, 2.78, 11.2, 3.33, 1.39, None),
    "McD DC-9/10":      (10.3, 2.87, 12.2, 4.00, 1.38, 1.55),
    "Boeing 737-100C":  (10.2, 2.76, 12.7, 3.77, 1.62, None),
    "Boeing 737-200":   (10.6, 2.72, 12.1, 4.35, 2.35, 1.39),
    "Caravelle XIR":    (13.4, 1.77, 10.5, 4.63, 1.87, 1.43),
    "Airbus A300B2":    (14.5, 1.95, 11.8, 4.47, 1.94, 2.32),
}
SLIDE41_GROUPS = ("wing", "tail", "fuselage", "landing_gear", "surface_controls", "nacelle")


def slide41_range(group):
    """(min, mean, max) % MTOW of one group over the twin-engine jets of slide 41."""
    i = SLIDE41_GROUPS.index(group)
    v = [row[i] for row in SLIDE41_TWIN.values() if row[i] is not None]
    return min(v), sum(v) / len(v), max(v)


# ================================================================== geometry helpers
def wing_strip_area_ratio(eta_a, eta_b, taper):
    """Area between spanwise stations eta_a..eta_b (both halves) / S, straight-tapered wing.
    Same relation as WP2 Eq. 4.16."""
    return (eta_b - eta_a) * (2 - (1 - taper) * (eta_a + eta_b)) / (1 + taper)


def fuselage_wetted_area(d, l_fus, nose_slender, tail_slender):
    """ADSEE / Torenbeek fuselage wetted area (nose ellipsoid + cylinder + tail cone)."""
    l_n, l_t = nose_slender * d, tail_slender * d
    l_c = l_fus - l_n - l_t
    return (math.pi * d / 4) * (1 / (3 * l_n ** 2) * ((4 * l_n ** 2 + d ** 2 / 4) ** 1.5 - d ** 3 / 8)
                                - d + 4 * l_c + 2 * math.sqrt(l_t ** 2 + d ** 2 / 4))


def derived(p):
    """Geometric / derived quantities in SI."""
    d = {}
    d["sweep"] = math.radians(p["sweep_c4_deg"])
    d["S_flap"] = p["flap_cf_c"] * wing_strip_area_ratio(*p["flap_eta"], p["taper"]) * p["S_w"]
    d["S_ail"] = p["ail_ca_c"] * wing_strip_area_ratio(*p["ail_eta"], p["taper"]) * p["S_w"]
    d["S_csw"] = d["S_flap"] + d["S_ail"]
    d["S_e"] = p["Se_Sht"] * p["S_ht"]
    d["S_r"] = p["Sr_Svt"] * p["S_vt"]
    d["S_cs"] = d["S_csw"] + d["S_e"] + d["S_r"]
    d["b_ht"] = math.sqrt(p["AR_ht"] * p["S_ht"])
    # fuselage width at the HT: linear tail-cone taper evaluated at the tail a.c.
    l_tc = p["tail_slender"] * p["d_fus"]
    d["F_w"] = (p["F_w"] if p.get("F_w") is not None
                else max(p["d_fus"] * (p["l_fus"] - p["x_ac_tail"]) / l_tc, 0.0))
    d["S_f"] = fuselage_wetted_area(p["d_fus"], p["l_fus"], p["nose_slender"], p["tail_slender"])
    d["L_struct"] = p["l_fus"] - p["l_radome_tailcap"]
    d["K_ws"] = 0.75 * ((1 + 2 * p["taper"]) / (1 + p["taper"])) * p["b_w"] * math.tan(d["sweep"]) / d["L_struct"]
    d["l_nac"] = p["l_nac"]
    d["w_nac"] = p["d_nac"]
    # nacelle wetted area per nacelle: frustum between D_n and D_ef over l_n
    d["S_n"] = math.pi * (p["d_nac"] + p["d_nac_exit"]) / 2 * p["l_nac"]
    d["V_pr"] = math.pi / 4 * p["d_fus_inner"] ** 2 * (p["l_cabin"] + p["l_cockpit"])
    d["L_ec"] = p["N_en"] * (p["x_engine_front"] - p["l_cockpit"])  # total, multi-engine
    d["L_a"] = p["x_engine_front"] - p["l_cockpit"]                  # electrical routing
    d["N_crew"] = p["N_flightcrew"] + p["N_cabincrew"]
    d["N_p"] = p["N_pax"] + d["N_crew"]
    d["n_ult"] = 1.5 * p["n_limit"]
    # landing gear load factor: kinetic energy of the sink rate absorbed over the stroke,
    # N_gear = V^2 / (2 g eta S)  (Raymer Ch. 11 stroke equation, tyre deflection neglected)
    d["N_gear"] = (p["N_gear"] if p.get("N_gear") is not None
                   else p["V_sink"] ** 2 / (2 * G * p["eta_oleo"] * p["stroke"]))
    d["N_l"] = 1.5 * d["N_gear"]
    return d


# ================================================================== Raymer equations
def raymer_weights(p):
    """Component weights in kg, Raymer Eqs. 15.25-15.44 (evaluated in imperial units)."""
    d = derived(p)
    # --- imperial conversions
    W_dg = p["MTOM"] * KG2LB
    W_l = p["f_MLM"] * p["MTOM"] * KG2LB
    Nz, Nl = d["n_ult"], d["N_l"]
    S_w = p["S_w"] * M22FT2
    S_csw = d["S_csw"] * M22FT2
    S_ht, S_vt = p["S_ht"] * M22FT2, p["S_vt"] * M22FT2
    S_e = d["S_e"] * M22FT2
    S_cs = d["S_cs"] * M22FT2
    L_t = p["L_t"] * M2FT
    K_y, K_z = 0.3 * L_t, L_t                                   # Raymer approximations
    F_w, B_h = d["F_w"] * M2FT, d["b_ht"] * M2FT
    L = d["L_struct"] * M2FT
    D = p["d_fus"] * M2FT
    S_f = d["S_f"] * M22FT2
    L_m, L_n = p["L_m"] * M2IN, p["L_n"] * M2IN
    V_stall = p["V_stall"] * MS2FTS
    W_en = p["m_en"] * KG2LB
    N_Lt, N_w, S_n = d["l_nac"] * M2FT, d["w_nac"] * M2FT, d["S_n"] * M22FT2
    W_ec = 2.331 * W_en ** 0.901 * 1.0 * p["K_tr"]              # K_p = 1.0 (jet)
    L_ec, L_a = d["L_ec"] * M2FT, d["L_a"] * M2FT
    V_t = p["V_fuel"] * M32GAL
    V_i = p["integral_fraction"] * V_t
    V_p = 0.0                                                    # no self-sealing tanks
    B_w, L_f = p["b_w"] * M2FT, p["l_fus"] * M2FT
    I_y = (((B_w + L_f) / 2) ** 2 * W_dg * p["R_zbar"] ** 2) / 4  # yaw, lb ft^2 (Roskam I_zz, mass in lb)
    W_uav = p["m_uav"] * KG2LB
    W_c = p["m_cargo_max"] * KG2LB
    V_pr = d["V_pr"] * M32FT3
    N_en, N_f, N_m = p["N_en"], p["N_f"], p["N_m"]
    N_c = p["N_flightcrew"]
    lam, A, tc = p["taper"], p["AR"], p["tc_root"]
    cosL = math.cos(d["sweep"])

    w = {}
    # --- structures
    w["Wing"] = (0.0051 * (W_dg * Nz) ** 0.557 * S_w ** 0.649 * A ** 0.5 * tc ** -0.4
                 * (1 + lam) ** 0.1 * cosL ** -1.0 * S_csw ** 0.1)                          # 15.25
    w["Horizontal tail"] = (0.0379 * p["K_uht"] * (1 + F_w / B_h) ** -0.25 * W_dg ** 0.639
                            * Nz ** 0.10 * S_ht ** 0.75 * L_t ** -1.0 * K_y ** 0.704
                            * math.cos(math.radians(p["sweep_ht_deg"])) ** -1.0
                            * p["AR_ht"] ** 0.166 * (1 + S_e / S_ht) ** 0.1)                   # 15.26
    w["Vertical tail"] = (0.0026 * (1 + p["Ht_Hv"]) ** 0.225 * W_dg ** 0.556 * Nz ** 0.536
                          * L_t ** -0.5 * S_vt ** 0.5 * K_z ** 0.875
                          * math.cos(math.radians(p["sweep_vt_deg"])) ** -1
                          * p["AR_vt"] ** 0.35 * p["tc_root_vt"] ** -0.5)                       # 15.27
    w["Fuselage"] = (0.3280 * p["K_door"] * p["K_Lg"] * (W_dg * Nz) ** 0.5 * L ** 0.25
                     * S_f ** 0.302 * (1 + d["K_ws"]) ** 0.04 * (L / D) ** 0.10)               # 15.28
    w["Main landing gear"] = (0.0106 * p["K_mp"] * W_l ** 0.888 * Nl ** 0.25 * L_m ** 0.4
                              * p["N_mw"] ** 0.321 * p["N_mss"] ** -0.5 * V_stall ** 0.1)      # 15.29
    w["Nose landing gear"] = (0.032 * p["K_np"] * W_l ** 0.646 * Nl ** 0.2 * L_n ** 0.5
                              * p["N_nw"] ** 0.45)                                             # 15.30
    w["Nacelle group"] = (0.6724 * p["K_ng"] * N_Lt ** 0.10 * N_w ** 0.294 * Nz ** 0.119
                          * W_ec ** 0.611 * N_en ** 0.984 * S_n ** 0.224)                      # 15.31
    # --- propulsion
    w["Engines (dry)"] = N_en * W_en
    w["Engine controls"] = 5.0 * N_en + 0.80 * L_ec                                           # 15.32
    w["Starter (pneumatic)"] = 49.19 * (N_en * W_en / 1000) ** 0.541                           # 15.33
    w["Fuel system"] = (2.405 * V_t ** 0.606 * (1 + V_i / V_t) ** -1.0
                        * (1 + V_p / V_t) * p["N_t"] ** 0.5)                                    # 15.34
    # --- fixed equipment
    w["Flight controls"] = (145.9 * N_f ** 0.554 * (1 + N_m / N_f) ** -1.0 * S_cs ** 0.20
                            * (I_y * 1e-6) ** 0.07)                                             # 15.35
    w["APU installed"] = 2.2 * p["m_apu_uninst"] * KG2LB                                      # 15.36
    w["Instruments"] = 4.509 * 1.0 * 1.0 * N_c ** 0.541 * N_en * (L_f + B_w) ** 0.5            # 15.37
    w["Hydraulics"] = 0.2673 * N_f * (L_f + B_w) ** 0.937                                      # 15.38
    w["Electrical"] = 7.291 * p["R_kva"] ** 0.782 * L_a ** 0.346 * p["N_gen"] ** 0.10         # 15.39
    w["Avionics"] = 1.73 * W_uav ** 0.983                                                      # 15.40
    w["Furnishings"] = 0.0577 * N_c ** 0.1 * W_c ** 0.393 * S_f ** 0.75                        # 15.41
    w["Air conditioning"] = 62.36 * d["N_p"] ** 0.25 * (V_pr / 1000) ** 0.604 * W_uav ** 0.10  # 15.42
    w["Anti-ice"] = 0.002 * W_dg                                                               # 15.43
    # Raymer 15.44 "handling gear" (3.0e-4 W_dg) is NOT included: cargo restraints and
    # handling provisions are part of Torenbeek's furnishing & equipment group (Table 8-12),
    # which Eq. 8-44 covers. Keeping both would count the same hardware twice.

    w = {k: v / KG2LB for k, v in w.items()}
    # optional calibration ("fudge") factors per group, Raymer Sec. 15.4
    fudge = p.get("fudge", {})
    for g, items in GROUPS.items():
        if g == "fixed_equipment":
            continue                     # no group factor on fixed equipment as a whole
        for k in items:
            if k in w:
                w[k] *= fudge.get(g, 1.0)
    # slide 41 calibrates "surface controls" separately -> Raymer 15.35 flight controls only
    w["Flight controls"] *= fudge.get("surface_controls", 1.0)
    return w, d


GROUPS = {
    "wing": ["Wing"],
    "tail": ["Horizontal tail", "Vertical tail"],
    "fuselage": ["Fuselage"],
    "landing_gear": ["Main landing gear", "Nose landing gear"],
    "nacelle": ["Nacelle group"],
    "propulsion": ["Engines (dry)", "Engine controls", "Starter (pneumatic)", "Fuel system"],
    "fixed_equipment": ["Flight controls", "APU installed", "Instruments", "Hydraulics",
                        "Electrical", "Avionics", "Furnishings", "Air conditioning",
                        "Anti-ice"],
}


def operational_items(p, d):
    """Torenbeek Table 8-13, operational items for transport aircraft (kg)."""
    n = p["N_pax"]
    v_fuel_l = p["V_fuel"] * 1000.0
    return {
        "Crew + baggage": 93.0 * p["N_flightcrew"] + 68.0 * p["N_cabincrew"],
        "Passenger cabin supplies": (8.62 if p["long_range"] else 6.35) * n,
        "Potable water & toilet chem.": (2.95 if p["long_range"] else 1.36) * n,
        "Safety equipment": (3.4 if p["overwater"] else 0.907) * n,
        "Residual fuel": 0.151 * v_fuel_l ** (2.0 / 3.0),
        "Cargo containers (ULD)": p.get("m_uld", 0.0),
    }


def torenbeek_wing(p, mzfm):
    """Wing group mass in kg, Torenbeek Eq. 8-12 (SI: kg, m), transport category.

        W_w / W_ZF = k_w * b_s^0.75 * (1 + sqrt(b_ref / b_s)) * n_ult^0.55
                     * ((b_s / t_r) / (W_ZF / S))^0.30
        k_w = 6.67e-3 (W_TO > 5670 kg), b_ref = 1.905 m, b_s = b / cos(sweep_c/2)

    Includes high-lift devices and ailerons. t_r is taken at the centreline root
    chord of the trapezoidal planform."""
    lam = p["taper"]
    tan_c2 = math.tan(math.radians(p["sweep_c4_deg"])) - (1 / p["AR"]) * (1 - lam) / (1 + lam)
    b_s = p["b_w"] * math.sqrt(1 + tan_c2 ** 2)
    c_root = 2 * p["S_w"] / (p["b_w"] * (1 + lam))
    t_r = p["tc_root"] * c_root
    n_ult = 1.5 * p["n_limit"]
    frac = (6.67e-3 * b_s ** 0.75 * (1 + math.sqrt(1.905 / b_s)) * n_ult ** 0.55
            * ((b_s / t_r) / (mzfm / p["S_w"])) ** 0.30)
    corr = 1.02 if p["wing_spoilers"] else 1.0
    corr *= {2: 0.95, 4: 0.90}.get(p["N_en_wing"], 1.0)
    if p["K_Lg"] != 1.0:                                          # main gear not on the wing
        corr *= 0.95
    return corr * frac * mzfm


def class2(p=INPUTS, tol=1.0, max_iter=50):
    """Full Class II result: components, groups, empty mass and OEM (kg).

    Torenbeek's furnishings and wing equations are functions of MZFM = OEM + max
    payload, so the OEM is found by fixed-point iteration."""
    comp, d = raymer_weights(p)
    raymer_furn = comp.pop("Furnishings")
    raymer_wing = comp["Wing"]
    ops = operational_items(p, d)
    m_ops = sum(ops.values())
    f_wing = p.get("fudge", {}).get("wing", 1.0)

    def assemble(mzfm):
        c = dict(comp)
        if p["wing_method"] == "torenbeek":
            c["Wing"] = f_wing * torenbeek_wing(p, mzfm)
        if p["furnishings_method"] == "torenbeek":
            c["Furnishings (Torenbeek 8-44)"] = 0.196 * mzfm ** 0.91
        else:
            c["Furnishings (Raymer 15.41)"] = raymer_furn
        return c

    oem = sum(comp.values()) + raymer_furn + m_ops
    for _ in range(max_iter):
        new = sum(assemble(oem + p["payload_max"]).values()) + m_ops
        converged = abs(new - oem) < tol
        oem = new
        if converged:
            break
    c = assemble(oem + p["payload_max"])
    m_empty = sum(c.values())
    oem = m_empty + m_ops
    mzfm = oem + p["payload_max"]
    groups = {g: sum(c[k] for k in c if k in items or (g == "fixed_equipment" and k.startswith("Furnishings")))
              for g, items in GROUPS.items()}
    return dict(components=c, groups=groups, derived=d, empty=m_empty, ops=ops,
                OEM=oem, MZFM=mzfm, raymer_furnishings=raymer_furn,
                wing_raymer=raymer_wing, wing_torenbeek=f_wing * torenbeek_wing(p, mzfm))


# ================================================================== reporting
def print_report(p=INPUTS):
    r = class2(p)
    mtom = p["MTOM"]
    print("=" * 60)
    print(" CLASS II WEIGHT ESTIMATION  (Raymer, cargo/transport)")
    print("=" * 60)
    print(f"{'Component':<30}{'Mass [kg]':>10}{'% MTOM':>10}{'% OEM':>10}")
    print("-" * 60)
    comps = r["components"]
    for g, items in GROUPS.items():
        print(f"[{g.replace('_', ' ')}]")
        names = [c for c in comps if c in items or (g == "fixed_equipment" and c.startswith("Furnishings"))]
        for c in names:
            m = comps[c]
            print(f"  {c:<28}{m:>10,.0f}{100 * m / mtom:>10.2f}{100 * m / r['OEM']:>10.2f}")
    print("-" * 60)
    print(f"{'Empty mass':<30}{r['empty']:>10,.0f}{100 * r['empty'] / mtom:>10.2f}")
    print("[operational items, Torenbeek Table 8-13]")
    for k, m in r["ops"].items():
        print(f"  {k:<28}{m:>10,.0f}{100 * m / mtom:>10.2f}{100 * m / r['OEM']:>10.2f}")
    print("-" * 60)
    print(f"{'OEM (Class II)':<30}{r['OEM']:>10,.0f}{100 * r['OEM'] / mtom:>10.2f}")
    print(f"{'MZFM (OEM + max payload)':<30}{r['MZFM']:>10,.0f}{100 * r['MZFM'] / mtom:>10.2f}")
    if p["furnishings_method"] == "torenbeek":
        print(f"  (Raymer 15.41 furnishings would give {r['raymer_furnishings']:,.0f} kg)")
    print(f"  (wing: Raymer 15.25 {r['wing_raymer']:,.0f} kg, Torenbeek 8-12 "
          f"{r['wing_torenbeek']:,.0f} kg; '{p['wing_method']}' used in the sum)")
    print(f"{'OEM (Class I)':<30}{p['OEM_class1']:>10,.0f}{100 * p['OEM_class1'] / mtom:>10.2f}")
    diff = r["OEM"] - p["OEM_class1"]
    print(f"{'Difference II - I':<30}{diff:>10,.0f}{100 * diff / p['OEM_class1']:>9.1f}%")

    print("\nValidation against ADSEE-II slide 41 (twin-engine jet transports, % MTOW)")
    print(f"{'Group':<18}{'Class II %':>11}{'min':>7}{'mean':>7}{'max':>7}  in range?")
    for g in SLIDE41_GROUPS:
        m = r["components"]["Flight controls"] if g == "surface_controls" else r["groups"][g]
        lo, mean, hi = slide41_range(g)
        pct = 100 * m / mtom
        print(f"{g:<18}{pct:>11.2f}{lo:>7.2f}{mean:>7.2f}{hi:>7.2f}  {'yes' if lo <= pct <= hi else 'NO'}")

    d = r["derived"]
    print("\nDerived inputs")
    for k, unit in [("S_flap", "m2"), ("S_ail", "m2"), ("S_csw", "m2"), ("S_cs", "m2"),
                    ("S_f", "m2"), ("L_struct", "m"), ("K_ws", "-"), ("F_w", "m"),
                    ("b_ht", "m"), ("l_nac", "m"), ("w_nac", "m"), ("S_n", "m2"),
                    ("V_pr", "m3"), ("L_ec", "m"), ("n_ult", "-"), ("N_l", "-")]:
        print(f"  {k:<10}{d[k]:>10.3f} {unit}")
    return r


if __name__ == "__main__":
    print_report()
