"""
AE2111-I  WP3.3b  -  Calibration of the Raymer Class II equations (Raymer Sec. 15.4)
Group A06

Raymer Sec. 15.4: "fudge factor" = actual group weight / weight calculated by the
equations, determined on an existing aircraft and then applied to the new design.
This keeps the method pure Raymer (the equations are not mixed with another method);
only a multiplicative correction per group is added.

Actual group weights : ADSEE-II weight slides, slide 41 (Torenbeek, Table 8-?):
                       Boeing 737-200 (MTOW 100,000 lb) and Airbus A300B2 (MTOW 304,000 lb)
Geometry             : Jenkinson, Simpkin & Rhodes, "Civil Jet Aircraft Design",
                       Aircraft Data File (Boeing 737-200; Airbus A300-600R used for the
                       A300B2 wing, fuselage and tails, which are the same size)
                       https://booksite.elsevier.com/9780340741528/appendices/data-a/
Engines              : JT8D-9 and CF6-50C dry weights, public type data (verify!)
Everything else      : [ASSUMED] with the same rules as used for our own aircraft,
                       so that an assumption that biases the reference also biases ours.

The groups calibrated are the six groups that slide 41 gives:
wing, tail, fuselage, landing gear, surface controls, nacelle.
Propulsion and the other systems are not in slide 41 -> stay pure Raymer (factor 1).
"""

import class2_weight as c2

LB = 1 / c2.KG2LB
IN = 1 / c2.M2IN

# ------------------------------------------------------------------ slide 41: actual group weights [lb]
SLIDE41_ACTUAL = {
    "B737-200": dict(MTOW=100_000, wing=10_613, tail=2_718, fuselage=12_108,
                     landing_gear=4_354, surface_controls=2_348, nacelle=1_392),
    "A300B2":   dict(MTOW=304_000, wing=44_131, tail=5_941, fuselage=35_820,
                     landing_gear=13_611, surface_controls=5_808, nacelle=7_039),
}

# ------------------------------------------------------------------ reference geometry (SI)
REF = {
    "B737-200": dict(
        MTOM=100_000 * LB,              # slide 41 (weights belong to this MTOW)
        f_MLM=0.89,                     # MLW/MTOW 46,720/52,390            [Jenkinson]
        S_w=91.04, b_w=28.35, AR=8.83, taper=0.266, sweep_c4_deg=25.0,  #   [Jenkinson]
        tc_root=0.1289,                 # Jenkinson average t/c (same definition as our 0.12)
        flap_eta=(0.10, 0.76), flap_cf_c=0.30,   # flap span ratio 0.74     [Jenkinson / ASSUMED chord]
        ail_eta=(0.76, 0.96), ail_ca_c=0.25,     #                          [ASSUMED]
        l_fus=29.54, d_fus=3.73, d_fus_inner=3.53, l_cabin=20.6,           # [Jenkinson / ASSUMED cabin]
        nose_slender=1.6, tail_slender=2.8,      #                          [ASSUMED]
        S_ht=31.31, AR_ht=5.15, sweep_ht_deg=30.0,                         # [Jenkinson]
        S_vt=19.70, AR_vt=1.74, sweep_vt_deg=35.0,                         # [Jenkinson]
        L_t=14.78, x_ac_tail=27.8,               # tail arm [Jenkinson], x [ASSUMED]
        L_m=1.60, L_n=1.30,                      # strut lengths            [ASSUMED]
        N_mw=4, N_mss=2, N_nw=2, V_stall=52.0,   # wheels [Jenkinson], V_s  [ASSUMED]
        N_en=2, m_en=3_300 * LB, T_max_en=14_500 * 4.44822,                # JT8D-9 [public data]
        l_en=123.5 * IN, d_en=43.0 * IN,
        l_nac=7.00, d_nac=1.50, d_nac_exit=1.10, x_engine_front=11.0,     # [Jenkinson / ASSUMED]
        V_fuel=19.5, N_pax=130, N_flightcrew=2, N_cabincrew=3, m_cargo_max=2_600.0,
        payload_max=15_500.0,
    ),
    "A300B2": dict(
        MTOM=304_000 * LB,
        f_MLM=0.92,                     # MLW/MTOW A300B2                   [ASSUMED]
        S_w=260.0, b_w=44.84, AR=7.73, taper=0.30, sweep_c4_deg=28.0,     # [Jenkinson -600R]
        tc_root=0.105,                  # Jenkinson average t/c
        flap_eta=(0.11, 0.80), flap_cf_c=0.30,
        ail_eta=(0.82, 0.98), ail_ca_c=0.25,
        l_fus=53.30, d_fus=5.64, d_fus_inner=5.28, l_cabin=39.5,
        nose_slender=1.5, tail_slender=2.6,
        S_ht=69.45, AR_ht=3.81, sweep_ht_deg=34.0,
        S_vt=45.2, AR_vt=1.64, sweep_vt_deg=40.0,
        L_t=25.6, x_ac_tail=48.6,
        L_m=2.50, L_n=2.00,
        N_mw=8, N_mss=2, N_nw=2, V_stall=53.5,
        N_en=2, m_en=8_731 * LB, T_max_en=51_000 * 4.44822,               # CF6-50C [public data]
        l_en=183.0 * IN, d_en=105.0 * IN,
        l_nac=6.70, d_nac=2.70, d_nac_exit=1.90, x_engine_front=17.0,
        V_fuel=43.9, N_pax=251, N_flightcrew=3, N_cabincrew=6, m_cargo_max=9_000.0,
        payload_max=31_000.0,
    ),
}

CAL_GROUPS = ["wing", "tail", "fuselage", "landing_gear", "surface_controls", "nacelle"]


def method_groups(r):
    g = dict(r["groups"])
    g["surface_controls"] = r["components"]["Flight controls"]
    return g


def calibrate(verbose=True):
    factors = {}
    for name, geo in REF.items():
        p = {**c2.INPUTS, **geo, "fudge": {}}
        r = c2.class2(p)
        calc = method_groups(r)
        act = SLIDE41_ACTUAL[name]
        factors[name] = {g: act[g] * LB / calc[g] for g in CAL_GROUPS}
        if verbose:
            print(f"\n{name}  (MTOW {act['MTOW']:,} lb = {p['MTOM']:,.0f} kg)")
            print(f"  {'group':<17}{'Raymer [kg]':>12}{'actual [kg]':>12}{'factor':>8}")
            for g in CAL_GROUPS:
                print(f"  {g:<17}{calc[g]:>12,.0f}{act[g] * LB:>12,.0f}{factors[name][g]:>8.3f}")
    mean = {g: sum(f[g] for f in factors.values()) / len(factors) for g in CAL_GROUPS}
    if verbose:
        print("\nMean fudge factors (applied to our design):")
        for g in CAL_GROUPS:
            spread = [f[g] for f in factors.values()]
            print(f"  {g:<17}{mean[g]:6.3f}   (737-200 {spread[0]:.3f}, A300B2 {spread[1]:.3f})")
    return mean, factors


if __name__ == "__main__":
    mean, _ = calibrate()
    print("\n=== Our aircraft with the calibrated Raymer equations ===")
    c2.print_report({**c2.INPUTS, "fudge": mean})
