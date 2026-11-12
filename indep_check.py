"""Independent check: predict the OEM of two aircraft NOT used in the calibration (Jenkinson data A)."""
import json, class2_weight as c2
LB, IN = 1 / c2.KG2LB, 1 / c2.M2IN
d = json.load(open("calib3.json"))
CHECK = {
 "A321-200": (48_000, dict(
    MTOM=89_000, f_MLM=73_500/89_000, S_w=122.4, b_w=33.91, AR=9.39, taper=0.240, sweep_c4_deg=25.0, tc_root=0.12,
    flap_eta=(0.10, 0.78+0.10), flap_cf_c=0.25, ail_eta=(0.80, 0.97), ail_ca_c=0.25,
    l_fus=44.51, d_fus=3.95, d_fus_inner=3.70, l_cabin=34.4, nose_slender=1.6, tail_slender=2.8,
    S_ht=31.0, AR_ht=5.0, sweep_ht_deg=29.0, S_vt=21.5, AR_vt=1.82, sweep_vt_deg=34.0,
    L_t=16.2, x_ac_tail=41.0, Ht_Hv=0.0, F_w=None,
    L_m=2.0, L_n=1.8, N_mw=4, N_mss=2, N_nw=2, V_stall=138*0.5144/1.23,
    N_en=2, N_en_wing=2, m_en=5_250*LB, T_max_en=142e3, l_en=102.4*IN, d_en=68.3*IN,
    l_nac=4.44, d_nac=2.37, d_nac_exit=1.60, x_engine_front=14.0,
    V_fuel=23.7, N_pax=186, N_cabincrew=4, m_cargo_max=3_000.0, m_uld=0.0, payload_max=22_780.0,
    long_range=False, overwater=False)),
 "757-200 (PW2037)": (58_040, dict(
    MTOM=115_900, f_MLM=95_450/115_900, S_w=185.25, b_w=38.05, AR=7.82, taper=0.243, sweep_c4_deg=25.0, tc_root=0.12,
    flap_eta=(0.10, 0.757+0.10), flap_cf_c=0.25, ail_eta=(0.80, 0.97), ail_ca_c=0.25,
    l_fus=46.96, d_fus=4.00, d_fus_inner=3.54, l_cabin=36.1, nose_slender=1.6, tail_slender=2.8,
    S_ht=50.35, AR_ht=4.59, sweep_ht_deg=27.5, S_vt=34.37, AR_vt=1.56, sweep_vt_deg=39.0,
    L_t=19.85, x_ac_tail=43.5, Ht_Hv=0.0, F_w=None,
    L_m=2.4, L_n=2.0, N_mw=8, N_mss=2, N_nw=2, V_stall=137*0.5144/1.23,
    N_en=2, N_en_wing=2, V_fuel=42.6, N_pax=186, N_cabincrew=4, m_cargo_max=4_000.0, m_uld=0.0,
    payload_max=25_690.0, long_range=True, overwater=True)),
}
print(f"{'aircraft':<18}{'actual OEM':>11}{'k(2 twins)':>12}{'err':>7}{'k(3 twins)':>12}{'err':>7}{'pure Raymer':>12}{'err':>7}")
for n, (oem, geo) in CHECK.items():
    row = []
    for f in (d["mean2"], d["mean3"], {}):
        r = c2.class2({**c2.INPUTS, **geo, "fudge": f}); row += [r["OEM"], 100 * (r["OEM"] / oem - 1)]
    print(f"{n:<18}{oem:>11,}" + "".join(f"{row[i]:>12,.0f}{row[i+1]:>+6.1f}%" for i in (0, 2, 4)))
