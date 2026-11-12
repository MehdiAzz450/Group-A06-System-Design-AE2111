"""Third reference twin (DC-9-10), leave-one-out check and independent OEM checks."""
import itertools, class2_weight as c2, calibrate_raymer as cr
LB, IN = cr.LB, cr.IN
# slide 41 % MTOW of the DC-9/10 (MTOW 77,700 lb) -> actual lb
# slide 41 absolute group weights [lb] of the DC-9/10 (MTOW 91,500 lb)
cr.SLIDE41_ACTUAL["DC-9-10"] = dict(MTOW=91_500, wing=9_470, tail=2_630, fuselage=11_206,
                                    landing_gear=3_660, surface_controls=1_264, nacelle=1_417)
cr.REF["DC-9-10"] = dict(
    MTOM=91_500 * LB, f_MLM=33_565 / (91_500 * LB),                  # slide 41 MTOW; MLM Jenkinson table 5
    S_w=86.77, b_w=27.25, AR=8.56, taper=0.194, sweep_c4_deg=24.5,   # Jenkinson; sweep: Douglas data [verify]
    tc_root=0.1165,
    flap_eta=(0.11, 0.78), flap_cf_c=0.30, ail_eta=(0.78, 0.96), ail_ca_c=0.25,   # flap span 0.67 Jenkinson
    l_fus=27.91, d_fus=3.60, d_fus_inner=3.10, l_cabin=17.5,
    nose_slender=1.6, tail_slender=2.8,
    S_ht=25.60, AR_ht=4.93, sweep_ht_deg=30.0, S_vt=14.96, AR_vt=1.02, sweep_vt_deg=42.0,
    L_t=14.35, x_ac_tail=26.5, Ht_Hv=1.0,                            # T-tail
    L_m=1.40, L_n=1.20, N_mw=4, N_mss=2, N_nw=2, V_stall=115 * 0.5144 / 1.3,
    N_en=2, N_en_wing=0, m_en=3_100 * LB, T_max_en=54.5e3,          # JT8D-5
    l_en=120.0 * IN, d_en=43.0 * IN,
    l_nac=5.30, d_nac=1.56, d_nac_exit=1.15, x_engine_front=20.5,
    V_fuel=14.0, N_pax=90, N_flightcrew=2, N_cabincrew=2, m_cargo_max=1_700.0, payload_max=9_570.0,
)
G = cr.CAL_GROUPS
mean, fac = cr.calibrate(verbose=False)
print(f"{'group':<17}" + "".join(f"{n:>10}" for n in fac) + f"{'mean3':>8}{'mean2':>8}{'CoV3':>7}")
import statistics as st
m2 = {g: (fac['B737-200'][g] + fac['A300B2'][g]) / 2 for g in G}
for g in G:
    v = [fac[n][g] for n in fac]
    print(f"{g:<17}" + "".join(f"{x:>10.3f}" for x in v) + f"{mean[g]:>8.3f}{m2[g]:>8.3f}{st.stdev(v)/st.mean(v):>7.2f}")
# leave-one-out: predict each reference with the mean of the other two
print("\nLeave-one-out: group mass error when the aircraft is NOT in the calibration set")
for out in fac:
    k = {g: st.mean([fac[n][g] for n in fac if n != out]) for g in G}
    r = c2.class2({**c2.INPUTS, **cr.REF[out], "fudge": {}})
    calc = cr.method_groups(r); act = cr.SLIDE41_ACTUAL[out]
    tot_p = sum(calc[g] * k[g] for g in G); tot_a = sum(act[g] * LB for g in G)
    print(out, " ".join(f"{g[:5]}{100*(calc[g]*k[g]/(act[g]*LB)-1):+6.1f}%" for g in G), f" six-group sum {100*(tot_p/tot_a-1):+.1f}%")
import json; json.dump(dict(mean3=mean, mean2=m2, fac=fac), open("calib3.json", "w"), indent=1)
