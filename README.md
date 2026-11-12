# Group-A06-System-Design-AE2111
Repo for Group A06 from the System Design Course of the Aerospace Engineering Bachelor at TU Delft. 2026-2027 academic year! 
AE2111-I WP3 - Group A06 - Class II drag and mass estimation and Class I/II iteration
Python 3, standard library only. Run the scripts from this folder.

Main result (Chapter 4 and 5)
  python iteration.py        Iteration history and converged design
                             (Table 5.1, Table 5.2, 5.1.3; MTOM 95 103 kg).
                             Recalculates the calibration factors first and warns
                             if the factors in class2_weight.py are out of date.

Modules used by the iteration (each can also be run on its own; when run on its own,
the input geometry is the WP1/WP2 design before the iteration, not the converged design)
  class1_weight.py           Class I mass and fuel fraction, eta_eng = V/(TSFC e_f)
  class2_drag.py             Class II zero-lift drag build-up, Oswald factor, wave drag (Table 4.1)
  class2_weight.py           Class II mass, Raymer 15.25-15.43 + Torenbeek 8-44 / Table 8-13
                             (Tables 4.3, 4.8, 4.9). All inputs with their sources are in INPUTS.
  calibrate_raymer.py        Calibration factors from the B737-200 and A300B2 (Table 4.4)

Checks of the calibration (Section 4.2.2 and 4.2.5)
  python calib3.py           Third reference aircraft, DC-9/10 (Table 4.5); writes calib3.json
  python indep_check.py      Validation on the A321-200 and B757-200 (Table 4.10); needs calib3.json
