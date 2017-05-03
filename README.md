# Chris Waudby's pulse program library

* Use at your own risk!
* Please acknowledge the use of sequences from this library.
* With thanks to John Kirkpatrick on whose work many sequences are based. 
* Contact: c.waudby@ucl.ac.uk


## Calibration

- [x] `zg490.cw` - 1H 90 calibration
- [x] `zgfd.cw` - 1H water flip-down
- [x] `zgfb.cw` - 1H water flip-back
- [ ] `calib_15n.cw` - 15N high-power 90 calibration
- [ ] `calib_13c.cw` - 13C high-power 90 calibration
- [x] `zg2h.cw` - 2H 1D/calibration

## Basic experiments

- [ ] `zgesgp.cw` - 1H 1D (excitation sculpting)
- [ ] `zgpr.cw` - 1H 1D (presaturation)
- [ ] `hmqcgpphpr.cw` - 13C HMQC

## 1H,15N correlation experiments

- [ ] `sfhmqcf3gpphpr1d.cw` - 1D 15N SOFAST with presaturation
- [ ] `sfhmqcf3gpph.nuws.cw` - 1H,15N SOFAST-HMQC with NUWS

## 1H,13C correlation experiments

- [ ] `sfhmqcf2gpph.nuws.cw` - 1H,13C SOFAST-HMQC with NUWS
- [x] `sfhmqcf2et.nuws.cw` - 1H,13C gradient-selected SOFAST-HMQC

## Diffusion

- [x] `stebpgp1s19.cw` - 1H STE (watergate)
- [ ] `stebpgp1spr.cw` - 1H STE (presaturation)
- [x] `sordid_pureZ.2pt.cw` - 15N SORDID


## Relaxation

### CPMG relaxation dispersion

- [ ] `hsqcrexfpf3gpsi_TS32.cw` - 15N SQ CW-CPMG
- [ ] `ch3_mq_cpmg.cw` - 13C MQ CPMG for 13CH3 methyls (thanks to Vitali Tugarinov)
- [ ] `CHD2_1H_CPMG.cw` - 1H SQ CPMG for 13CHD2 methyls
- [ ] `C13_CHD2_METmethylCPMG.cw` - 13C SQ CPMG for 13CHD2 methyls (thanks to Goran Carlstrom, Lund)

### CEST

- [ ] `cestN15.cw` - 15N CEST
- [x] `13Cm_CEST.cw` - 13C methyl CEST
- [x] `13Cm_CEST_B1cal.cw` - B1 calibration for 13C methyl CEST

### CCR

- [ ] `hmqcgpphpr_HHdipCC.cw`


## Misc

- [x] `zgesgp.kinetics.cw` - 1H 1D for real-time measurements

