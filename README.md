# Chris Waudby's pulse program library

* These sequences have been uploaded for my own convenience in synchronising work across multiple spectrometers. Some experiments remain under development. Use at your own risk!
* Please acknowledge the use of sequences from this library.
* With thanks to John Kirkpatrick on whose work many sequences are based. 
* Contact: c.waudby@ucl.ac.uk

## Table of Contents
* [git quickstart](#git-quickstart)
* [syntax highlighting in vim](#bruker-syntax-highlighting-in-vim)
* [pulse programs](#pulse-programs)


# git quickstart

## Creating a local clone of the pp repository

```
git clone https://github.com/chriswaudby/pp.git
# move `pp` directory to required location
# then within `pp`, configure git:
git config user.name "Chris Waudby"
git config user.email "c.waudby@ucl.ac.uk"
git config color.ui true
git config color.status auto
git config color.branch auto
git config core.editor vim
git config merge.tool vimdiff
git config --list
```

## Updating local repo

```
git pull
```

This will *fetch* and *merge* remote changes.


## Making changes
Based on http://rogerdudler.github.io/git-guide/

0. Check current repo status:

```
git status
```

1. Add new files to the index:

```
git add <filename>
git add *
```

2. Commit changes:

```
git commit -m "Commit message"
```

Now file is committed in local repo but not in remote.

3. Push changes to remote:

```
git push origin master
```

Here `origin` refers to the remote repo (from which we originally cloned the repository) and
`master` is the branch to which changes should be pushed.





# Bruker syntax highlighting in vim

* Place `syntax` folder in `~/.vim`
* In vim: `set syntax=bruker`





# Pulse programs

Checkmarks indicate a tested sequence.

## Calibration

- [x] `zg490.cw` - 1H 90 calibration
- [x] `zgfd.cw` - 1H water flip-down
- [x] `zgfb.cw` - 1H water flip-back
- [x] `calib_15n.cw` - 15N high-power 90 calibration
- [x] `calib_13c.cw` - 13C high-power 90 calibration
- [x] `zg2h.cw` - 2H 1D/calibration

## Basic experiments

- [x] `zgesgp.cw` - 1H 1D (excitation sculpting)
- [ ] `zgpr.cw` - 1H 1D (presaturation)

## 1H,15N correlation experiments

- [x] `hsqcetfpf3gpsi2.cw` - 1H,15N sensitivity-improved HSQC
- [x] `fhsqcf3gpph.cw` - 1H,15N FHSQC
- [ ] `sfhmqcf3gpph.cw` - 1H,15N SOFAST-HMQC
- [ ] `sfhmqcf3gpphpr1d.cw` - 15N filtered/edited 1D SOFAST with presaturation
- [ ] `sfhmqcf3gpph.nuws.cw` - 1H,15N SOFAST-HMQC with NUWS
- [x] `hetsofast.cw` - based on IBS-HETSOFAST
- [ ] `hzdqcf3.cw`

## 1H,13C correlation experiments

- [ ] `hmqcgpphpr.cw` - 1H,13C HMQC
- [ ] `hmqcgpphpr_2.cw` - 1H,13C HMQC with purge element
- [ ] `hmqcgpphpr_2_wet.cw` - with WET solvent suppression
- [x] `hsqcctetgpsp.cw` - 1H,13C constant-time HSQC
- [x] `hisqcctetgpsp.cw` 1H,13C constant-time refocused (in-phase) HSQC
- [ ] `sfhmqcf2gpph.nuws.cw` - 1H,13C SOFAST-HMQC with NUWS
- [x] `sfhmqcf2et.nuws.cw` - 1H,13C gradient-selected SOFAST-HMQC

## Diffusion

- [x] `stebpgp1s19.cw` - 1H STE (watergate)
- [ ] `stebpgp1spr.cw` - 1H STE (presaturation)
- [ ] `stebpgp1s19xn.4.cw` - 15N XSTE
- [x] `sordid_pureZ.2pt.cw` - 15N SORDID


## Relaxation

### CPMG relaxation dispersion

- [ ] `hsqcrexfpf3gpsi_TS32.cw` - 15N SQ CW-CPMG
- [x] `ch3_mq_cpmg.cw` - 13C MQ CPMG for 13CH3 methyls (thanks to Vitali Tugarinov)
- [ ] `CHD2_1H_CPMG.cw` - 1H SQ CPMG for 13CHD2 methyls
- [ ] `C13_CHD2_METmethylCPMG.cw` - 13C SQ CPMG for 13CHD2 methyls (thanks to Goran Carlstrom, Lund)

### CEST

- [x] `cest_15N.cw` - 15N CEST
- [x] `cest_15N_B1cal.cw` - B1 calibration for 15N CEST
- [x] `13Cm_CEST.cw` - 13C methyl CEST
- [x] `13Cm_CEST_B1cal.cw` - B1 calibration for 13C methyl CEST
- [ ] `cest_zz_H.cw` - 1HN CEST (zz saturation)
- [ ] `cest_zz_H_B1cal.cw`
- [ ] `cest_zz_H_B1cal.1d.cw`
- [x] `cest_1H_zzstart.cw` - 1HN CEST
- [x] `cest_13Cali.cw` - 13Cali CEST
- [x] `cest_1HA.cw` - 1HA CEST
- [x] `cest_caconh.cw` - 13CA CEST
- [ ] `sfhmqcf3gpph_1H_cest.cw` - crude 1H CEST

### CCR

- [ ] `tract.cw` - TRACT
- [ ] `fab_nnh_tr_noSE.cw` - N/HN CCR (symmetric reconversion method)
- [x] `hmqcgpphpr_HHHH.cw` - HH/HH CCR (2Q or 3Q) for methyl order parameters
- [ ] `hsqcphpr_1Hcoupled.cw` - 1H,13C HSQC (F1-coupled)
- [ ] `adaptive_hsqcphpr.2d.cw`
- [ ] `adaptive_hsqcphpr.cw`
- [x] `best_trosy_ccr.cw` - formerly v.7


## Misc

- [x] `zgesgp.kinetics.cw` - 1H 1D for real-time measurements
- [ ] `sfhmqcf3gpph.3d.cw` - 1H,15N SOFAST-HMQC for real-time measurements



