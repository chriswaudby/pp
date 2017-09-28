# Chris Waudby's pulse program library

* These sequences have been uploaded for my own convenience in synchronising work across multiple spectrometers. Some experiments remain under development. Use at your own risk!
* Please acknowledge the use of sequences from this library.
* With thanks to John Kirkpatrick on whose work many sequences are based. 
* Contact: c.waudby@ucl.ac.uk



# Table of Contents
* [git quickstart](#git-quickstart)
* [syntax highlighting in vim](#bruker-syntax-highlighting-in-vim)
* [my shell settings](#my-shell-settings)
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

* Copy contents of `syntax-highlighting` folding into `~/.vim`
* Pulse programs will be automatically detected by the presence of `#include <Avance.incl>`
* Alteratively, in vim: `set syntax=bruker`



# My shell settings

## .bashrc

```
# history searching
bind '"\e[A": history-search-backward'
bind '"\e[B": history-search-forward'

# append history (in case of multiple sessions)
shopt -s histappend

# nice prompt
export PS1="\[\033[38;5;15m\][\[$(tput sgr0)\]\[\033[38;5;11m\]\u@\h\[$(tput sgr0)\]\[\033[38;5;15m\] \[$(tput sgr0)\]\[\033[38;5;13m\]\W\[$(tput sgr0)\]\[\033[38;5;15m\]]\\$ \[$(tput sgr0)\]"
```

## .cshrc

```
set savehist = (5000 merge )
set filec
set nobeep
set autolist=true

bindkey -k up history-search-backward # Use up and down arrow to search
bindkey -k down history-search-forward

# enable colouring for ls
setenv CLICOLOR

# custom prompt
set prompt="%B%{\033[36m%}[%c2]%b % "

# add utils to path
set path =($path ~/pp/util)
```



# Pulse programs

Checkmarks indicate a tested sequence.

* [Calibration](#calibration)
* [1Ds](#1Ds)
* [1H experiments](#1H-experiments)
* [1H,15N correlation experiments](#1H,15N-correlation-experiments)
* [1H,13C correlation experiments](#1H,13C-correlation-experiments)
* [Assignment](#assignment)
* [Diffusion](#diffusion)
* [Relaxation and dynamics](#relaxation-and-dynamics)
* [Miscellaneous](#misc)


## Calibration

- [x] `zg490.cw` - 1H 90 calibration
- [x] `zgfd.cw` - 1H water flip-down
- [x] `zgfb.cw` - 1H water flip-back
- [x] `calib_15n.cw` - 15N high-power 90 calibration
- [x] `calib_13c.cw` - 13C high-power 90 calibration
- [ ] `calib_1h_nut.cw` - 1H calibration (2d array of pulse lengths)
- [x] `zg2h.cw` - 2H 1D/calibration


## 1Ds

- [x] `zgesgp.cw` - 1H 1D (excitation sculpting)
- [x] `zgpr.cw` - 1H 1D (presaturation)
- [x] `zg_doublesolventsuppression.cw` - 1H 1D (2x solvent flipdowns)
- [x] `cpmgpr1d` - 1H CPMG with presat (TS3.5pl5 library)
- [x] `zgpg30` - 13C 1D (two-level cpd during d1 and acquisition, TS3.5pl5 library)


## 1H experiments

- [x] `project.cw` - PROJECT (T2 relaxation via perfect echo)
- [ ] `stddiffesgp.cw` - STD (saturation transfer difference)
- [x] `wlogsy.cw` - WaterLOGSY


## 1H,15N correlation experiments

- [x] `hsqcetfpf3gpsi2.cw` - 1H,15N sensitivity-improved HSQC
- [ ] `hsqcfpf3gpphwg.cw` - 1H,15N HSQC (NB try this instead of fhsqc)
- [x] `fhsqcf3gpph.cw` - 1H,15N FHSQC
- [ ] `sfhmqcf3gpph.cw` - 1H,15N SOFAST-HMQC
- [ ] `trosyetf3gpsi.cw` - 1H,15N TROSY
- [ ] `sfhmqcf3gpphpr1d.cw` - 15N filtered/edited 1D SOFAST with presaturation
- [ ] `sfhmqcf3gpph.nuws.cw` - 1H,15N SOFAST-HMQC with NUWS
- [x] `hetsofast.cw` - based on IBS-HETSOFAST
- [x] `hzdqcf3.cw`
- [x] `sfhzdqcf3.cw`


## 1H,13C correlation experiments

- [x] `hmqcgpphpr.cw` - 1H,13C HMQC
- [x] `hmqcgpphpr_2.cw` - 1H,13C HMQC with purge element
- [ ] `hmqcgpphpr_2_wet.cw` - with WET solvent suppression
- [x] `hsqcctetgpsp.cw` - 1H,13C constant-time HSQC
- [x] `hisqcetgpsp.cw` - 1H,13C non-CT refocused (in-phase) HSQC
- [x] `hisqcctetgpsp.cw` 1H,13C constant-time refocused (in-phase) HSQC
- [ ] `sfhmqcf2gpph.nuws.cw` - 1H,13C SOFAST-HMQC with NUWS
- [x] `sfhmqcf2et.nuws.cw` - 1H,13C gradient-selected SOFAST-HMQC


## Diffusion

- [x] `stebpgp1s.cw` - 1H/19F STE (no solvent suppression)
- [x] `stebpgp1s19.cw` - 1H STE (watergate)
- [x] `stebpgp1spr.cw` - 1H STE (presaturation)
- [x] `stebpgp1spr_13c.cw` - 1H STE (presat) with 13C decoupling during acquisition
- [x] `stebpgp1s19xn.4.cw` - 15N XSTE
- [x] `sordid_pureZ.2pt.cw` - 15N SORDID
- [x] `stebpgp1sprxc_dec.3.cw` - 13C XSTE


## Assignment
 - [x] `b_hncogp3d.2` - BEST-HNCO (TS3.5pl6 library)
 - [x] `b_hncacogp3d` - BEST-HNCACO (TS3.5pl6 library)
 - [x] `b_hncocacbgp3d.2` - BEST-HNCOCACB (TS3.5pl6 library)
 - [x] `b_hncacbgp3d.2` - BEST-HNCACB (TS3.5pl6 library)
 - [x] `ccconhgp3d.2` - CCCONH TOCSY (TS3.5pl6 library)
 - [x] `hcchcogp3d` - HC(C)H COSY (TS3.5pl6 library)
 - [ ] `hbhaconhgp3d.cw` - HBHACONH
 - [x] `hccconhgpwg3d2.cw` - H(CCCO)NH TOCSY
 - [x] `hcchdigp3d.cw` - HC(C)H TOCSY
 - [x] `hmcmcgcbca.linear.2h.cw` - Cali assignment for linearised methyl backbones
 - [ ] `hmcmcgcbcaco.linear.2h.cw` - CO assignment for linearised methyl backbones


## Relaxation and dynamics

* [15N](#15n-relaxation)
* [methyl](#methyl)
* [CPMG relaxation dispersion](#cpmg-relaxation-dispersion)
* [CEST](#cest)
* [EXSY](#exsy)

### 15N relaxation

- [ ] `T1N15_tr.cw` - 15N T1 (TROSY)
- [ ] `T1rhoN15_tr.cw` - 15N T1rho (TROSY)
- [ ] `T1rhoN15_tr_calib.cw`
- [ ] `N15HetNoe_tr.cw` - 15N heteronuclear NOE (TROSY)
- [ ] `hsqcfpf3gpphwg_T2.cw` - amide 1H T2 measurement (for PREs)
- [x] `tract.cw` - TRACT
- [ ] `fab_nnh_tr_noSE.cw` - N/HN CCR (symmetric reconversion method)
- [x] `best_trosy_ccr.cw` - formerly v.7

### methyl 

- [x] `hmqcgpphpr_1HT2.cw` - methyl 1H relaxation (L2 line, Tugarinov 2006)
- [x] `hmqcgpphpr_1HT2.1.cw` - methyl 1H relaxation (all lines)
- [x] `hmqcgpphpr_1HT2.2.cw` - methyl 1H relaxation (L1+L3 lines, Tugarinov 2006)
- [x] `hmqcgpphpr_1HT2.3.cw` - methyl 1H relaxation (L1-L3 lines, Tugarinov 2006)
- [x] `hmqcgpphpr_HHHH.cw` - HH/HH CCR (2Q or 3Q) for methyl order parameters (Tugarinov 2007, Sun 2011)
- [ ] `hmqcgpphpr_HHHC.cw` - HH/HC CCR for methyl order parameters
- [x] `hsqcphpr_1Hcoupled.cw` - 1H,13C HSQC (F1-coupled)
- [x] `hsqcphpr_1Hcoupled.2.cw` - 1H,13C HSQC (F1-coupled) with multiplet selection
- [x] `hsqcphpr_1Hcoupled.3.cw` - 1H,13C HSQC (F1-coupled) with multiplet selection (pseudo-3D)
- [x] `13C_T1.cw` - 13C T1
- [x] `13C_NOE.cw` - 13C NOE
- [x] `hsqcphpr_1Hcoupled_T1p.cw` - 13C T1p (with multiplet editing)
- [ ] `adaptive_hsqcphpr.2d.cw`
- [ ] `adaptive_hsqcphpr.cw`

### CPMG relaxation dispersion

- [ ] `hsqcrexfpf3gpsi_TS32.cw` - 15N SQ CW-CPMG
- [x] `ch3_mq_cpmg.cw` - 13C MQ CPMG for 13CH3 methyls (thanks to Vitali Tugarinov)
- [ ] `CHD2_1H_CPMG.cw` - 1H SQ CPMG for 13CHD2 methyls
- [ ] `C13_CHD2_METmethylCPMG.cw` - 13C SQ CPMG for 13CHD2 methyls (thanks to Goran Carlstrom, Lund)
- [x] `TQ_CPMG.cw` - 1H 3Q CPMG for methyls

### CEST

- [x] `cest_15N.cw` - 15N CEST
- [x] `cest_15N_B1cal.cw` - B1 calibration for 15N CEST
- [x] `13Cm_CEST.cw` - 13C methyl CEST
- [x] `13Cm_CEST_B1cal.cw` - B1 calibration for 13C methyl CEST
- [ ] `cest_zz_H.cw` - 1HN CEST (zz saturation)
- [ ] `cest_zz_H_B1cal.cw`
- [ ] `cest_zz_H_B1cal.1d.cw`
- [x] `cest_1H_zzstart.cw` - 1HN CEST
- [x] `cest_13Cali.cw` - 13Cali CEST (constant time t1 evolution)
- [x] `cest_13Cali_nonCT.cw` - 13Cali CEST (non-constant time)
- [x] `cest_13Cali_B1cal.cw`
- [x] `cest_1HA.cw` - 1HA CEST
- [x] `cest_caconh.cw` - 13CA CEST
- [ ] `sfhmqcf3gpph_1H_cest.cw` - crude 1H CEST

### EXSY

- [ ] `hsqc_Nz_exsy.cw` - 15N EXSY


## Misc

- [x] `zgesgp.kinetics.cw` - 1H 1D for real-time measurements
- [ ] `sfhmqcf3gpph.3d.cw` - 1H,15N SOFAST-HMQC for real-time measurements
- [x] `19F_CEST.cw` - 19F CEST (pseudo-2D)
- [x] `19F_xQ.cw` - 19F DD/DD CCR measurements (pseudo-2D)
- [x] `19F_T2.cw` - 19F hahn-echo
- [x] `hsqcphpr_19Fcoupled.cw` - 13C-19F HSQC
- [ ] `imaging_zgesgp.cw` - 1D CSI

