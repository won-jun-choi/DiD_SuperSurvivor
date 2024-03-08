

# Description of the Files 

## 1. Main File
- `SampleSimulationCode.R`: a sample code to run a simulation. Please modify DGP in this file.
- ```simulation1.R``` is a (temporary) master file to run a simulation. It first generates data using `simDGP1.R` and then estimates the ATTs using `DDsurv.R`.

## 2. DGP
- log-linear DGP + (continuous version) estimation: ```log_linear_dgp.R```
- log-linear (OR Weibull) DGP: ```simDGP1.R```
- Weibull DGP: ```simDGP2.R```

## 3. DiD estimators
- `DiD_supersurvivor.R`
- `DiD_supersurvivor_noControl.R`: same as the above, but no control variables in Y regression.
- `DiD_CS.R`: ATT_RI of Callaway & Sant'anna (2020).
- `DiD_infeasible.R`: estimates the ATT using the true supersurvivor as a control group

## 3. Cured Model Estimation
- ```supersurvivor.R```: cured model estimation (allow all "Weibull", "LogNormal", "LogNormal_discrete")
