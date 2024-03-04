

# Description of the Files 

## 1. Main File 
- ```simulation1.R``` is a (temporary) master file to run a simulation. It first generates data using `simDGP1.R` and then estimates the ATTs using `DDsurv.R`.

## 2. DGP
- log-linear DGP + (continuous version) estimation: ```log_linear_dgp.R```
- log-linear (OR Weibull) DGP: ```simDGP1.R```
- Weibull DGP: ```simDGP2.R```

## 3. DiD estimators
- ```DDsurv.R``` contains a function `DiD_supersurvivor` which first estimates the probability of being a supersurvivor using `supersurvivor.R` and then estimates the ATT using regression imputation. The function also calcuates ATT_reweight which is a reweight version of ATT using the estimated probability of being a supersurvivor as a weight.

-```InfeasibleDD.R``` contains a function `DD_infeasible` which estimates the ATT using the true supersurvivor as a control group, which is not available to empirical researcher.

## 3. Cured Model Estimation
- ```supersurvivor.R```: cured model estimation (allow all "Weibull", "LogNormal", "LogNormal_discrete")

## 4. Ignore the following files
- Callaway & Sant'anna: ```cs_estimator.R```
- ```InfeasibleDD.R```
- ```did_example.ipynb```
- ```LATE.R``` 
- ```dddddd.ipynb```