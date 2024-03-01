

# Description of the Files 

## 1. Main File 
- ```DDsurv.R```: uses ```supersurvivor.R```. Seems to do both infeasible DD and our estimator? 
- ```simulation1.R```: uses ```simDGP1.R```, ```DDsurv.R```, ```InfeasibleDD.R``` (```DDsurv.R```uses ```supersurvivor.R```) Seems to do both infeasible DD and our estimator? and reweighted ols??? (seems to have some overlaps with ```DDsurv.R```)
- ```simulation2.R```: uses ```simDGP2.R```, ```supersurvivor.R```, ```DDsurv.R```, ```InfeasibleDD.R```. 


## 2. DGP
- log-linear DGP + (continuous version) estimation: ```log_linear_dgp.R```
- log-linear (OR Weibull) DGP: ```simDGP1.R```
- Weibull DGP: ```simDGP2.R```

## 3. Cured Model Estimation
- ```supersurvivor.R```: cured model estimation (allow all "Weibull", "LogNormal", "LogNormal_discrete")

## 3. ATT Estimation
- Callaway & Sant'anna: ```cs_estimator.R```
- Infeasible DD (When we can correctly identify super-survivor): ```InfeasibleDD.R```
- Our Estimator: ???? 

## Don't know what it is (PLEASE ADD DESCRIPTIONS)
- ```did_example.ipynb```
- ```LATE.R``` 
- ```dddddd.ipynb```


--- 

# Which files are we using in the end? 
- ```DDsurv.R```
- ```InfeasibleDD.R```
- ```simDGP1.R```
- ```simulation1.R```
- ```supersurvivor.R```

Please add the difference between ```DDsurv.R``` and ```simulation1.R```. Seems to have some overlaps. 
```DDsurv.R``` shows error message like: object 'DDsurv' not found. Seems to have deleted the object recently. 
