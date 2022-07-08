*-------------------------------------------------------------------------------
**# Instalation of ado file
net install github, from("https://haghish.github.io/github/")
github search qresid
github install psotob91/qresid

*-------------------------------------------------------------------------------

**# Normal example simulation
**## Simulation 


**## command: regress

**## command: glm

*-------------------------------------------------------------------------------
**# Bernoulli example simulation:
clear all
set seed 123
set obs 100
gen y = rbinomial(1, 0.3)
gen x = _n

**## command: glm
glm y x, family(binomial) link(logit)
predict dev, dev

qresid res1
qresid res2, nqres(1)
qresid res3, standardized nqres(1)

qnorm dev, name(qq1, replace) nodraw
qnorm res11, name(qq2, replace) nodraw
qnorm res21, name(qq3, replace) nodraw
qnorm res3_std1, name(qq4, replace) nodraw
graph combine qq1 qq2 qq3 qq4

**## command: logit
logit y x

capture drop res*
qresid res1
qresid res2, nqres(1)
qresid res3, standardized nqres(1)

qnorm dev, name(qq1, replace) nodraw
qnorm res11, name(qq2, replace) nodraw
qnorm res21, name(qq3, replace) nodraw
qnorm res3_std1, name(qq4, replace) nodraw
graph combine qq1 qq2 qq3 qq4


*-------------------------------------------------------------------------------
**# Binomial example simulation:
clear all
set seed 123
set obs 30
*gen n = runiformint(100, 200)
gen n = 20
gen y = rbinomial(n, 0.05)
gen x = _n

**## command: glm
glm y x, family(binomial n) link(logit)
predict dev, dev

qresid res1
qresid res2, nqres(1)
qresid res3, standardized nqres(1)

qnorm dev, name(qq1, replace) nodraw
qnorm res11, name(qq2, replace) nodraw
qnorm res21, name(qq3, replace) nodraw
qnorm res3_std1, name(qq4, replace) nodraw
graph combine qq1 qq2 qq3 qq4

*-------------------------------------------------------------------------------
**# Poisson example simulation:
**## Simulation 
clear all 
set seed 123
set obs 50
gen y = rpoisson(1)
gen x = _n

**## command: glm
glm y x, family(poisson)
predict dev, dev

qresid res1
qresid res2, nqres(1)
qresid res3, standardized nqres(1)

qnorm dev, name(qq1, replace) nodraw
qnorm res11, name(qq2, replace) nodraw
qnorm res21, name(qq3, replace) nodraw
qnorm res3_std1, name(qq4, replace) nodraw
graph combine qq1 qq2 qq3 qq4

**## command: poisson
poisson y x

capture drop res*
qresid res1
qresid res2, nqres(1)
qresid res3, standardized nqres(1)

qnorm dev, name(qq1, replace) nodraw
qnorm res11, name(qq2, replace) nodraw
qnorm res21, name(qq3, replace) nodraw
qnorm res3_std1, name(qq4, replace) nodraw
graph combine qq1 qq2 qq3 qq4


*-------------------------------------------------------------------------------
**# Gamma example simulation:



*------------------------------------------------------------------------------
**# Negative binomial simulation: 