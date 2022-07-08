*! version 0.0.5 Percy Soto-Becerra 06ago2021

/***
{* *! version 1.0  5 Jul 2022}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "Install command2" "ssc install command2"}{...}
{vieweralsosee "Help command2 (if installed)" "help command2"}{...}
{viewerjumpto "Syntax" "qresid##syntax"}{...}
{viewerjumpto "Description" "qresid##description"}{...}
{viewerjumpto "Options" "qresid##options"}{...}
{viewerjumpto "Remarks" "qresid##remarks"}{...}
{viewerjumpto "Examples" "qresid##examples"}{...}
{viewerjumpto "References" "qresid##references"}{...}

Title
====== 

{bf:qresid} {hline 2} (Randomized) quantile residuals for diagnostic of models 


Syntax
------ 

> __qresid__ _newvarname_ [_{help if}_] [_{help in}_] [, _options_]

Options
-------

| _options_         |  Description                                                             |
|:------------------|:-------------------------------------------------------------------------|
| standardized      | compute leverage standardized quantile residuals                        |
| nqres(#)             | number of randomized quantile residuals; # default is 4 for discrete outcomes; not apply in continuous outcomes                   |


Description
-----------

__qresid__ calculates (randomized) quantile residuals for the diagnostic of a great variety of models. 

This version of __qresid__ only allows to create quantile residuals for _generalized linear models_. Further updates will include _zero inflated models_, _beta regression_ and others _related models_. 


Remarks
--------

__qresid__ implements the computation of quantile residuals proposed by [Peter K. Dunn and Gordon K. Smyth (1996)](http://www.jstor.org/stable/1390802) for generalized linear models (GLM). 

The residuals and their plots help to evaluate the assumptions of the statistical models. In the general linear model with normal outcome distribution, the residuals (raw, standardized, and studentized) follow a normal distribution when the model is right. 

In the case of generalized linear models (GLMs) for non-normal outcome distributions, the deviance and Pearson residuals (and their standardized forms) do not necessarily follow a normal distribution. This makes model diagnostic difficult tin many practical scenarios. As if that were not enough, the GLM residuals for discrete outcomes can present patterns of parallel bands that artifact their proper visual inspection. 


Examples
--------

	__1) Outcome with binomial distribution__
	
	~~~
		. clear all
		
		. set seed 123
		
		. set obs 30
		
		. gen n = 20
		
		. gen y = rbinomial(n, 0.05)
		
		. gen x = _n

		. glm y x, family(binomial n) link(logit)
		
		. predict dev, dev

		. qresid res1
		
		. qresid res2, nqres(1)
		
		. qresid res3, standardized nqres(1)

		. qnorm dev, name(qq1, replace) nodraw
		
		. qnorm res11, name(qq2, replace) nodraw
		
		. qnorm res21, name(qq3, replace) nodraw
		
		. qnorm res3_std1, name(qq4, replace) nodraw
		
		. graph combine qq1 qq2 qq3 qq4
	~~~
		
	__2) Outcome with Poisson distribution__
	
	~~~
		. clear all 
		
		. set seed 123
		
		. set obs 50
		
		. gen y = rpoisson(1)
		
		. gen x = _n

		. glm y x, family(poisson)
		
		. predict dev, dev

		. qresid res1
		
		. qresid res2, nqres(1)
		
		. qresid res3, standardized nqres(1)

		. qnorm dev, name(qq1, replace) nodraw
		
		. qnorm res11, name(qq2, replace) nodraw
		
		. qnorm res21, name(qq3, replace) nodraw
		
		. qnorm res3_std1, name(qq4, replace) nodraw
		
		. graph combine qq1 qq2 qq3 qq4

		. poisson y x

		. capture drop res*
		
		. qresid res1
		
		. qresid res2, nqres(1)
		
		. qresid res3, standardized nqres(1)

		. qnorm dev, name(qq1, replace) nodraw
		
		. qnorm res11, name(qq2, replace) nodraw
		
		. qnorm res21, name(qq3, replace) nodraw
		
		. qnorm res3_std1, name(qq4, replace) nodraw
		
		. graph combine qq1 qq2 qq3 qq4
		~~~
		
		
Author
------

Percy Soto-Becerra  
InkaStats Data Science Solutions - Medical Branch   
percys1991@gmail.com


References
----------
{pstd}

Peter K. Dunn & Gordon K. Smyth. 1996. 
[Randomized Quantile Residuals](http://www.jstor.org/stable/1390802). Journal of Computational and Graphical Statistics, 5:3, 236-244.

- - -

This help file was dynamically produced by 
[MarkDoc Literate Programming package](http://www.haghish.com/markdoc/) 
***/

program define qresid
	version 15.0
	syntax newvarname(max=1) [if] [in] [, standardized nqres(int 4)] 
	
	/*Quantile residual for normal linear regression (by Ordinary Least Square 
	or Maximum Likelihood estimation*/
	
	if "`e(cmd)'" == "regress" {
		tempvar residuals
		predict double `residuals', residuals
		gen `typlist' `varlist' = `residuals' / `e(rmse)' `if' `in'
		label variable `varlist' "Quantile residuals"
		
		/*Creating standardized quantile residual*/
		
		if "`standardized'" == "standardized" {
			tempvar hat
			local stand  std
			predict double `hat', hat
			gen `typlist' `varlist'_`stand' = `varlist' / sqrt(1 - `hat') 
			label variable `varlist'_`stand' "Standardized quantile residuals"
		}
	
	}
	else if "`e(cmd)'" == "glm" & "`e(varfunct)'" == "Gaussian" {
		tempvar deviance_res
		predict double `deviance_res', deviance
		gen `typlist' `varlist' = `deviance_res' / `e(dispers)' `if' `in'
		label variable `varlist' "Quantile residuals"
		
		/*Creating standardized quantile residual*/
		
		if "`standardized'" == "standardized" {
			tempvar hat
			local stand std
			predict double `hat', hat
			gen `typlist' `varlist'_`stand' = `varlist' / sqrt(1 - `hat') `if' `in'
			label variable `varlist'_`stand' "Standardized quantile residuals"
		}
	}
	
	/*Quantile residual for bernoulli/binomial logistic regression 
	(GLM, family = Bernoulli/Binomial, link = any admissible)*/
	
	else if "`e(cmd)'" == "logit" | "`e(cmd)'" == "logistic" {
		tempvar n p y a b u
		predict double `p', pr
		gen double `n' = 1 /*Prior weights*/
		gen double `y' = `n' * `p'
		gen double `a' = binomial(`n', `y' - 1, `p')
		gen double `b' = binomial(`n', `y', `p')
		foreach i of numlist 1 / `nqres' {
			gen double `u'`i' = runiform(`a', `b')
			gen `typlist' `varlist'`i' = invnormal(`u'`i') `if' `in'
			label variable `varlist'`i' "Randomized quantile residuals `i'"
		}
		
		/*Creating standardized quantile residual*/
		
		if "`standardized'" == "standardized" {
			local oldformula = "`e(cmdline)'"
			if "`e(cmd)'" == "logit" {
			    local preformula = substr("`oldformula'", 6, .)
			} 
			else if "`e(cmd)'" == "logistic" {
			    local preformula = substr("`oldformula'", 9, .)
			}
			local pos = strpos("`preformula'", ",") - 1
			if "`pos'" == "-1" {
				local new_formula = "`preformula'"
			} 
			else {
				local new_formula = substr("`preformula'", 1, `pos')
			}

			quietly glm `new_formula', family(binomial) link(logit)
			tempvar hat
			local stand std
			predict double `hat', hat
			
			quietly `oldformula'
			
			foreach i of numlist 1 / `nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
				label variable `varlist'_`stand'`i' "Stand. rand. quantile residuals `i'"
			}

		}
	} 
	else if "`e(cmd)'" == "glm" & ("`e(varfunct)'" == "Bernoulli" | /// 
		"`e(varfunct)'" == "Binomial") {
		tempvar y n mu p a b u
		gen double `n' = `e(m)' /*Prior weights*/
		predict double `mu', mu
		gen double `p' = `mu' / `n'
		gen double `y' = `n' * `p'
		gen double `a' = binomial(`n', `y' - 1, `p')
		gen double `b' = binomial(`n', `y', `p')
		foreach i of numlist 1 / `nqres' {
			gen double `u'`i' = runiform(`a', `b')
			gen `typlist' `varlist'`i' = invnormal(`u'`i') `if' `in'
			label variable `varlist'`i' "Randomized quantile residuals `i'"
		}
		
		/*Creating standardized quantile residual*/
		
		if "`standardized'" == "standardized" {
			tempvar hat
			local stand std
			predict double `hat', hat
			foreach i of numlist 1 / `nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
				label variable `varlist'_`stand'`i' "Stand. rand. quantile residuals `i'"
			}
		}
	} 
	
	/*Quantile residual for Poisson regression 
	(GLM, family = Poisson, link = any admissible)*/
	
	else if "`e(cmd)'" == "poisson" {
		tempvar y n mu p a b u
		gen double `y' = `e(depvar)'
		predict double `mu', n /*Fitted values in count data n = exp(ir + exposure)*/
		gen double `a' = cond(`y' > 0, poisson(`mu', `y' - 1), 0)
		gen double `b' = poisson(`mu', `y')
		foreach i of numlist 1 / `nqres' {
			gen double `u'`i' = runiform(`a', `b')
			gen `typlist' `varlist'`i' = invnormal(`u'`i') `if' `in'
			label variable `varlist'`i' "Randomized quantile residuals `i'"
		}		
		
		/*Creating standardized quantile residual*/
		
		if "`standardized'" == "standardized" {
		    
			local oldformula = "`e(cmdline)'"
			local preformula = substr("`oldformula'", 8, .)
			local pos = strpos("`preformula'", ",") - 1
			if "`pos'" == "-1" {
				local new_formula = "`preformula'"
			} 
			else {
				local new_formula = substr("`preformula'", 1, `pos')
			}

			quietly glm `new_formula', family(poisson) link(log)
			tempvar hat
			local stand std
			predict double `hat', hat
			
			quietly `oldformula'
			
			foreach i of numlist 1 / `nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
				label variable `varlist'_`stand'`i' "Stand. rand. quantile residuals `i'"
			}
			
		}
	}
	else if "`e(cmd)'" == "glm" & "`e(varfunct)'" == "Poisson" {
		tempvar y n mu p a b u
		gen double `y' = `e(depvar)'
		predict double `mu', mu /*Fitted values in count data n = exp(ir + exposure)*/
		gen double `a' = cond(`y' > 0, poisson(`mu', `y' - 1), 0)
		gen double `b' = poisson(`mu', `y')
		foreach i of numlist 1 / `nqres' {
			gen double `u'`i' = runiform(`a', `b')
			gen `typlist' `varlist'`i' = invnormal(`u'`i') `if' `in'
			label variable `varlist'`i' "Randomized quantile residuals `i'"
		}
		
		/*Creating standardized quantile residual*/
		
		if "`standardized'" == "standardized" {
			tempvar hat
			local stand std
			predict double `hat', hat
			foreach i of numlist 1 / `nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
				label variable `varlist'_`stand'`i' "Stand. rand. quantile residuals `i'"
			}
		}
	}
	
	/*Quantile residual for Negative Binomial Regression 
	(GLM, family = nbinomial, link = any admissible)*/
	
	else if "`e(cmd)'" == "glm" & "`e(varfunct)'" == "Neg. Binomial" {
		tempvar y alphaCam size mu p ones max a b u
		gen double `y' = `e(depvar)'
		predict double `mu', mu /*Fitted values in count data n = exp(ir + exposure)*/
		gen `ones' = 1
		egen  `max' = rowmax(`y' `ones')
		
		local oldformula = "`e(cmdline)'"
		local preformula = substr("`oldformula'", 4, .)
		local pos = strpos("`preformula'", ",") - 1
		local new_formula = substr("`preformula'", 1, `pos')
		
		quietly nbreg `new_formula'
		gen double `alphaCam' = `e(alpha)'
		gen double `size' = 1 / `alphaCam'
		gen double `p' = `size' / (`mu' + `size')
		
		quietly `oldformula'
		
		gen double `a' = cond(`y' > 0, ibeta(`size', `max', `p'), 0)
		gen double `b' = ibeta(`size', `y' + 1, `p')
		foreach i of numlist 1 / `nqres' {
			gen double `u'`i' = runiform(`a', `b')
			gen `typlist' `varlist'`i' = invnormal(`u'`i') `if' `in'
			label variable `varlist'`i' "Randomized quantile residuals `i'"
		}
		
		/*Creating standardized quantile residual*/
		
		if "`standardized'" == "standardized" {
			tempvar hat
			local stand std
			predict double `hat', hat
			foreach i of numlist 1 / `nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
				label variable `varlist'_`stand'`i' "Stand. rand. quantile residuals `i'"
			}
		}
	}
	else if "`e(cmd)'" == "nbreg" {
		tempvar y alphaCam size mu p ones max a b u
		gen double `y' = `e(depvar)'
		gen double `alphaCam' = `e(alpha)'
		gen double `size' = 1 / `alphaCam'
		predict double `mu', n /*Fitted values in count data n = exp(ir + exposure)*/
		gen double `p' = `size' / (`mu' + `size')
		gen `ones' = 1
		egen  `max' = rowmax(`y' `ones')
		gen double `a' = cond(`y' > 0, ibeta(`size', `max', `p'), 0)
		gen double `b' = ibeta(`size', `y' + 1, `p')
		foreach i of numlist 1 / `nqres' {
			gen double `u'`i' = runiform(`a', `b')
			gen `typlist' `varlist'`i' = invnormal(`u'`i') `if' `in'
			label variable `varlist'`i' "Randomized quantile residuals `i'"
		}
		
		/*Creating standardized quantile residual*/
		
		if "`standardized'" == "standardized" {
		    
			local oldformula = "`e(cmdline)'"
			local preformula = substr("`oldformula'", 6, .)
			local pos = strpos("`preformula'", ",") - 1
			if "`pos'" == "-1" {
				local new_formula = "`preformula'"
			} 
			else {
				local new_formula = substr("`preformula'", 1, `pos')
			}
			
			quietly glm `new_formula', family(nbinomial ml) link(log)
			tempvar hat
			local stand std
			predict double `hat', hat
			
			quietly `oldformula'
			
			foreach i of numlist 1 / `nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
				label variable `varlist'_`stand'`i' "Stand. rand. quantile residuals `i'"
			}
		}
	}
	else {
		display "It is not a valid GLM"
	}

end
