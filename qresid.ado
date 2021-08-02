*! version 0.0.2 Percy Soto-Becerra 31jul2021
program qresid
	version 17.0
	syntax newvarname(max=1) [if] [in] [, standardized nqres(int 4)]
	/*Quantile residual for normal linear regression (by Ordinary Least Square 
	or Maximum Likelihood estimation*/
	if "`e(cmd)'" == "regress" {
		tempvar residuals
		predict double `residuals', residuals
		gen `typlist' `varlist' = `residuals' / `e(rmse)' `if' `in'
		if "`standardized'" == "standardized" {
			tempvar hat
			local stand  std
			predict double `hat', hat
			gen `typlist' `varlist'_`stand' = `varlist' / sqrt(1 - `hat') 
		}
	}
	else if "`e(cmd)'" == "glm" & "`e(varfunct)'" == "Gaussian" {
		tempvar deviance_res
		predict double `deviance_res', deviance
		gen `typlist' `varlist' = `deviance_res' / `e(dispers)' `if' `in'
		if "`standardized'" == "standardized" {
			tempvar hat
			local stand std
			predict double `hat', hat
			gen `typlist' `varlist'_`stand' = `varlist' / sqrt(1 - `hat') `if' `in'
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
		foreach i of numlist 1/`nqres' {
			gen double `u'`i' = runiform(`a', `b')
			gen `typlist' `varlist'`i' = invnormal(`u'`i') `if' `in'
		}
		if "`standardized'" == "standardized" {
			tempvar hat
			local stand std
			predict double `hat', hat
			foreach i of numlist 1/`nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
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
		foreach i of numlist 1/`nqres' {
			gen double `u'`i' = runiform(`a', `b')
			gen `typlist' `varlist'`i' = invnormal(`u'`i') `if' `in'
		}
		if "`standardized'" == "standardized" {
			tempvar hat
			local stand std
			predict double `hat', hat
			foreach i of numlist 1/`nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
			}
		}
	} 
	/*Quantile residual for Poisson regression 
	(GLM, family = Poisson, link = any admissible)*/
	else if "`e(cmd)'" == "poisson" {
		tempvar y n mu p a b u
		gen double `y' = `e(depvar)'
		predict double `mu', n /*Fitted values in count data n = exp(ir + exposure)*/
		gen double `a' = poisson(`mu', `y' - 1)
		gen double `b' = poisson(`mu', `y')
		foreach i of numlist 1/`nqres' {
			gen double `u'`i' = runiform(`a', `b')
			gen `typlist' `varlist'`i' = invnormal(`u'`i') `if' `in'
		}		
		
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
			
			foreach i of numlist 1/`nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
			}
		}
	}
	else if "`e(cmd)'" == "glm" & "`e(varfunct)'" == "Poisson" {
		tempvar y n mu p a b u
		gen double `y' = `e(depvar)'
		predict double `mu', mu /*Fitted values in count data n = exp(ir + exposure)*/
		gen double `a' = poisson(`mu', `y' - 1)
		gen double `b' = poisson(`mu', `y')
		foreach i of numlist 1/`nqres' {
			gen double `u'`i' = runiform(`a', `b')
			gen `typlist' `varlist'`i' = invnormal(`u'`i') `if' `in'
		}
		if "`standardized'" == "standardized" {
			tempvar hat
			local stand std
			predict double `hat', hat
			foreach i of numlist 1/`nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
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
		foreach i of numlist 1/`nqres' {
			gen double `u'`i' = runiform(`a', `b')
			gen `typlist' `varlist'`i' = invnormal(`u'`i') `if' `in'
		}
		if "`standardized'" == "standardized" {
			tempvar hat
			local stand std
			predict double `hat', hat
			foreach i of numlist 1/`nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
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
		foreach i of numlist 1/`nqres' {
			gen double `u'`i' = runiform(`a', `b')
			gen `typlist' `varlist'`i' = invnormal(`u'`i') `if' `in'
		}
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
			
			foreach i of numlist 1/`nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
			}
		}
	}
	else {
		display "It is not a valid GLM"
	}
end