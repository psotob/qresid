*! version 0.0.4 Percy Soto-Becerra 06ago2021
program define qresid
	version 17.0
	syntax newvarname(max=1) [if] [in] [, standardized nqres(int 4) diag_glm]
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
			
			if "`diag_glm'" == "diag_glm" {
			    display "regress is not supported yet, use glm ..., family(gaussian) link(identity)"
			}
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
			
			if "`diag_glm'" == "diag_glm" {
			    tempvar mu working
				predict double `mu' if e(sample), mu
				gen `typlist' mu_scaled = `mu'
				predict double xb if e(sample), xb
				predict double `working' if e(sample), working
				gen `typlist' workresp = `working' + xb
				predict `typlist' cooksd if e(sample), cooksd
				predict `typlist' dev_student if e(sample), deviance studentized 
				predict `typlist' hat_lever if e(sample), hat
				lab var workresp "Working responses"
				lab var mu_scaled "Fitted values (on const. inf. scale)"
				lab var hat_lever "Leverage (h_i of Hat matrix's diagonal)"
			}
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
			
			if "`diag_glm'" == "diag_glm" {
			    tempvar mu working
				predict double `mu' if e(sample), mu
				gen `typlist' mu_scaled = sin(sqrt(`mu')) ^ (-1)
				predict double xb if e(sample), xb
				predict double `working' if e(sample), working
				gen `typlist' workresp = `working' + xb
				predict `typlist' cooksd if e(sample), cooksd
				predict `typlist' dev_student if e(sample), deviance studentized 
				predict `typlist' hat_lever if e(sample), hat
				lab var workresp "Working responses"
				lab var mu_scaled "Fitted values (on const. inf. scale)"
				lab var hat_lever "Leverage (h_i of Hat matrix's diagonal)"
			}
			
			quietly `oldformula'
			
			foreach i of numlist 1 / `nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
				label variable `varlist'_`stand'`i' "Standardized randomized quantile residuals `i'"
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
				label variable `varlist'_`stand'`i' "Standardized randomized quantile residuals `i'"
			}
			if "`diag_glm'" == "diag_glm" {
			    tempvar mu working
				predict double `mu' if e(sample), mu
				gen `typlist' mu_scaled = sin(sqrt(`mu')) ^ (-1) 
				predict double xb if e(sample), xb
				predict double `working' if e(sample), working
				gen `typlist' workresp = `working' + xb
				predict `typlist' cooksd if e(sample), cooksd
				predict `typlist' dev_student if e(sample), deviance studentized 
				predict `typlist' hat_lever if e(sample), hat
				lab var workresp "Working responses"
				lab var mu_scaled "Fitted values (on const. inf. scale)"
				lab var hat_lever "Leverage (h_i of Hat matrix's diagonal)"
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
			
			if "`diag_glm'" == "diag_glm" {
			    tempvar mu working
				predict double `mu' if e(sample), mu
				gen `typlist' mu_scaled = sqrt(`mu')
				predict double xb if e(sample), xb
				predict double `working' if e(sample), working
				gen `typlist' workresp = `working' + xb
				predict `typlist' cooksd if e(sample), cooksd
				predict `typlist' dev_student if e(sample), deviance studentized 
				predict `typlist' hat_lever if e(sample), hat
				lab var workresp "Working responses"
				lab var mu_scaled "Fitted values (on const. inf. scale)"
				lab var hat_lever "Leverage (h_i of Hat matrix's diagonal)"
			}
			
			quietly `oldformula'
			
			foreach i of numlist 1 / `nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
				label variable `varlist'_`stand'`i' "Standardized randomized quantile residuals `i'"
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
				label variable `varlist'_`stand'`i' "Standardized randomized quantile residuals `i'"
			}
			if "`diag_glm'" == "diag_glm" {
			    tempvar mu working
				predict double mu if e(sample), mu
				gen `typlist' mu_scaled = sqrt(`mu')
				predict double xb if e(sample), xb
				predict double `working' if e(sample), working
				gen `typlist' workresp = `working' + xb
				predict `typlist' cooksd if e(sample), cooksd
				predict `typlist' dev_student if e(sample), deviance studentized 
				predict `typlist' hat_lever if e(sample), hat
				lab var workresp "Working responses"
				lab var mu_scaled "Fitted values (on const. inf. scale)"
				lab var hat_lever "Leverage (h_i of Hat matrix's diagonal)"
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
				label variable `varlist'_`stand'`i' "Standardized randomized quantile residuals `i'"
			}
			if "`diag_glm'" == "diag_glm" {
			    tempvar mu working
				predict double `mu' if e(sample), mu
				gen `typlist' mu_scaled = `mu'
				predict double xb if e(sample), xb
				predict double `working' if e(sample), working
				gen `typlist' workresp = `working' + xb
				predict `typlist' cooksd if e(sample), cooksd
				predict `typlist' dev_student if e(sample), deviance studentized 
				predict `typlist' hat_lever if e(sample), hat
				lab var workresp "Working responses"
				lab var mu_scaled "Fitted values (on const. inf. scale)"				
				lab var hat_lever "Leverage (h_i of Hat matrix's diagonal)"
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
			
			if "`diag_glm'" == "diag_glm" {
			    tempvar mu working
				predict double `mu' if e(sample), mu
				gen `typlist' mu_scaled = `mu'
				predict double xb if e(sample), xb
				predict double `working' if e(sample), working
				gen `typlist' workresp = `working' + xb
				predict `typlist' cooksd if e(sample), cooksd
				predict `typlist' dev_student if e(sample), deviance studentized 
				predict `typlist' hat_lever if e(sample), hat
				lab var workresp "Working responses"
				lab var mu_scaled "Fitted values (on const. inf. scale)"
				lab var hat_lever "Leverage (h_i of Hat matrix's diagonal)"
			}
			
			quietly `oldformula'
			
			foreach i of numlist 1 / `nqres' {
			    gen `typlist' `varlist'_`stand'`i' = `varlist'`i' / sqrt(1 - `hat') `if' `in'
				label variable `varlist'_`stand'`i' "Standardized randomized quantile residuals `i'"
			}
		}
	}
	else {
		display "It is not a valid GLM"
	}

end