{smcl}
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


{title:Title}

{p 4 4 2}
{bf:qresid} {hline 2} (Randomized) quantile residuals for diagnostic of models 



{title:Syntax}

{p 8 8 2} {bf:qresid} {it:newvarname} [{it:{help if}}] [{it:{help in}}] [, {it:options}]


{title:Options}

{col 5}{it:options}{col 24}Description
{space 4}{hline}
{col 5}standardized{col 24}compute leverage standardized quantile residuals
{col 5}nqres(#){col 24}number of randomized quantile residuals; # default is 4 for discrete outcomes; not apply in continuous outcomes
{space 4}{hline}


{title:Description}

{p 4 4 2}
{bf:qresid} calculates (randomized) quantile residuals for the diagnostic of a great variety of models. 

{p 4 4 2}
This version of {bf:qresid} only allows to create quantile residuals for {it:generalized linear models}. Further updates will include {it:zero inflated models}, {it:beta regression} and others {it:related models}. 



{title:Remarks}

{p 4 4 2}
{bf:qresid} implements the computation of quantile residuals proposed by  {browse "http://www.jstor.org/stable/1390802":Peter K. Dunn and Gordon K. Smyth (1996)} for generalized linear models (GLM). 

{p 4 4 2}
The residuals and their plots help to evaluate the assumptions of the statistical models. In the general linear model with normal outcome distribution, the residuals (raw, standardized, and studentized) follow a normal distribution when the model is right. 

{p 4 4 2}
In the case of generalized linear models (GLMs) for non-normal outcome distributions, the deviance and Pearson residuals (and their standardized forms) do not necessarily follow a normal distribution. This makes model diagnostic difficult tin many practical scenarios. As if that were not enough, the GLM residuals for discrete outcomes can present patterns of parallel bands that artifact their proper visual inspection. 



{title:Examples}

{p 4 4 2}
	{bf:1) Outcome with binomial distribution}

{asis}
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
{smcl}
		
{p 4 4 2}
	{bf:2) Outcome with Poisson distribution}

{asis}
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
{smcl}
		
		

{title:Author}

{p 4 4 2}
Percy Soto-Becerra    {break}
InkaStats Data Science Solutions - Medical Branch     {break}
percys1991@gmail.com



{title:References}
{pstd}

{p 4 4 2}
Peter K. Dunn & Gordon K. Smyth. 1996. 
{browse "http://www.jstor.org/stable/1390802":Randomized Quantile Residuals}. Journal of Computational and Graphical Statistics, 5:3, 236-244.

{space 4}{hline}

{p 4 4 2}
This help file was dynamically produced by 
{browse "http://www.haghish.com/markdoc/":MarkDoc Literate Programming package} 


