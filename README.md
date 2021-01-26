# likelihoodR Package

The likelihood approach is one of several approaches to making
inferences from data. The best description and justification for the
approach is given by [Edwards, A.W.F. (1992) Likelihood, Johns Hopkins
Press](https://www.amazon.co.uk/Likelihood-W-F-Edwards/dp/0801844436/). Others, such as [R. Royall](https://www.amazon.co.uk/Statistical-Evidence-Likelihood-Monographs-Probability/dp/0412044110/), [S. Goodman](https://ajph.aphapublications.org/doi/abs/10.2105/AJPH.78.12.1568), [Z. Dienes](https://www.amazon.co.uk/Understanding-Psychology-Science-Introduction-Statistical/dp/023054231X/), [S. Glover & P. Dixon](https://link.springer.com/article/10.3758/BF03196706) have subsequently made
important contributions.

The likelihood approach focusses on the observed data, using maximum
likelihood for estimates, and calculates likelihood ratios for specific
parameter values given the collected data. The log of a likelihood ratio
is known as the *support*. This statistic has distinct advantages. It 
represents the weight of evidence with a scale that ranges
from positive to negative (indicating support for or against a hypothetical 
parameter value). Support values obtained from independent studies
can simply be added together to give their combined support. Unlike *p* values, 
support values are insensitive to transforms. 

Support values give the relative strength of evidence for one hypothesis 
value versus another. They range from negative infinity to positive infinity. 
A support of 0 shows no support for either hypothesis, while a support of 1 
represents weak evidence for the first hypothesis versus the second. 
A support of 2 represents moderate evidence, and is roughly equivalent to the 
5% significance level in the frequentist approach. A support of 3 represents
strong evidence and 4 represents extremely strong evidence, and so on. 
Negative values represent the same strength of evidence, except for the 
second hypothesis value versus the first. Support values for a given likelihood
function can be compared with each other simply by subtracting them. In the 
example below, the support for d versus the 2nd hypothesis can be obtained 
by subtracting the support for the 2nd hypothesis versus null (4.115) from the support
for d versus the null (2.191) giving -1.924, as listed in the 3rd line of output.
The support for d versus the observed mean would be 2.191 - 4.352 = -2.161, negative
as expected since d lower on the likelihood function curve from the maximum
represented by the observed mean (MLE). 

The likelihood interval can be obtained for many statistics. This represents 
the range of values that are consistent with the observed data for a given level
of support. For example, the S-2 likelihood interval includes all values which are 
not more different from the MLE by more than a support
value of 2. As noted for the figure below, the S-2 interval is numerically closely equivalent 
to the 95% confidence interval, as the support of 2 is similar to the 5% significance 
level.
An S-3 likelihood interval would include all values which are not more different 
from the MLE by more than a support value of 3, and so on. The stronger 
the interval the wider it will be. The interpretation of the likelihood interval 
is distinct from a confidence interval.
The confidence interval represents the long run probability of capturing the
population parameter and may need to be corrected for multiple testing, stopping
rule, etc. The likelihood interval is also distinct from the Bayesian credibility
interval that represents the subjective probability for a population value occurring
within it.

There are few statistical packages that implement the likelihood
approach and which calculate support. I have created an R package called
likelihoodR which calculates support for a range of statistical
analyses: *t* tests, ANOVA, regression and categorical analyses. In R it
can be installed from my GitHub repository using the **devtools** command:
``` r
install.packages("devtools")
devtools::install_github("PeterC-alfaisal/likelihoodR")
``` 
The functions in the package complement my recent book: [Cahusac, P.M.B.
(2020) Evidence-Based Statistics, Wiley](https://onlinelibrary.wiley.com/doi/book/10.1002/9781119549833)  
[Amazon.co.uk](https://www.amazon.co.uk/Evidence-Based-Statistics-Introduction-Evidential-Statistical/dp/1119549809/)  

I would be interested in feedback <pcahusac@alfaisal.edu>          Peter Cahusac  

Once installed, full details about each function, including arguments, outputs, and examples, can be obtained in R by typing
?function_name. For example, for the 2 independent samples *t* test function: 
``` r
?L_2S_ttest
```
### **Example for this function (2 independent samples *t* test), using a variation on Gosset's original additional hours of sleep data**
``` r
mysample <- c(0.7, -1.6, -0.2, -1.2, -0.1, 3.4, 3.7, 0.8, 0.0, 2.0)
treat <- rep(1:0,each=5)
L_2S_ttest(mysample, treat, veq=0, null=0, d=0.5, alt.2=2, L.int=2)
```
    ## Maximum support for the observed mean 2.46 (dashed line) against the null 0 (black line) = 4.352
    ##  Support for d of 0.5 (0.9245269, blue line) versus null = 2.191
    ##  Support for d versus 2nd alt Hypothesis 2 (green line) = -1.924
    ##  Support for 2nd alt Hypothesis versus null = 4.115
    ##
    ##  S-2 likelihood interval (red line) is from 0.99557 to 3.92443
<figure>
<img src="https://github.com/PeterC-alfaisal/likelihoodR/blob/master/Rplot.jpeg" id="id" class="class" style="width:60.0%;height:60.0%" />
</figure> Plot showing likelihood function around the observed mean difference (dashed line). Two additional coloured lines are shown for other possible hypotheses. The horizontal red line indicates the likelihood interval for a support value of 2 (closely equivalent to the 95% confidence interval).  


### **All the available functions**
``` r
# One sample and related samples t test
L_ttest(data1, data2, null=0, d=0.5, alt.2=NULL, L.int=2) 

# Independent samples t test
L_2S_ttest(data, group, veq=0, null=0, d=0.5, alt.2=NULL, L.int=2)

# Sample size calculation using the evidential approach for t tests
L_t_test_sample_size(MW = 0.05, sd = 1, d = 1.2, S = 3, paired = FALSE)

# One-way independent samples ANOVA
L_1way_ANOVA(data, group, cont1=NULL, cont2=NULL)

# Two-way independent samples factorial ANOVA
L_2way_Factorial_ANOVA(data, factor1, factor2, contrast1=NULL, contrast2=NULL)

# One-way repeated measures ANOVA
L_1way_RM_ANOVA(dat, group, ID)

# Correlation
L_corr(xv, yv, null=0, exp.r=NULL, L.int=2, alpha=.05)

# Regression, comparing linear, quadratic and cubic fits
L_regress(y, x)

# Multiple logistic regression
L_logistic_regress(yv, p1, p2=NULL, p3=NULL, p4=NULL, p5=NULL, p6=NULL)

# One-way categorical data, goodness of fit
L_1way_cat(obs, exp.p=NULL, L.int=2, alpha=0.05, toler=0.0001)

# Two-way categorical data, tests for interaction and main effects
L_2way_cat(table)

# Odds Ratio
L_OR(table, null=1, exp.OR=NULL, L.int=2, alpha=0.05, cc=FALSE, toler=0.0001)

# Relative Risk
L_RR(table, null=1, exp.RR=NULL, L.int=2, alpha=0.05, cc=FALSE, toler=0.0001)

# Efficacy
L_efficacy(a, n, null=0, exp.eff=NULL, L.int=2, alpha=0.05, toler=0.0001)

# Likelihood-based confidence interval for the binomial
binpL(a, n, alpha, toler)
```
