# Dialysis Treatment Analysis in R

## Introduction

The objective of this analysis is to investigate how dialysis treatment duration and average weight gain between treatments influence the number of days hospitalized for kidney-failure
patients during one year. Using data from a sample of patients treated at a dialysis facility, the number of days was compared across predefined categories of treatment length (short, long)
and weight gain (slight, moderate, substantial).

## Exploratory Data Analysis

A contingency table of treatment duration by weight gain confirmed that the design is balanced, with 10 observations per cell. Boxplots suggest that days increase with greater weight
gain and are higher for shorter treatments than longer with a main-effects plot supporting
these patterns. An interaction plot suggested the presence of some interaction effect(lines
not completely parallel), but the limited sample size requires the need to formally test for its
presence.

<img width="1247" height="372" alt="image" src="https://github.com/user-attachments/assets/14d5a0c8-5db7-48b1-9297-3d127e0f4641" />

## Two-Way ANOVA With Interaction

The results indicate that weight gain and treatment duration have statistically significant
effects at the 5% level, whereas the interaction term does not.

To understand the difference between the three weight gain categories, a Tukey comparison procedure was performed and the results show that the mean number of hospitalization
days does not differ significantly between slight and moderate weight gain categories. However,substantial weight gain is associated with a significantly higher mean than both slight and
moderate gain.

<img width="948" height="292" alt="image" src="https://github.com/user-attachments/assets/ea844a76-b9ab-4b63-88d0-cfda7f93d610" />

## Model diagnostics and assumption checking

The validity of ANOVA relies on several assumptions: homoscedasticity, normality of the
residuals, independence, and absence of influential outliers. To evaluate any deviations from
the model assumptions, the residuals of the 2-way ANOVA were analyzed.

When plotted against fitted values, the residuals showed an increasing spread as fitted
values increased, indicating non-constant variance. The residuals were also more dispersed
for long than for short durations, and their spread increased with weight gain, violating the
homoscedasticity assumption.

A sequence plot of studentized residuals at lag 1 showed no systematic pattern, with suggested no autocorrelation, consistent with the assumption of independent errors. Normality
of the residuals was assessed via a normal Q-Q plot, which revealed some deviations in the
tails, indicating that the normal-error assumption might not be respected.

<img width="1138" height="451" alt="image" src="https://github.com/user-attachments/assets/934b4fb8-2715-480c-b65d-6a73d4223d55" />

## Alternative Modeling Approaches

**Weighted least squares (WLS regression)** using the LM method with weights constructed
as the inverse of the cell-specific residual variances also helps to correct for heteroscedasticity.
The non-normality is primarily due to deviations at the tails which does not significantly harm
WLS. The model explains 29%-35% of the variation in the data (adjusted R² = 0.29), and the
WLS model ensures reliability by accounting for differences in variability between weight and
treatment groups.

Finally, a general linear models (GLM) approach was explored. A Poisson GLM with log link
was initially fitted, as this is a standard approach for count outcomes such as days hospitalised
([Tamborrino, 2015](https://www.hilarispublisher.com/open-access/count-data-analysis-in-randomised-clinical-trials-2155-6180-1000227.pdf)). However, the residual deviance divided by the residual degrees of freedom
was 3.65, indicating substantial overdispersion, which violates a key Poisson model assumption
([Walker, 2017](https://www.middleprofessor.com/files/applied-biostatistics_bookdown/_book/generalized-linear-models-i-count-data)). A negative binomial GLM, which fits over-dispersed data better (Walker, 2017)
was used instead, with the reference group being moderate weight gain and long treatment
duration.

Since WLS offers superior interpretability over the transformed model (effects are multiplicative) and NB (added dispersion parameters), we can opt for it to explain the predictor
effects.

<img width="1216" height="566" alt="image" src="https://github.com/user-attachments/assets/6e7303d8-5daf-48e4-8f37-9b6cb5283369" />

## Conclusions

The results show that average weight gain is a strong indicator of number of days hospitalized `(F = 11.01, p < 0.001)` with patients having substantial weight gain spending 3.8 days
longer when treatment duration is held constant, on average. The expected number of days
hospitalized was 3.7 for baseline weight and treatment duration (moderate and long). Treatment duration had a mostly non-significant effect `(F = 2.92, p = 0.093)`, indicating that it
did not meaningfully influence the number of days hospitalized once weight was accounted
for on average. There is also no statistically significant interaction effect in the picture here
`(F = 2.22, p = 0.12)`, highlighting that the effect of duration of treatment does not differ
between weight categories.

