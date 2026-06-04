# Event Prediction

Predicts enrollment and events at the design stage using assumed
enrollment and treatment-specific time-to-event models, or at the
analysis stage using blinded or unblinded data and specified enrollment
and time-to-event models through simulations.

## Details

Accurately predicting the date at which a target number of subjects or
events will be achieved is critical for the planning, monitoring, and
execution of clinical trials. The `eventPred` package provides
enrollment and event prediction capabilities using assumed enrollment
and treatment-specific time-to-event models at the design stage, using
blinded or unblinded data and specified enrollment and time-to-event
models at the analysis stage.

At the design stage, enrollment is often specified using a piecewise
Poisson process with a constant enrollment rate during each specified
time interval. At the analysis stage, before enrollment completion, the
`eventPred` package considers several models, including the homogeneous
Poisson model, the time-decay model with an enrollment rate function
\\\lambda(t) = (\mu/\delta) (1 - \exp(-\delta t))\\, the B-spline model
with the daily enrollment rate \\\lambda(t) = \exp(B(t)'\theta)\\, and
the piecewise Poisson model. If prior information exists on the model
parameters, it can be combined with the likelihood to yield the
posterior distribution.

The `eventPred` package also offers several time-to-event models,
including exponential, Weibull, log-logistic, log-normal, piecewise
exponential, model averaging of Weibull and log-normal, spline, and cox.
For time to dropout, the same set of model options are considered. If
enrollment is complete, ongoing subjects who have not had the event of
interest or dropped out of the study before the data cut contribute
additional events in the future. Their event times are generated from
the conditional distribution given that they have survived at the data
cut. For new subjects that need to be enrolled, their enrollment time
and event time can be generated from the specified enrollment and
time-to-event models with parameters drawn from the posterior
distribution. Time-to-dropout can be generated in a similar fashion.

The `eventPred` package displays the Akaike Information Criterion (AIC),
the Bayesian Information Criterion (BIC) and a fitted curve overlaid
with observed data to help users select the most appropriate model for
enrollment and event prediction. Prediction intervals in the prediction
plot can be used to measure prediction uncertainty, and the simulated
enrollment and event data can be used for further data exploration.

The most useful function in the `eventPred` package is `getPrediction`,
which combines model fitting, data simulation, and a summary of
simulation results. Other functions perform individual tasks and can be
used to select an appropriate prediction model.

The `eventPred` package implements a model parameterization that
enhances the asymptotic normality of parameter estimates. Specifically,
the package utilizes the following parameterization to achieve this
goal:

- Enrollment models

  - Poisson: \\\theta = \log(\lambda)\\.

  - Time-decay: \\\theta = (\log(\mu), \log(\delta))'\\.

  - B-spline: \\\log(\lambda(t)) = B(t)' \theta\\, \\B(t) = (B_1(t),
    \ldots, B\_{k+4}(t))'\\ are the B-spline basis with \\k\\ inner
    knots.

  - Piecewise Poisson: \\\theta_j = \log(\lambda_j)\\ for the \\j\\th
    time interval. The left endpoints of time intervals, denoted as
    `accrualTime`, are considered fixed.

- Event or dropout models

  Let \\x\\ denote the covariates for a subject. Let \\\beta\\ denote
  the regression coefficients and \\\sigma\\ denote the scale parameter
  of the AFT model, \$\$\log(T) = \beta' x + \sigma \epsilon.\$\$

  - Exponential: \\\log(\lambda) = \theta' x\\. In other words, \\\theta
    = -\beta\\.

  - Weibull: \\\log(\texttt{weibull scale}) = \theta_1' x\\,
    \\\log(\texttt{weibull shape}) = -\theta_2\\. In other words,
    \\\theta = (\beta', \log(\sigma))'\\.

  - Log-logistic: For the logistic distribution of \\\log(T)\\,
    \\\texttt{location} = \theta_1' x\\, \\\log(\texttt{scale}) =
    \theta_2\\. In other words, \\\theta = (\beta', \log(\sigma))'\\.

  - Log-normal: For the normal distribution of \\\log(T)\\,
    \\\texttt{mean} = \theta_1' x\\, \\\log(\texttt{sd}) = \theta_2\\.
    In other words, \\\theta = (\beta', \log(\sigma))'\\.

  - Piecewise exponential: \\\log(\lambda_j) = \theta\_{1j} + \theta_2'
    x\\ for the \\j\\th time interval, \\\theta = (\theta_1',
    \theta_2')'\\. The left endpoints of time intervals, denoted as
    `piecewiseSurvivalTime` for event model and `piecewiseDropoutTime`
    for dropout model, are considered fixed.

  - Model averaging: \\\theta = (\theta\_{\texttt{weibull}}',
    \theta\_{\texttt{lnorm}}')'\\. The covariance matrix for \\\theta\\
    is structured as a block diagonal matrix, with the upper-left block
    corresponding to the Weibull component and the lower-right block
    corresponding to the log-normal component. In other words, the
    covariance matrix is partitioned into two distinct blocks, with no
    off-diagonal elements connecting the two components. The weight
    assigned to the Weibull component, denoted as \\w_1\\, is considered
    fixed.

  - Spline: Let \\S(t\|x)\\ denote the survival function given
    covariates \\x\\. We model a transformation of the survival function
    as a cubic spine: \$\$g(S(t\|x)) = c(u) + \theta_2' x,\$\$ where
    \$\$c(u) = \gamma_1 + \gamma_2 u + \gamma_3 v_1(u) + \cdots +
    \gamma\_{k+2} v_k(u)\$\$ is the cubic spline in \\u=\log(t)\\,
    \\\theta = (\theta_1', \theta_2')'\\, \\\theta_1 = (\gamma_1,
    \ldots, \gamma\_{k+2})'\\, assuming \\k\\ inner knots (\\k =
    \texttt{knots}\\), and \\v_1(u), \ldots, v_k(u)\\ are the basis of
    the Royston/Parmar spline. The transformation is given as follows:

    - For `scale = "hazard"`, \\g(S(t)) = \log(-\log(S(t)))\\.

    - For `scale = "odds"`, \\g(S(t)) = \log(1/S(t)-1)\\.

    - For `scale = "normal"`, \\g(S(t)) = -\Phi^{-1}(S(t))\\.

    The hazard, odds, and normal scales correspond to extensions of the
    Weibull, log-logistic, and log-normal distributions, respectively.

  - Cox: Let \\t_1 \< \cdots \< t_M\\ denote the distinct observed event
    times, \\\lambda_j\\ denote the estimated baseline hazard rate in
    the \\j\\th time interval, \\(t\_{j-1}, t_j\]\\, and \\\beta\\
    denote the regression coefficients (log hazard ratios) from the Cox
    model. The model parameters including the baseline hazards are
    \\\theta = (\log(\lambda_1),\ldots,\log(\lambda_M), \beta^T)^T\\.

The `eventPred` package uses days as its primary time unit. If you need
to convert enrollment or event rates per month to rates per day, simply
divide by 30.4375.

## References

Emilia Bagiella and Daniel F. Heitjan. Predicting analysis times in
randomized clinical trials. Stat in Med. 2001; 20:2055-2063.

Gui-shuang Ying and Daniel F. Heitjan. Weibull prediction of event times
in clinical trials. Pharm Stat. 2008; 7:107-120.

Xiaoxi Zhang and Qi Long. Stochastic modeling and prediction for accrual
in clinical trials. Stat in Med. 2010; 29:649-658.

Patrick Royston and Mahesh K. B. Parmar. Flexible parametric
proportional-hazards and proportional-odds models for censored survival
data, with application to prognostic modelling and estimation of
treatment effects. Stat in Med. 2002; 21:2175-2197.

## Author

Kaifeng Lu, <kaifenglu@gmail.com>
