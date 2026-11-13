# Simulation protocol

## Design and generating mechanism

The unit of sampling is an individual time-stratified case-crossover matched set. Candidate dates are in the same year, month and day of week. A single event date is drawn from conditional logistic probabilities proportional to exp{h(X)+g(W)} within each set. The generated set, its case indicator and all of its reference dates remain linked through pair expansion, cross-fitting, estimation and covariance aggregation.

Eight active exposures are PM2.5, coarse particulate matter, O3, NO2, SO2, CO, temperature and relative humidity. `config/mechanism.R` specifies all exposure-distribution, response, confounding, scaling and reference constants. `R/02_calendar.R` generates temporally correlated candidate measurements and fixed-lag summaries; `R/10_mechanism.R` maps concentrations to the specified ranges and samples the event date. Fixed lags are encoded per variable as 0:1, 0:1, 0:1, 1, 0:3, 1, 0:3 and 0:3 days, respectively.

Pollutant main responses use the specified concentration/percentage anchor points, monotone Hermite interpolation between anchors and smooth tail continuation, followed by the log1p transformation and fixed coefficients. Temperature uses separate quadratic cold and heat responses around 23 degrees C. Humidity acts through its interaction with temperature.

Let Z_j=(X_j-center_j)/scale_j, using the constants in the mechanism file. The interaction response is:

    0.135743242415130250 Z_temperature Z_humidity
  + 0.084016234171294737 Z_temperature Z_PM2.5
  + 0.085641763170191634 Z_temperature Z_O3
  + 0.069229207135094992 Z_PM2.5 Z_O3.

The confounder response is 0.35 W + 0.8 W^2. Null exposures are generated independently of the response mechanism, with variance one, a shared matched-set component of variance 0.20 and a within-set component of variance 0.80. Increasing candidate dimension preserves the active measurements and existing null columns within a repetition. Module C uses 8,16,40,80 null exposures, giving total dimensions 16,24,48,88.

## Estimation

An independent learning sample determines standardization, coverage-enhanced Gaussian kernel centers and outcome-supervised distance parameters. The fitted exposure surface is linear in whitened kernel features. Each case/reference pair is represented in both orientations with matched-set-normalized weights. Nuisance functions g and f are cross-fitted by complete matched sets.

The coefficient estimator minimizes 0.5 ||T U_N(beta;g,f)||^2 + 0.5 lambda beta'P beta. T is a fixed learning-sample inverse-Jacobian transformation and P is the normalized integrated mixed-derivative penalty. A damped Gauss–Newton algorithm uses backtracking and both gradient and scaled-step stopping criteria. Lambda remains fixed during an N comparison. The finite-basis regularized target and sampling accuracy are evaluated jointly in simulation.

Covariance uses exact objective curvature, matched-set score aggregation and fitted-nuisance influence contributions, conditional on the independently learned kernel. Model SE, empirical SD across fresh repetitions, bias and interval coverage characterize its empirical uncertainty calibration.

## Scenarios and comparators

| Scenario | g working family | f working family |
|---|---|---|
| G1 | Contains the response and W,W² terms | Nine ordered confounder-pair cells |
| G2 | Contains the response and W,W² terms | Pooled cells |
| G3 | W² omitted | Nine ordered confounder-pair cells |
| G4 | W² omitted | Pooled cells |

The coefficients of each nuisance model are estimated from data. Linear-CLR uses linear exposure terms; Spline-CLR uses additive natural cubic splines with internal learning-sample quantiles 0.25,0.50,0.75 and boundary quantiles 0.005,0.995. Each comparator is fitted with W,W² and with W adjustment. Since CLR uses the outcome adjustment alone, the appropriate same fit supplies G1/G2 or G3/G4 displays. KDR-CC permits interactions; additive CLR produces an exactly zero four-point log-OR interaction contrast. Interpretation of method comparisons distinguishes response representation from nuisance robustness.

## Effect targets

A fixed independent evaluation sample of 5,000 matched sets supplies marginal concentration quantiles and the median reference vector b. The observed-combination distribution uses one candidate row from each of the first 1,000 sets.

- Joint contrast: h(Q_1(q),...,Q_p(q))-h(b), over q=0.10 to 0.90.
- Observed-combination contrast: h(x)-h(b), at the saved sampled exposure combinations.
- Single-exposure contrast: h(a,b_-j)-h(b), with other exposures at their reference.
- Pairwise interaction: h(a,c,b_-jk)-h(a,b_k,b_-jk)-h(b_j,c,b_-jk)+h(b).

All effects are log OR contrasts; exponentiation gives OR. Identically zero anchor contrasts are excluded from aggregate error and significance summaries. The single and interaction grids span the saved joint-quantile ranges. Constructed reference-background contrasts and observed-combination contrasts use explicitly different evaluation coordinates.

For error e_rj=estimate_rj-truth_j, summarize RMS bias=sqrt(mean_j(mean_r(e_rj)^2)), RMSE=sqrt(mean_rj(e_rj^2)), RMS empirical SD, RMS model SE and their ratio. Absolute OR relative error is 100*abs(exp(e_rj)-1); report mean, median, quartiles and the 95th percentile. Correct-direction significant detection requires the nominal interval to exclude zero on the true-effect side; wrong-direction detection is also recorded. Pointwise coverage and interval width are secondary uncertainty measures. Monte Carlo uncertainty uses whole independent repetitions.

Single-exposure figures show repetition-mean predictions. Light shading connects the separately averaged lower and upper pointwise 95% confidence limits. Darker shading is the empirical 2.5th–97.5th percentile range of repeated point estimates. These summarize typical model uncertainty and repeated-study variability. Interval availability is tabulated separately.

Quantitative OR summaries include joint P25/P75/P90 versus P50 and P75 versus P25, single-exposure P75 versus P25, and temperature P10/P90 versus the specified optimum. Other exposures are fixed at b for single-variable contrasts. The full contrast covariance determines the interval, including covariance between the two predictions.

## Contribution shares

Let A contain the eight prespecified active exposures and b denote the marginal-median reference. Baseline Shapley values allocate Delta_A(x)=h(x_A,b_other)-h(b) over A and satisfy sum_{j in A} phi_j(x)=Delta_A(x). Additional coordinates stay at b throughout every coalition. Importance is A_j=mean_x |phi_j(x)| on the active coordinates of the fixed observed-combination sample; shares normalize these quantities within A. This sign-free allocation includes main and interaction contributions at the specified background.

For a separable Gaussian kernel, factors for the additional coordinates are evaluated at b and absorbed into the center coefficients. Interpolating each active-coordinate factor between baseline and observed levels gives a degree-seven polynomial for each Shapley marginal contribution. Four Gauss–Legendre nodes integrate it exactly up to floating-point error. Prefix and suffix products avoid division by small coordinate factors. Additive CLR contributions follow coordinate-wise contrasts. The known main-plus-bilinear truth has an analytic Shapley allocation, verified against subset enumeration.

Contribution accuracy compares the eight fitted shares with the eight true shares. Each repetition yields a share RMSE over those eight coordinates, followed by aggregation over independent repetitions. All contribution bars contain the same eight active exposures and sum to 100%. The active set is specified before fitting and remains the same in every repetition and dimension setting.

## Interaction detection

The four specified pairs form a Holm family at alpha=0.05. For a pair, D is the four-point contrast operator across its non-anchor grid and the overall statistic is mean((D beta_hat)^2). Under a centered Gaussian coefficient approximation, its null reference is a weighted sum of independent chi-square(1) variables with weights from the contrast covariance. Davies integration evaluates the upper tail. Ambiguous numerical decisions are refined; unresolved tests remain unavailable.

Two-sided pointwise Wald tests assess direction. The corrected version applies Holm across all non-anchor grid points of all four pairs. A correct significant direction requires a nonzero true contrast, rejection and a matching sign. Rates are averaged within each pair and repetition, then across repetitions. Overall pair detection and direction recovery have separate tables and figures; wrong-direction rates are retained. Approximate conditional inference is evaluated through empirical operating characteristics under the specified regularization.

## Exposure identification

The derivative-design quadratic statistic beta_hat' Q_j beta_hat tests each coordinate's participation in the response. Its numerical chi-square-mixture reference uses the coefficient covariance. Holm correction covers all p candidate exposures within a fit. Interaction-only humidity is active because its response derivative is nonzero over the evaluation distribution. The numeric-tail calculation uses accuracy 1e-8, with refined integration and a seeded dimension-adaptive fallback where required. The fallback settings and calculation status are saved.

C reports the selected fraction of eight active exposures, the selected fraction of q null exposures, and the identification frequency of each active variable. Whole-repetition rates determine Monte Carlo intervals; per-variable frequencies use Wilson intervals. Missing tests have explicit availability counts.

## Modules and reproducibility

A uses 100 repetitions for each of 50 configurations. Its assessment rule is given in PARAMETERS.md. B uses a frozen configuration and 100 repetitions under each of eight distinct fits. C uses that configuration under G1 at four null counts, with 100 repetitions each and the complete B-style evaluation outputs.

A's random namespace is independent of B/C. N and L comparisons use prefixes of estimation and independent learning pools within a repetition. Kernels and nuisance estimates are shared only for identical inputs. Semantic seeds are independent of worker scheduling. Every fit stores numerical convergence, coefficient covariance, effects, contribution shares, interactions and test availability. Reporting reads these records separately from fitting.

## Fixed-configuration correlation assessment

Module D reuses completed B fits as the original-correlation reference and adds moderate/high dependence among eight active exposures. Mean absolute Pearson targets are 0.60 and 0.85; all 28 pairs must meet absolute lower bounds of 0.45 and 0.75. Each new setting has 100 repetitions with eight B-style fits, giving 1,600 new fits. Fixed empirical marginal maps and Gaussian common scores recouple final lag-defined exposures. Independent samples calibrate and validate achieved correlations before fitting. Marginal distributions are approximately preserved; temporal and exposure–confounder dependence change under the specified recoupling. Response coefficients and the matched-set structure remain fixed.

B supplies the exact reference concentrations and constructed joint, single-exposure, OR and interaction contrasts. Observed-combination effects and contribution shares use each setting's own evaluation distribution. B and D batches are independent, while D settings share random inputs within a repetition. Complete B snapshots, integrity hashes, calibration parameters and evaluation coordinates are saved. Details and submission commands are in MODULE_D.md.

## Data-driven configuration selection

Module E uses seven exposure settings and 100 independent repetitions per setting. An independent learning sample builds candidate kernel representations. Three-fold conditional likelihood validation on a separate configuration pool compares 18 full KDR-CC fits, using only exposure, covariate and case labels. Complete matched sets remain within folds. Nuisance exposure features comprise learning-data natural splines and retained kernel columns; confounder cells provide the specified G1 adjustment. This is an estimated nuisance approximation under full adjustment. All-fold availability is required for candidate eligibility; exact ties use a deterministic order. The fixed configuration is the explicit fallback when all candidates fail.

Each chosen configuration and the original fixed configuration are fitted to the same fresh final sample. Simulation truth enters only subsequent evaluation. Adaptive results are aggregated over the configurations selected in all repetitions, including fallback repetitions, with availability reported separately. Joint curves use common eight-active-exposure coordinates; observed-combination effects use scenario-specific evaluation distributions. Module E is independent of A–D random batches. MODULE_E.md provides the complete role allocation, candidate table, evaluation and computational plan.
