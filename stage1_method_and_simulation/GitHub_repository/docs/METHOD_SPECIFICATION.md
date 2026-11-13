# KDR-CC method specification

Version 1.0.0. The mathematical definitions below specify the distributed implementation and the A–E study. The Chinese appendix contains the derivations, proofs, response constants and numbered equations. The generating-mechanism identity is `9d3750065660d0322c16b0fd90ff134c`.

The appendix uses symbolic dimensions and tuning settings in S.1–S.9. The finite-dimensional auxiliary formulation includes the implemented conditional-likelihood regression and ordered-cell conditional means as special cases. Simulation-specific dimensions, numerical settings and A–E designs appear in S.10–S.12.

## 1. Observation and response

The independent unit is a complete individual matched set with one event date and four or five same-year, same-month, same-weekday candidate dates. Conditional on the candidate covariates,

\[
P(J_i=j\mid C_i)=\frac{\exp\{h_0(A_{ij})+q_0(W_{ij})\}}{\sum_k\exp\{h_0(A_{ik})+q_0(W_{ik})\}}.
\]

An effect is the conditional contrast `h(a)-h(b)` on the log OR scale; exponentiation gives the OR. Every split, validation loss and covariance contribution preserves complete matched sets. The simulation's eight active exposures are PM2.5, coarse PM, O3, NO2, SO2, CO, temperature and humidity. Fixed lag windows are respectively 0:1, 0:1, 0:1, 1, 0:3, 1, 0:3, 0:3 days. `config/mechanism.R` defines the complete concentration maps, response anchors, temperature function, four bilinear interactions and confounding function.

## 2. Ordered pairs and base moment

Each event/reference pair is included in both orientations, with weight `1/[2(n_i-1)]`. The weights sum to one within each original matched set. Let D indicate event-first orientation, Z be the first-minus-second feature difference, and X the ordered covariate pair. Symmetric pairing preserves the conditional direction odds `exp(beta0'Z+g0(X))` when the true response is in the fixed feature space.

\[
\psi(\beta,g,f)=(Z-f(X))\{D\exp(-Z^T\beta-g(X))-(1-D)\}.
\]

All expectations refer to the probability measure formed by matched-set-normalized pair weights. With `f0(X)=E[Z|D=0,X]` and `p0(X)=P(D=0|X)`,

\[
E\psi(\beta_0,g,f)=E[p_0(X)\{\exp(g_0(X)-g(X))-1\}\{f_0(X)-f(X)\}].
\]

This is the double-robust identity. A correct finite response representation, integrability and a locally invertible moment Jacobian connect its zero to an identified coefficient. At both correct nuisances, the Jacobian is the negative weighted conditional residual covariance of Z. A finite representation error enters the moment through its exponential remainder. Independent-learning and finite-sample performance are evaluated against the generating response in simulation.

## 3. Independent kernel learning

An independent sample of L matched sets supplies geometry, kernel coefficients used for distance learning, the feature transform and the final penalty/weight matrices.

- Geometry weights are `1/(L*n_i)` per row. Coordinates use weighted medians and weighted standard deviations.
- The first approximately 75% of M centers use weighted k-means++ and weighted Lloyd iterations (two starts, 25 iterations, tolerance 1e-5). Remaining centers use farthest-point coverage among observations inside coordinate-wise 0.5%–99.5% quantiles when that candidate set is large enough.
- Features are `K_rho(z,u_m)-K_rho(0,u_m)`, with `K=exp(-sum_q rho_q*(z_q-u_mq)^2)`.
- The geometry scale is the reciprocal median positive squared distance to centers. All distances start at 0.0625 times this scale and have bounds 0.01–100 times the same scale.
- For fixed distances, minimize matched-set conditional negative log likelihood divided by L plus `kappa/2` times the squared norm of all kernel and W/W² coefficients. Newton/backtracking solves the inner problem. Profile gradients and L-BFGS-B learn log distances (80 iterations per block, one continuation block, projected-gradient tolerance 1e-5).
- Within-set-centered features undergo full-rank pivoted QR whitening. Rank tolerance is 1e-12. The same invertible transform applies to coefficients, gradients and penalties.
- On 256 learning rows, average outer products of all pairwise mixed coordinate derivatives. Transform to the whitened basis, symmetrize and normalize its trace to M, giving Ω.
- The learning-sample pair Jacobian at zero final coefficients, with full learned g and nine-cell f, supplies T=J_L(0)^{-1}. The fixed moment metric is W0=T'T.

Kernel parameters stay strictly positive in this implementation. Variable participation is evaluated by the fitted-function derivative test in Section 8.

## 4. Auxiliary estimation

Three-fold cross-fitting uses whole matched sets. `g(X)=dW'gamma`, where dW comprises W differences and squared-W differences. Gamma is extracted from a joint conditional-likelihood fit of exposure nuisance features and covariates; its covariance calculation retains all fitted nuisance coefficients.

A–D's complete exposure dictionary contains the six generating pollutant bases, separate cold/hot temperature squares, humidity square, standardized linear coordinates, all 28 bilinear pairs and linear/quadratic null coordinates. All coefficients are estimated. G1/G2 include W²; G3/G4 omit W². The f estimate is the pair-weighted mean of Z among D=0 records: nine ordered W-pair cells in G1/G3, one pooled cell in G2/G4.

E constructs nuisance exposure features from learning-data natural splines (interior quantiles .25/.5/.75, boundaries .005/.995) and up to 24 whitened kernel columns. It uses W/W² and nine-cell f. Both E strategies use this estimated approximation. Its G1 denotes complete covariate adjustment with the same f specification.

## 5. Implemented coefficient estimator

Let R=Z-f, e=D exp(-Z'beta-dW'gamma), t=e-(1-D). With cross-fitted nuisances,

\[
U_N=N^{-1}\sum_{ir}w_{ir}R_{ir}t_{ir},\qquad
J_N=-N^{-1}\sum_{ir}w_{ir}e_{ir}R_{ir}Z_{ir}^T.
\]

The coefficient criterion and stationary equation are

\[
Q_N(\beta)=\tfrac12U_N^TW_0U_N+\tfrac\lambda2\beta^T\Omega\beta,
\qquad F_N=J_N^TW_0U_N+\lambda\Omega\beta=0.
\]

Each job sets `reference_N` to its own N. Thus `lambda_eff=lambda*reference_N/N` equals the configured lambda, including N comparisons. The fixed-penalty population target is the corresponding locally unique minimizer βλ. Under the base DR conditions, `F(beta0)=lambda*Omega*beta0`. Locally,

\[
\beta_\lambda-\beta_0=-\lambda(J_0^TW_0J_0)^{-1}\Omega\beta_0+o(\lambda).
\]

The base double-robust zero and fixed-penalty target are separate mathematical objects. Consistency for βλ follows from local identification and uniform convergence; consistency for β0 additionally requires appropriate penalty decay or a zero penalty gradient at the truth. The simulation evaluates the complete fixed-penalty procedure against h0.

Optimization starts at zero. A damped Gauss–Newton direction uses `J'W0J+lambda*Omega`, with damping candidates 0, 1e-10, 1e-8, 1e-6, 1e-4, .01 times its scale. Armijo backtracking uses the actual criterion, coefficient 1e-4 and up to 20 halvings. The maximum is 70 iterations; gradient and scaled-step tolerances are 1e-7 and 1e-6. Objective, moment residual, gradient, step and damping are saved.

## 6. Conditional covariance

The exact criterion Hessian includes curvature from the moment residual:

\[
H=J^TW_0J+N^{-1}\sum_{ir}w_{ir}e_{ir}\{R_{ir}^TW_0U\}Z_{ir}Z_{ir}^T+\lambda\Omega.
\]

The sandwich score aggregates the direct empirical-moment and empirical-Jacobian contributions with the influence of every fold's fitted g and f. Covariate scores use the complete auxiliary conditional-likelihood information matrix before extracting gamma. Cell means use ratio-estimator influences. The factor `N/N_train` maps each training-fold influence to the final scale. All contributions of one original matched set are summed before taking outer products.

\[
\widehat V=H^{-1}\left\{\frac{\sum_i(S_i-\bar S)(S_i-\bar S)^T}{N(N-1)}\right\}H^{-T}.
\]

With fixed dimension, balanced independent folds, positive cell masses, nonsingular auxiliary information, locally positive exact curvature, integrable derivatives and numerical residual o_p(N^(-1/2)), finite-dimensional smooth-M-estimation yields asymptotic normality around βλ conditional on independently learned objects. The appendix gives the full sensitivity matrices and proof. Positive exact curvature and numerical conditioning are checked. A linear contrast c uses variance c'Vc and pointwise normal intervals; OR intervals exponentiate the limits. Cross-repetition empirical variability includes the learning stage, whereas each model covariance conditions on that learned kernel and, in E, selected configuration. Truth coverage jointly reflects sampling, representation and regularization errors.

## 7. Scientific outputs

- **Joint curve:** all candidate exposures move together to their marginal q quantiles. True response depends on eight active coordinates; null coordinates still move in the fitted input.
- **Observed-combination effects:** h(a)-h(b) on 1,000 fixed independent observed combinations. Both the point and reference enter the contrast covariance.
- **Single-exposure curves:** replace one coordinate of median b, holding the others at b. Show mean estimates, separately averaged pointwise CI limits, and empirical 2.5%–97.5% point-estimate quantiles.
- **Quantitative OR:** joint P25/P75/P90 versus P50 and P75 versus P25; single P75 versus P25; temperature P10/P90 versus 23 °C.
- **Interaction:** the four-point difference `h(x,y)-h(x,bk)-h(bj,y)+h(b)`. Positive and negative contrasts denote departure above and below additivity on the log OR scale.
- **Contribution:** define the active set A as the eight prespecified generating exposures. Baseline Shapley values allocate h(a_A,b_other)-h(b), with all additional coordinates fixed at their marginal-median reference. Average absolute allocations over the active coordinates of the observed-combination evaluation sample, then normalize within A. The coalition set is fixed by the simulation design. Gaussian factors of the reference-background coordinates multiply each kernel-center coefficient; four Gauss–Legendre nodes exactly integrate the remaining degree-seven polynomial. Allocations sum to the active-exposure contrast. Variable-identification and full-input prediction targets are evaluated separately.

## 8. Tests and numerical tails

Variable participation uses the mean squared derivative statistic `beta'Q_j beta`, with Qj formed from feature derivatives with respect to standardized coordinate j. Its working null is a zero derivative in the evaluated function space. Under the null's centered Gaussian approximation, the reference is a weighted chi-square mixture with weights from `V^(1/2) Qj V^(1/2)`. Holm correction covers every candidate variable. Its centered reference distribution requires `Qj^(1/2) beta_lambda=0`. With strictly positive distances, a finite centered Gaussian expansion whose coordinate derivative vanishes on an open set is identically zero: analyticity extends the zero derivative globally, and sending that coordinate to infinity leaves a constant, which centering makes zero. On a finite evaluation grid, full-rank Qj likewise makes the working null equivalent to beta_lambda=0. Accordingly, the derivative rule is an empirical working-response screening procedure. Recognition and false inclusion against generating active/null labels are simulation performance measures. Holm's error-control theorem applies to valid working-null p values; variable-specific error control for generating nulls requires a separate null-compatible model or test. The appendix gives the proof and distinguishes these targets.

For each of four specified interaction pairs, mean squared grid contrast has the analogous covariance-spectrum mixture. Holm covers the four overall p values. Pointwise two-sided Wald tests have a second Holm family spanning all non-anchor grid points of all four pairs. Correct significant direction requires both rejection and the true sign. Additive comparison models have identically zero four-point operators and p=1.

Davies integration uses absolute accuracy 1e-8 and 100,000 terms, retried with 1,000,000 terms. Decision sensitivity is checked using the returned p value plus/minus the requested numerical tolerance and refined up to twice by factors of 100. Variable-tail fallback samples the mixture with `B=max(199999,ceiling(200*p/.05)-1)` and the plus-one exceedance probability. Available tests and numerical status are retained separately.

## 9. A–E roles and fixed configuration

The reviewed configuration is **A035: N=20,000, L=20,000, M=144, kappa=1e-4, lambda=1**. A's grid reference uses lambda=10. The source contains the published fixed file; the work-directory freeze file supplies execution settings.

| Module | Design and independent repetitions |
|---|---|
| A | 50 unique parameter settings, 100 repeats each; G1 KDR-CC. Truth-based weighted performance score and eligibility review. |
| B | 100 repeats; four KDR fits and two adjustment specifications per CLR method per repeat. Three methods displayed under G1–G4. |
| C | 8/16/40/80 null coordinates, 100 repeats each; fixed G1 KDR-CC; full function, interpretation and selection evaluation. |
| D | Mean absolute correlations .60/.85, pairwise floors .45/.75, tolerance .03; 100 repeats per level; B-style methods/scenarios. Original B supplies the reference. |
| E | Seven settings; 100 paired repeats each. Independent L=20,000 learning, 30,000 selection and N=20,000 final samples. E-specific data-derived nuisance dictionary for both strategies. |

E compares six (M,kappa) pairs crossed with lambda .3/1/3. Three-fold validation uses two folds for the full estimator and the third for unpenalized conditional negative log likelihood, preserving groups. An eligible candidate has three successful scores. The minimum weighted mean wins; deterministic ordering breaks ties. All-candidate failure invokes the recorded fixed fallback. Final fixed/adaptive fits share fresh estimation data; truth is used in evaluation after selection. Summaries average the selected procedure over all repeats.

## 10. Performance and reproducibility

Joint and observed-combination targets each report RMS Bias, RMSE, normalized errors, RMS empirical SD, RMS model SE, their ratio, pointwise coverage, mean interval width, absolute OR relative errors `100*abs(exp(error)-1)`, and correct/wrong significant directions. Bias is computed pointwise over independent repetitions before its squared aggregation. Evaluation coordinates and truths are verified identical within each summary cell. Replicate-level summaries define Monte Carlo units; RMS metrics use delete-one-repeat jackknife uncertainty.

Variable reporting includes individual active-variable recognition, mean true-positive rate and null false-positive rate. Every summary includes planned, completed and usable denominators. Scientific configurations, task seeds, matching, saved coefficients, covariances, learning objects, numerical statuses and source snapshots form the reproduction record. Independent PSOCK workers use deterministic task streams and single-threaded numerical libraries.

## 11. Implementation map

| Contract | Source |
|---|---|
| Calendar, history, event generation | `R/02_calendar.R`, `R/10_mechanism.R`, `config/mechanism.R` |
| Nuisance dictionaries | `R/11_nuisance_design.R`, `R/57_E_data.R` |
| Gaussian basis, geometry, whitening, learning and penalty | `R/20_kernel.R` |
| Cross-fitting and regularized moment solution | `R/30_estimation.R` |
| Exact curvature and full auxiliary sandwich | `R/31_covariance.R` |
| Function outputs and variable tests | `R/40_evaluation.R`, `R/41_response_grid.R` |
| Shapley allocation and interaction inference | `R/42_attribution.R`, `R/43_interaction_inference.R`, `R/45_interpretation.R` |
| A–E plans, data and execution | `R/50_design.R` through `R/66_E_execution.R` |
| Summary definitions and publication reports | `R/80_summary.R` through `R/92_E_publication.R` |
| Independent submission and compatibility | `R/95_submission.R`, `R/96_compatibility.R`, `05_submit.R` |
| Statistical identity and derivative checks | `tests/method_contract_checks.R` |

Primary sources for the underlying moment identity, matching, M-estimation, numerical optimization, quadratic-form tails and Shapley allocation are listed with full bibliographic records in the appendix.
