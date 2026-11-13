# Figure and table captions

Matched independent repetitions: 100

**Figure_G1_joint.png. Recovery of the joint exposure effect.** All eight exposures follow the archived joint quantile path. Curves show generating truth and means over the same matched independent repetitions. The blue light band averages pointwise 95% confidence limits separately; the darker band is the empirical central 95% range of package point estimates.

**Figure_G1_individual.png. Recovery at observed exposure combinations.** Each point represents one archived evaluation combination. Orange open circles and blue points give the archived and package repetition-mean log ORs; the diagonal is equality with generating truth.

**Figure_G1_single.png. Single-exposure effect recovery.** All other exposures remain at the archived marginal-median reference. Concentration values, references, contrasts and repetition IDs are identical in the two implementations. Blue bands have the same definitions as the joint-curve figure.

**Figure_G1_interactions.png. Interaction-effect recovery.** Four-point log-OR contrasts measure departures from additivity. Points compare truth with matched-cohort mean predictions. Blue vertical segments show central 95% repetition ranges. Anchor contrasts are excluded.

**Figure_G1_contribution_shares.png. Exposure contribution shares.** Baseline Shapley allocations use the same reference and evaluation combinations. Mean absolute allocations are normalized within the eight active exposures for each fit; bar heights are generating shares and averages of fitted shares.

**Figure_G1_interaction_detection.png. Interaction significance and direction recovery.** Overall tests use Holm correction over the four prespecified pairs. Correct significant directions use Holm over all eligible grid contrasts and require the true sign. Grid rates are averaged within each repetition before averaging across repetitions.

**Figure_G1_quantitative_OR.png. Quantitative odds-ratio recovery.** Crosses show true ORs. Orange and blue points are averages of OR estimates. Horizontal segments connect separately averaged pointwise interval endpoints. Joint quantile contrasts, single-exposure interquartile contrasts and temperature-versus-optimum contrasts follow the archived definitions.

**Figure_G1_OR_errors.png. Absolute relative OR error distributions.** Errors equal 100*abs(exp(estimated log OR minus true log OR)-1), pooled over matched repetitions and eligible evaluation points. Boxplots show quartiles and conventional whiskers; all numerical values enter the summary tables.

**Figure_G1_variable_identification.png. Identification of active exposures.** Derivative-based participation tests use Holm correction. Bars give the proportion of repetitions selecting each active exposure; error bars give Wilson Monte Carlo intervals.

**Model diagnostics.** The dedicated model_diagnostics directory compares final iteration counts, maximum gradients, regularized objective values and exact-curvature condition numbers.

**Tables.** CSV files preserve numerical columns for further analysis. Matching HTML files include a title, header and note and use three-line table rules. Availability, matched cohort, pointwise effects, errors, interval performance, contributions, quantitative ORs, participation and interaction detection are reported separately.
