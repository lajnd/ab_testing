function fetchFormAndCalculate() {
    try {
        // Fetch and parse input values from the form
        const controlVisitors = Number(document.getElementById('control-visitors').value);
        const controlConversions = Number(document.getElementById('control-conversions').value);
        const variantVisitors = Number(document.getElementById('variant-visitors').value);
        const variantConversions = Number(document.getElementById('variant-conversions').value);
        const confidenceLevelInput = Number(document.getElementById('confidence-level').value);
        const statisticalPowerInput = 80; // Default statistical power for now

        // Perform A/B test calculations
        const results = calculateAbTest(
            controlVisitors,
            controlConversions,
            variantVisitors,
            variantConversions,
            confidenceLevelInput,
            statisticalPowerInput
        );

        console.log(results);

        // Helper function to safely format values
        function formatValue(value) {
            return isFinite(value) && !isNaN(value) ? value.toFixed(2) : 'N/A';
        }

        // Multiply relevant values by 100 to convert to percentages
        const relativeDifferencePercent = results.lift;
        const relativeConfidenceIntervalPercent = results.relativeConfidenceInterval;
        const rightSidedIntervalPercent = results.rightSidedIntervalPP;
        const leftSidedIntervalPercent = results.leftSidedIntervalPP;        ;

        // Update the results in the DOM
        document.getElementById('conversion-rate-control').innerText = `${formatValue(results.conversionRateControl)}%`;
        document.getElementById('conversion-rate-variant').innerText = `${formatValue(results.conversionRateVariant)}%`;
        document.getElementById('lift').innerText = `${formatValue(results.lift)}%`;
        document.getElementById('absolute-difference').innerText = `${formatValue(results.differencePercentagePoints)}%`;
        document.getElementById('confidence-interval-diff').innerText =
            `${formatValue(results.confidenceIntervalDifferencePP[0])}% to ${formatValue(results.confidenceIntervalDifferencePP[1])}%`;
        document.getElementById('right-sided-interval').innerText =
            `${formatValue(rightSidedIntervalPercent[0])}% to ${rightSidedIntervalPercent[1] === Infinity ? '∞' : formatValue(rightSidedIntervalPercent[1])}%`;
        document.getElementById('left-sided-interval').innerText =
            `${leftSidedIntervalPercent[0] === -Infinity ? '-∞' : formatValue(leftSidedIntervalPercent[0])}% to ${formatValue(leftSidedIntervalPercent[1])}%`;
        document.getElementById('value-plus-minus-SE').innerText =
            `${formatValue(relativeDifferencePercent)} ± ${formatValue(results.valuePlusMinusSEPP[1])} %`;
        document.getElementById('p-value').innerText = results.pValue.toFixed(6);
        document.getElementById('z-score').innerText = formatValue(results.zScore);
        document.getElementById('significance').innerText = results.pValueOneSidedSignificance;
        document.getElementById('control-ci').innerText =
            `${formatValue(results.confidenceIntervalControl[0])}% to ${formatValue(results.confidenceIntervalControl[1])}%`;
        document.getElementById('variant-ci').innerText =
            `${formatValue(results.confidenceIntervalVariant[0])}% to ${formatValue(results.confidenceIntervalVariant[1])}%`;
        // Update Bayesian results
        document.getElementById('bayesian-variant-wins').innerText = `${formatValue(results.bayesianVariantWins)}%`;
        document.getElementById('bayesian-control-wins').innerText = `${formatValue(results.bayesianControlWins)}%`;
        document.getElementById('bayes-factor').innerText = formatValue(results.bayesFactorH1H0);
        // Update relative difference results
        const lowerCI = formatValue(relativeConfidenceIntervalPercent[0]);
        const upperCI = formatValue(relativeConfidenceIntervalPercent[1]);
        document.getElementById('relative-confidence-interval').innerText = `${lowerCI}% to ${upperCI}%`;
        // TODO, there is no good way to calculate the p-value right now
        // document.getElementById('relative-confidence-interval-right').innerText =
        //     `${formatValue(rightSidedIntervalPercent[0])}% to ${rightSidedIntervalPercent[1] === Infinity ? '∞' : formatValue(rightSidedIntervalPercent[1])}%`;
        // document.getElementById('relative-confidence-interval-left').innerText =
        //     `${leftSidedIntervalPercent[0] === -Infinity ? '-∞' : formatValue(leftSidedIntervalPercent[0])}% to ${formatValue(leftSidedIntervalPercent[1])}%`;
        // // fill in relative-difference-plus-minus-SE
        // document.getElementById('relative-difference-plus-minus-SE').innerText = results.relativeDifferencePlusMinusSE;
        // document.getElementById('relative-p-value').innerText = formatValue(results.pValueIterative);
        // document.getElementById('relative-z-score').innerText = formatValue(results.relativeZScore);
    } catch (error) {
        console.error(error);
        alert(`An error occurred: ${error.message}`);
    }
}

/**
 * Calculates various statistical metrics for A/B testing.
 * @param {number} controlVisitors - Number of visitors in the control group.
 * @param {number} controlConversions - Number of conversions in the control group.
 * @param {number} variantVisitors - Number of visitors in the variant group.
 * @param {number} variantConversions - Number of conversions in the variant group.
 * @param {number} confidenceLevelInput - Desired confidence level (e.g., 95 for 95% confidence).
 * @param {number} statisticalPowerInput - Desired statistical power (e.g., 80 for 80% power).
 * @returns {Object} - An object containing calculated metrics.
 */
function calculateAbTest(
    controlVisitors,
    controlConversions,
    variantVisitors,
    variantConversions,
    confidenceLevelInput,
    statisticalPowerInput
) {
    const confidenceLevel = confidenceLevelInput / 100;
    const statisticalPower = statisticalPowerInput / 100;

    // Calculate conversion rates
    const conversionRateControl = controlConversions / controlVisitors;
    const conversionRateVariant = variantConversions / variantVisitors;

    // Calculate p-value and z-score
    const pValueResult = calculatePValue(
        controlConversions,
        controlVisitors,
        variantConversions,
        variantVisitors
    );

    const zScore = pValueResult.zScore;
    const pValue = pValueResult.pValue;
    const pValueOneSidedSignificance = pValue < (1 - confidenceLevel) ? 'Significant' : 'Not Significant';

    // Calculate lift and difference in percentage points
    const absoluteDifference = conversionRateVariant - conversionRateControl;
    const lift = (absoluteDifference / conversionRateControl) * 100;
    const differencePercentagePoints = absoluteDifference * 100;

    // Standard error for the difference
    const standardErrorDifference = Math.sqrt(
        (conversionRateControl * (1 - conversionRateControl)) / controlVisitors +
        (conversionRateVariant * (1 - conversionRateVariant)) / variantVisitors
    );

    // Z-scores for confidence intervals
    const zAlphaTwoSided = jStat.normal.inv(1 - (1 - confidenceLevel) / 2, 0, 1);
    const zAlphaOneSided = jStat.normal.inv(confidenceLevel, 0, 1);

    // Confidence intervals for the difference
    const lowerCIDiff = absoluteDifference - zAlphaTwoSided * standardErrorDifference;
    const upperCIDiff = absoluteDifference + zAlphaTwoSided * standardErrorDifference;
    const confidenceIntervalDifferencePP = [lowerCIDiff * 100, upperCIDiff * 100];

    // One-sided intervals
    const lowerCIRightSidedPP = (absoluteDifference - zAlphaOneSided * standardErrorDifference) * 100;
    const upperCILeftSidedPP = (absoluteDifference + zAlphaOneSided * standardErrorDifference) * 100;

    // Value ± SE
    const SEPP = zAlphaTwoSided * standardErrorDifference * 100;
    const valuePlusMinusSEPP = [differencePercentagePoints, SEPP];

    const relativeABResult = calculateRelativeDifference(
        controlVisitors,
        controlConversions,
        variantVisitors,
        variantConversions,
        confidenceLevel,
        statisticalPower
    );

    const relativePValue = relativeABResult.pValue;
    const relativeZScore = relativeABResult.zScore;
    const relativeConfidenceInterval = relativeABResult.relativeConfidenceInterval;

    const CV1 = Math.sqrt((conversionRateControl * (1 - conversionRateControl)) / controlVisitors); // Coefficient of variation for control group
    const CV2 = Math.sqrt((conversionRateVariant * (1 - conversionRateVariant)) / variantVisitors); // Coefficient of variation for variant group
    const relativeDifference = (conversionRateVariant - conversionRateControl) / conversionRateControl; // Relative difference

    // Now we call the calculatePValueIteratively function
    const pValueIterative = calculatePValueIteratively(relativeDifference, CV1, CV2);


    // Confidence intervals for control and variant
    const zAlpha = jStat.normal.inv(confidenceLevel, 0, 1);
    const controlCI = calculateConfidenceInterval(conversionRateControl, controlVisitors, zAlpha);
    const variantCI = calculateConfidenceInterval(conversionRateVariant, variantVisitors, zAlpha);



    // Bayesian probability
    const bayesianResults = calculateBayesianProbability(
        controlConversions,
        controlVisitors,
        variantConversions,
        variantVisitors
    );

    return {
        conversionRateControl: conversionRateControl * 100,
        conversionRateVariant: conversionRateVariant * 100,
        lift,
        differencePercentagePoints,
        confidenceIntervalDifferencePP,
        rightSidedIntervalPP: [lowerCIRightSidedPP, Infinity],
        leftSidedIntervalPP: [-Infinity, upperCILeftSidedPP],
        valuePlusMinusSEPP,
        pValue,
        zScore,
        pValueOneSidedSignificance,
        relativeConfidenceInterval,
        relativeConfidenceIntervalRightSided: relativeABResult.relativeConfidenceIntervalRightSided,
        relativeConfidenceIntervalLeftSided: relativeABResult.relativeConfidenceIntervalLeftSided,
        relativeDifferencePlusMinusSE: relativeABResult.relativeDifferencePlusMinusSE,
        pValueIterative,
        relativeZScore,
        confidenceIntervalControl: [controlCI[0] * 100, controlCI[1] * 100],
        confidenceIntervalVariant: [variantCI[0] * 100, variantCI[1] * 100],
        bayesianVariantWins: bayesianResults.probabilityVariantWins * 100,
        bayesianControlWins: bayesianResults.probabilityControlWins * 100,
        bayesFactorH1H0: bayesianResults.bayesFactorH1H0
    };
}

/**
 * Calculates the confidence interval for a given conversion rate.
 * @param {number} conversionRate - The conversion rate (proportion between 0 and 1).
 * @param {number} visitors - The number of visitors.
 * @param {number} zAlpha - The z-score corresponding to the desired confidence level.
 * @returns {number[]} - An array containing the lower and upper bounds of the confidence interval.
 */
function calculateConfidenceInterval(conversionRate, visitors, zAlpha) {
    const marginOfError = zAlpha * Math.sqrt((conversionRate * (1 - conversionRate)) / visitors);
    return [conversionRate - marginOfError, conversionRate + marginOfError];
}

/**
 * Calculates the relative Minimum Detectable Effect (MDE) as a percentage.
 * @param {number} controlSampleSize - Sample size of the control group.
 * @param {number} variantSampleSize - Sample size of the variant group.
 * @param {number} baselineConversionRate - Baseline conversion rate (proportion between 0 and 1).
 * @param {number} confidenceLevel - Desired confidence level (proportion between 0 and 1).
 * @param {number} statisticalPower - Desired statistical power (proportion between 0 and 1).
 * @returns {number} - Relative MDE as a percentage.
 */
function calculateRelativeMDE(controlSampleSize, variantSampleSize, baselineConversionRate, confidenceLevel, statisticalPower) {
    const zAlpha = jStat.normal.inv(confidenceLevel, 0, 1);
    const zBeta = jStat.normal.inv(statisticalPower, 0, 1);

    const standardError = Math.sqrt(
        (baselineConversionRate * (1 - baselineConversionRate) / controlSampleSize) +
        (baselineConversionRate * (1 - baselineConversionRate) / variantSampleSize)
    );

    const absoluteMDE = (zAlpha + zBeta) * standardError;
    return (absoluteMDE / baselineConversionRate) * 100;
}

/**
 * Calculates the p-value and z-score for an A/B test.
 * @param {number} controlConversions - Number of conversions in the control group.
 * @param {number} controlVisitors - Number of visitors in the control group.
 * @param {number} variantConversions - Number of conversions in the variant group.
 * @param {number} variantVisitors - Number of visitors in the variant group.
 * @param {boolean} isOneSided - Whether the test is one-sided (default is true).
 * @returns {Object} - An object containing the z-score, p-value, and isSignificant (boolean).
 */
function calculatePValue(
    controlConversions,
    controlVisitors,
    variantConversions,
    variantVisitors,
    isOneSided = true
) {
    const conversionRateControl = controlConversions / controlVisitors;
    const conversionRateVariant = variantConversions / variantVisitors;
    const pooledConversionRate = (controlConversions + variantConversions) / (controlVisitors + variantVisitors);

    const pooledSE = Math.sqrt(
        pooledConversionRate * (1 - pooledConversionRate) * (1 / controlVisitors + 1 / variantVisitors)
    );

    if (pooledSE === 0) {
        throw new Error("Pooled standard error is zero. Check your input values.");
    }

    const zScore = (conversionRateVariant - conversionRateControl) / pooledSE;

    let pValue;
    if (isOneSided) {
        pValue = 1 - jStat.normal.cdf(zScore, 0, 1);
    } else {
        pValue = 2 * (1 - jStat.normal.cdf(Math.abs(zScore), 0, 1));
    }

    pValue = Math.min(Math.max(pValue, 0), 1);

    return {
        zScore,
        pValue,
        isSignificant: pValue < 0.05
    };
}

/**
 * Calculates the required sample size per group for an A/B test.
 * @param {number} p1 - The baseline conversion rate (proportion between 0 and 1).
 * @param {number} liftPercentage - The expected percentage increase in conversion rate.
 * @param {number} confidenceLevel - Desired confidence level (e.g., 95 for 95% confidence).
 * @param {number} power - Desired statistical power (e.g., 80 for 80% power).
 * @param {number} numVariants - Number of variants being tested (default is 2).
 * @param {boolean} isOneSided - Whether the test is one-sided (default is true).
 * @returns {number} - The required sample size per group.
 */
function calculateSampleSize(p1, liftPercentage, confidenceLevel, power, numVariants = 2, isOneSided = true) {
    let p2 = p1 * (1 + liftPercentage / 100);
    p2 = Math.min(p2, 0.9999);

    const alpha = 1 - confidenceLevel / 100;
    let zAlpha = isOneSided ? jStat.normal.inv(1 - alpha, 0, 1) : jStat.normal.inv(1 - alpha / 2, 0, 1);
    const zBeta = jStat.normal.inv(power / 100, 0, 1);

    if (numVariants > 2) {
        const alphaAdjusted = alpha / (numVariants - 1);
        zAlpha = jStat.normal.inv(1 - alphaAdjusted, 0, 1);
    }

    const pAvg = (p1 + p2) / 2;
    const sigmaPooled = Math.sqrt(2 * pAvg * (1 - pAvg));
    const standardError = Math.sqrt(p1 * (1 - p1) + p2 * (1 - p2));
    const deltaP = p2 - p1;

    const numerator = Math.pow(zAlpha * sigmaPooled + zBeta * standardError, 2);
    const denominator = Math.pow(deltaP, 2);

    return Math.ceil(numerator / denominator);
}

/**
 * Calculates the relative difference (lift) between two groups along with confidence intervals, Z-score, and P-value.
 *
 * @param {number} controlVisitors - Number of visitors in the control group (n1)
 * @param {number} controlConversions - Number of conversions in the control group (x1)
 * @param {number} variantVisitors - Number of visitors in the variant group (n2)
 * @param {number} variantConversions - Number of conversions in the variant group (x2)
 * @param {number} confidenceLevelInput - Confidence level (e.g., 0.95 for 95%)
 * @returns {Object} - An object containing the relative difference, confidence intervals, SE, P-value, and Z-score
 */
function calculateRelativeDifference(controlVisitors, controlConversions, variantVisitors, variantConversions, confidenceLevelInput) {
    // Input Validation
    if (controlVisitors === 0 || variantVisitors === 0) {
        throw new Error("Number of visitors must be greater than zero for both groups.");
    }
    if (controlConversions > controlVisitors || variantConversions > variantVisitors) {
        throw new Error("Number of conversions cannot exceed number of visitors in any group.");
    }

    // Ensure confidence level is between 0 and 1
    if (confidenceLevelInput <= 0 || confidenceLevelInput >= 1) {
        throw new Error("Confidence level must be between 0 and 1 (e.g., 0.95 for 95%).");
    }

    // 1. Calculate Proportions (Conversion Rates)
    const p1 = controlConversions / controlVisitors; // Control conversion rate
    const p2 = variantConversions / variantVisitors; // Variant conversion rate

    // Safeguard small or zero conversion rates to avoid invalid calculations
    const minValue = 1e-9; // Small value to prevent division by zero or invalid values
    const p1_safe = Math.max(p1, minValue);
    const p2_safe = Math.max(p2, minValue);

    // 2. Calculate Relative Difference (Lift)
    const lift = (p2_safe / p1_safe) - 1;

    // 3. Calculate Coefficients of Variation (CV1 & CV2)
    const CV1 = Math.sqrt((1 - p1_safe) / (p1_safe * controlVisitors));
    const CV2 = Math.sqrt((1 - p2_safe) / (p2_safe * variantVisitors));

    // 4. Correct Z-score Calculation for the Given Confidence Level
    const alpha = 1 - confidenceLevelInput; // confidenceLevelInput is already between 0 and 1
    const Z_twoSided = jStat.normal.inv(1 - alpha / 2, 0, 1); // Should give ~1.96 for 95% confidence

    // 5. Compute the Two-Sided Confidence Interval for Relative Difference
    const sqrtTermRaw = CV1 ** 2 + CV2 ** 2 - (Z_twoSided ** 2) * (CV1 ** 2) * (CV2 ** 2);
    let sqrtTerm = NaN;
    if (sqrtTermRaw >= 0) {
        sqrtTerm = Math.sqrt(sqrtTermRaw);
    }

    const numeratorLower = 1 - (Z_twoSided * CV1 ** 2);
    const numeratorUpper = 1 - (Z_twoSided * CV1 ** 2);

    const factorLower = 1 - (Z_twoSided * sqrtTerm);
    const factorUpper = 1 + (Z_twoSided * sqrtTerm);

    const relativeCILower = numeratorLower !== 0 && !isNaN(factorLower) ? ((1 + lift) * factorLower) / numeratorLower - 1 : NaN;
    const relativeCIUpper = numeratorUpper !== 0 && !isNaN(factorUpper) ? ((1 + lift) * factorUpper) / numeratorUpper - 1 : NaN;


    // 6. Compute One-Sided Confidence Intervals for Relative Difference

    // Right-Sided Interval Calculation
    const Z_oneSided = jStat.normal.inv(confidenceLevelInput, 0, 1); // e.g., 1.645 for 95% one-sided
    const sqrtTermOneSided = Math.sqrt(CV1 ** 2 + CV2 ** 2 - (Z_oneSided ** 2) * (CV1 ** 2) * (CV2 ** 2));
    
    const rightSidedFactor = 1 + Z_oneSided * sqrtTermOneSided;
    const lowerCIRightSided = (1 + lift) * rightSidedFactor - 1;  // Lower bound for right-sided CI

    // Left-Sided Interval Calculation
    const leftSidedFactor = 1 - Z_oneSided * sqrtTermOneSided;
    const upperCILeftSided = (1 + lift) * leftSidedFactor - 1;  // Upper bound for left-sided CI

    // Right-Sided Interval should be [lowerCI, +∞]
    const rightSidedInterval = [lowerCIRightSided, Infinity];

    // Left-Sided Interval should be [-∞, upperCI]
    const leftSidedInterval = [-Infinity, upperCILeftSided];

    // 7. Calculate Standard Error (SE)
    const SE_relative = Math.sqrt(CV1 ** 2 + CV2 ** 2);

    // 8. Calculate Z-score for Hypothesis Testing (H0: B ≤ A vs. Ha: B > A)
    const zScoreFromCI = lift / SE_relative;

    // 9. Now we call the p-value calculation function and pass the confidence intervals
    const pValueIterative = calculatePValueIteratively(relativeCILower, relativeCIUpper, Z_oneSided);

    // 10. Format Value ± SE (Convert lift and SE to percentages)
    const SEPP = (SE_relative * 100).toFixed(3); // Rounded to 3 decimal places

    // 11. Assemble the Results
    return {
        relativeDifference: lift * 100,  // Convert to percentage
        relativeConfidenceInterval: [relativeCILower * 100, relativeCIUpper * 100],  // Convert to percentage
        relativeConfidenceIntervalRightSided: rightSidedInterval.map(val => val === Infinity ? val : val * 100),  // Right-sided interval
        relativeConfidenceIntervalLeftSided: leftSidedInterval.map(val => val === -Infinity ? val : val * 100),  // Left-sided interval
        relativeDifferencePlusMinusSE: `${(lift * 100).toFixed(2)} ± ${SEPP}%`,
        pValue: pValueIterative,  // Use the iterative method to calculate the p-value
        zScore: zScoreFromCI
    };
}

/**
 * Function to calculate a p-value by iterating over confidence intervals.
 * @param {number} relativeCILower - The lower bound of the relative confidence interval.
 * @param {number} relativeCIUpper - The upper bound of the relative confidence interval.
 * @returns {number} - Approximate p-value for the relative difference.
 */
function calculatePValueIteratively(relativeCILower, relativeCIUpper, Z) {
    let iteration = 0;
    const step = 0.001;  // Increment Z-value gradually
    const max_iterations = 10000;

    // Iterate over Z-values to find the one where the CI excludes the null hypothesis
    while (iteration < max_iterations) {
        // Check if the confidence interval excludes 0 (null hypothesis)
        if (relativeCILower > 0 || relativeCIUpper < 0) {
            // Found the Z-value that excludes the null hypothesis
            break;
        }
        
        // Adjust Z-value and increment iteration count
        Z += step;
        iteration++;
    }

    // Calculate the p-value corresponding to the Z-value found
    const p_value = 1 - cumulativeStandardNormal(Z); // One-tailed p-value

    // Return the approximate p-value
    return p_value;
}


/**
 * Approximation for the error function (erf) in JavaScript
 * @param {number} x - The value to calculate erf for
 * @returns {number} - The approximated value of erf(x)
 */
function erf(x) {
    // Constants for approximation
    const a1 =  0.254829592;
    const a2 = -0.284496736;
    const a3 =  1.421413741;
    const a4 = -1.453152027;
    const a5 =  1.061405429;
    const p  =  0.3275911;

    // Save the sign of x
    const sign = (x >= 0) ? 1 : -1;
    x = Math.abs(x);

    // A&S formula 7.1.26 approximation for erf
    const t = 1.0 / (1.0 + p * x);
    const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);

    return sign * y;
}
/**
 * Cumulative standard normal distribution function to convert Z to p-value.
 * @param {number} Z - The Z-score.
 * @returns {number} - The cumulative probability corresponding to the Z-score.
 */
function cumulativeStandardNormal(Z) {
    // Approximation of the cumulative distribution function for a standard normal distribution
    return (1.0 + erf(Z / Math.sqrt(2))) / 2;
}

/**
 * Calculates the confidence interval bounds for a relative difference.
 * @param {number} RelDiff - The relative difference between two groups.
 * @param {number} Z - Z-value from a standard normal distribution (related to confidence level).
 * @param {number} CV1 - Coefficient of variation for the first group.
 * @param {number} CV2 - Coefficient of variation for the second group.
 * @returns {Object} - An object containing the lower and upper confidence interval bounds.
 */
function computeRelativeCIBounds(RelDiff, Z, CV1, CV2) {
    // Calculate the inner square root part of the equation
    const sqrtTerm = Math.sqrt(CV1**2 + CV2**2 - Z**2 * CV1**2 * CV2**2);
    
    // Calculate the denominator part
    const denominator = 1 - Z * CV1**2;
    
    // Calculate the positive and negative bounds of the CI
    const positiveBound = (1 + Z * sqrtTerm) / denominator;
    const negativeBound = (1 - Z * sqrtTerm) / denominator;
    
    // Calculate the final CI bounds for both positive and negative cases
    const CI_positive = (RelDiff + 1) * positiveBound - 1;
    const CI_negative = (RelDiff + 1) * negativeBound - 1;
    
    return {
        CI_negative,
        CI_positive
    };
}

/**
 * Calculates the Bayesian probability that the variant is better than the control.
 * @param {number} controlConversions - Number of conversions in the control group.
 * @param {number} controlVisitors - Number of visitors in the control group.
 * @param {number} variantConversions - Number of conversions in the variant group.
 * @param {number} variantVisitors - Number of visitors in the variant group.
 * @param {number} simulations - Number of Monte Carlo simulations to run (default is 100,000).
 * @returns {Object} - An object containing probabilities and Bayes Factor.
 */
function calculateBayesianProbability(
    controlConversions,
    controlVisitors,
    variantConversions,
    variantVisitors,
    simulations = 100000
) {
    const alphaPrior = 1;
    const betaPrior = 1;

    const controlPosteriorAlpha = alphaPrior + controlConversions;
    const controlPosteriorBeta = betaPrior + (controlVisitors - controlConversions);
    const variantPosteriorAlpha = alphaPrior + variantConversions;
    const variantPosteriorBeta = betaPrior + (variantVisitors - variantConversions);

    let variantWins = 0;

    for (let i = 0; i < simulations; i++) {
        const controlSample = jStat.beta.sample(controlPosteriorAlpha, controlPosteriorBeta);
        const variantSample = jStat.beta.sample(variantPosteriorAlpha, variantPosteriorBeta);

        if (variantSample > controlSample) {
            variantWins++;
        }
    }

    const probabilityVariantWins = variantWins / simulations;
    const probabilityControlWins = 1 - probabilityVariantWins;

    const lnMlH1 =
        jStat.betaln(controlConversions + alphaPrior, controlVisitors - controlConversions + betaPrior) -
        jStat.betaln(alphaPrior, betaPrior) +
        jStat.betaln(variantConversions + alphaPrior, variantVisitors - variantConversions + betaPrior) -
        jStat.betaln(alphaPrior, betaPrior);

    const totalConversions = controlConversions + variantConversions;
    const totalVisitors = controlVisitors + variantVisitors;
    const lnMlH0 =
        jStat.betaln(totalConversions + alphaPrior, totalVisitors - totalConversions + betaPrior) -
        jStat.betaln(alphaPrior, betaPrior);

    const lnBF10 = lnMlH1 - lnMlH0;
    let BF10 = Math.exp(lnBF10);

    if (!isFinite(BF10)) {
        BF10 = lnBF10 > 0 ? Infinity : 0;
    }

    return {
        probabilityVariantWins,
        probabilityControlWins,
        bayesFactorH1H0: BF10
    };
}