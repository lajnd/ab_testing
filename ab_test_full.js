/**
 * Fetches input values from the form, performs A/B test calculations, and updates the results in the DOM.
 */
function fetchFormAndCalculate() {
    try {
        // Fetch and parse input values from the form
        const controlVisitors = Number(document.getElementById('control-visitors').value);
        const controlConversions = Number(document.getElementById('control-conversions').value);
        const variantVisitors = Number(document.getElementById('variant-visitors').value);
        const variantConversions = Number(document.getElementById('variant-conversions').value);
        const confidenceLevelInput = Number(document.getElementById('confidence-level').value);
        // const statisticalPowerInput = Number(document.getElementById('statistical-power').value);
        const statisticalPowerInput = 80; // Default statistical for now

        // Perform A/B test calculations
        const results = calculateAbTest(
            controlVisitors,
            controlConversions,
            variantVisitors,
            variantConversions,
            confidenceLevelInput,
            statisticalPowerInput
        );

        // log all values to console
        console.log(results);

        // Update the results
        document.getElementById('conversion-rate-control').innerText = `${results.conversionRateControl.toFixed(2)}%`;
        document.getElementById('conversion-rate-variant').innerText = `${results.conversionRateVariant.toFixed(2)}%`;
        document.getElementById('lift').innerText = `${results.lift.toFixed(2)}%`;
        document.getElementById('absolute-difference').innerText = `${results.differencePercentagePoints.toFixed(2)}%`;
        document.getElementById('confidence-interval-diff').innerText =
            `${results.confidenceIntervalDifferencePP[0].toFixed(2)}% to ${results.confidenceIntervalDifferencePP[1].toFixed(2)}%`;
        document.getElementById('right-sided-interval').innerText =
            `${results.rightSidedIntervalPP[0].toFixed(2)}% to ${results.rightSidedIntervalPP[1] === Infinity ? '∞' : `${results.rightSidedIntervalPP[1].toFixed(2)}%`}`;
        document.getElementById('left-sided-interval').innerText =
            `${results.leftSidedIntervalPP[0] === -Infinity ? '-∞' : results.leftSidedIntervalPP[0].toFixed(2)}% to ${results.leftSidedIntervalPP[1].toFixed(2)}%`;
        document.getElementById('value-plus-minus-SE').innerText =
            `${results.valuePlusMinusSEPP[0].toFixed(2)} ± ${results.valuePlusMinusSEPP[1].toFixed(2)} %`;
        document.getElementById('p-value').innerText = results.pValue.toFixed(4);
        document.getElementById('z-score').innerText = results.zScore.toFixed(2);
        document.getElementById('significance').innerText = results.pValueOneSidedSignificance;
        document.getElementById('sample-size').innerText = results.sampleSizePerGroup;
        document.getElementById('control-ci').innerText =
            `${results.confidenceIntervalControl[0].toFixed(2)}% to ${results.confidenceIntervalControl[1].toFixed(2)}%`;
        document.getElementById('variant-ci').innerText =
            `${results.confidenceIntervalVariant[0].toFixed(2)}% to ${results.confidenceIntervalVariant[1].toFixed(2)}%`;
        document.getElementById('relative-mde').innerText = `${results.relativeMDE.toFixed(2)}%`;
        document.getElementById('bayesian-variant-wins').innerText = `${results.bayesianVariantWins.toFixed(2)}%`;
        document.getElementById('bayesian-control-wins').innerText = `${results.bayesianControlWins.toFixed(2)}%`;
        document.getElementById('bayes-factor').innerText = results.bayesFactorH1H0.toFixed(2);

/*         // Update relative difference results
        console.log(results.relativeConfidenceInterval);
        const lowerCI = results.relativeConfidenceInterval[0].toFixed(2);
        const upperCI = results.relativeConfidenceInterval[1].toFixed(2);
        console.log(lowerCI, upperCI);
        document.getElementById('relative-confidence-interval').innerText = `${lowerCI}% to ${upperCI}%`;

        document.getElementById('relative-confidence-interval-right').innerText =
            `${results.relativeConfidenceIntervalRightSided[0].toFixed(2)}% to ${results.relativeConfidenceIntervalRightSided[1] === Infinity ? '∞' : `${results.relativeConfidenceIntervalRightSided[1].toFixed(2)}%`}`;
        document.getElementById('relative-confidence-interval-left').innerText =
            `${results.relativeConfidenceIntervalLeftSided[0] === -Infinity ? '-∞' : results.relativeConfidenceIntervalLeftSided[0].toFixed(2)}% to ${results.relativeConfidenceIntervalLeftSided[1].toFixed(2)}%`;
        
            console.log(results.valuePlusMinusSEPP[1].toFixed(2),)
        //document.getElementById('value-plus-minus-SE').innerText =
        //    `${results.valuePlusMinusSEPP[0].toFixed(2)} ± ${results.valuePlusMinusSEPP[1].toFixed(2)} %`;
        console.log(results.pValue.toFixed(4));
        console.log(results.zScore.toFixed(2));
        document.getElementById('relative-p-value').innerText = results.relativePValue.toFixed(4);
        document.getElementById('relative-z-score').innerText = results.relativeZScore.toFixed(2); */
    } catch (error) {
        console.error(error);
        alert(`An error occurred: ${error.message}`);
    };
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
        confidenceLevelInput,
        statisticalPowerInput
    );

    const relativePValue = relativeABResult.pValue;
    const relativeZScore = relativeABResult.zScore;
    const relativeConfidenceInterval = relativeABResult.relativeConfidenceInterval;

    // Sample size calculation
    const sampleSizePerGroup = calculateSampleSize(
        conversionRateControl,
        lift,
        confidenceLevelInput,
        statisticalPowerInput
    );

    // Confidence intervals for control and variant
    const zAlpha = jStat.normal.inv(confidenceLevel, 0, 1);
    const controlCI = calculateConfidenceInterval(conversionRateControl, controlVisitors, zAlpha);
    const variantCI = calculateConfidenceInterval(conversionRateVariant, variantVisitors, zAlpha);

    // Relative MDE
    const relativeMDE = calculateRelativeMDE(
        controlVisitors,
        variantVisitors,
        conversionRateControl,
        confidenceLevel,
        statisticalPower
    );

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
        relativePValue,
        relativeZScore,
        sampleSizePerGroup,
        confidenceIntervalControl: [controlCI[0] * 100, controlCI[1] * 100],
        confidenceIntervalVariant: [variantCI[0] * 100, variantCI[1] * 100],
        relativeMDE,
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
 * Calculates various statistical metrics for A/B testing.
 * @param {number} controlVisitors - Number of visitors in the control group.
 * @param {number} controlConversions - Number of conversions in the control group.
 * @param {number} variantVisitors - Number of visitors in the variant group.
 * @param {number} variantConversions - Number of conversions in the variant group.
 * @param {number} confidenceLevelInput - Desired confidence level (e.g., 95 for 95% confidence).
 * @param {number} statisticalPowerInput - Desired statistical power (e.g., 80 for 80% power).
 * @returns {Object} - An object containing calculated metrics.
 */
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

    // 1. Calculate Proportions
    const p1 = controlConversions / controlVisitors; // Control conversion rate
    const p2 = variantConversions / variantVisitors; // Variant conversion rate

    // 2. Calculate Relative Difference (Lift)
    const lift = (p2 / p1) - 1;

    // 3. Calculate Coefficients of Variation (CV1 & CV2)
    const CV1 = Math.sqrt((1 - p1) / (p1 * controlVisitors));
    const CV2 = Math.sqrt((1 - p2) / (p2 * variantVisitors));

    // 4. Determine Z-value for the given confidence level
    // For two-sided confidence interval
    const alpha = 1 - confidenceLevelInput;
    const Z = jStat.normal.inv(1 - alpha / 2, 0, 1); // e.g., 1.96 for 95% confidence

    // 5. Compute the Confidence Interval for Relative Difference
    const sqrtTerm = Math.sqrt(CV1 ** 2 + CV2 ** 2 - (Z ** 2) * (CV1 ** 2) * (CV2 ** 2));
    const numeratorLower = 1 - (Z * CV1 ** 2);
    const numeratorUpper = 1 - (Z * CV1 ** 2);

    const factorLower = 1 - (Z * sqrtTerm);
    const factorUpper = 1 + (Z * sqrtTerm);

    const relativeCILower = ((1 + lift) * factorLower) / numeratorLower - 1;
    const relativeCIUpper = ((1 + lift) * factorUpper) / numeratorUpper - 1;

    const relativeCI = [relativeCILower, relativeCIUpper];

    // 6. Compute One-Sided Confidence Intervals
    const Z_oneSided = jStat.normal.inv(confidenceLevelInput, 0, 1); // e.g., 1.645 for 95% one-sided

    // Right-Sided Confidence Interval
    const sqrtTermOneSided = Math.sqrt(CV1 ** 2 + CV2 ** 2 - (Z_oneSided ** 2) * (CV1 ** 2) * (CV2 ** 2));
    const factorRight = 1 + (Z_oneSided * sqrtTermOneSided);
    const lowerCIRightSided = ((1 + lift) * factorRight) / (1 - (Z_oneSided * CV1 ** 2)) - 1;

    // Left-Sided Confidence Interval
    const factorLeft = 1 - (Z_oneSided * sqrtTermOneSided);
    const upperCILeftSided = ((1 + lift) * factorLeft) / (1 - (Z_oneSided * CV1 ** 2)) - 1;

    // 7. Calculate Standard Error (SE)
    // Rearranging the two-sided CI formula to solve for SE is not straightforward.
    // Instead, we approximate SE using the formula for relative difference:
    const SE_relative = Math.sqrt(CV1 ** 2 + CV2 ** 2);

    // 8. Calculate Z-score for Hypothesis Testing (H0: B ≤ A vs. Ha: B > A)
    const zScoreFromCI = lift / SE_relative;

    // 9. Calculate P-value for One-Tailed Test
    const pValueFromCI = 1 - jStat.normal.cdf(zScoreFromCI, 0, 1);

    // 10. Format Value ± SE
    const SEPP = SE_relative.toFixed(3); // Rounded to 3 decimal places

    // 11. Assemble the Results
    return {
        relativeDifference: lift,
        relativeConfidenceInterval: relativeCI,
        relativeConfidenceIntervalRightSided: [lowerCIRightSided, Infinity],
        relativeConfidenceIntervalLeftSided: [-Infinity, upperCILeftSided],
        relativeDifferencePlusMinusSE: `${lift.toFixed(4)} ±${SEPP}`,
        pValue: pValueFromCI,
        zScore: zScoreFromCI
    };
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