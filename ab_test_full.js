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
        const statisticalPowerInput = Number(document.getElementById('statistical-power').value);

        // Validate input values
        if (
            [controlVisitors, controlConversions, variantVisitors, variantConversions, confidenceLevelInput, statisticalPowerInput].some(isNaN)
        ) {
            throw new Error('Please enter valid numerical values.');
        }

        // Perform A/B test calculations
        const results = calculateAbTest(
            controlVisitors,
            controlConversions,
            variantVisitors,
            variantConversions,
            confidenceLevelInput,
            statisticalPowerInput
        );

        // Update the results
        document.getElementById('conversion-rate-control').innerText = `${results.conversionRateControl.toFixed(2)}%`;
        document.getElementById('conversion-rate-variant').innerText = `${results.conversionRateVariant.toFixed(2)}%`;
        document.getElementById('lift').innerText = `${results.lift.toFixed(2)}%`;
        document.getElementById('difference-pp').innerText = `${results.differencePercentagePoints.toFixed(2)} pp`;
        document.getElementById('confidence-interval-diff').innerText =
            `${results.confidenceIntervalDifferencePP[0].toFixed(2)}% to ${results.confidenceIntervalDifferencePP[1].toFixed(2)}%`;
        document.getElementById('right-sided-interval').innerText =
            `${results.rightSidedIntervalPP[0].toFixed(2)}% to ${results.rightSidedIntervalPP[1] === Infinity ? '∞' : `${results.rightSidedIntervalPP[1].toFixed(2)}%`}`;
        document.getElementById('left-sided-interval').innerText =
            `${results.leftSidedIntervalPP[0] === -Infinity ? '-∞' : results.leftSidedIntervalPP[0].toFixed(2)}% to ${results.leftSidedIntervalPP[1].toFixed(2)}%`;
        document.getElementById('value-plus-minus-SE').innerText =
            `${results.valuePlusMinus95SEPP[0].toFixed(2)} ± ${results.valuePlusMinus95SEPP[1].toFixed(2)} pp`;
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

    // Validate conversion rates
    if (
        [conversionRateControl, conversionRateVariant].some(
            rate => rate <= 0 || rate >= 1
        )
    ) {
        throw new Error('Conversion rates must be between 0% and 100% (exclusive).');
    }

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

    // Value ± 95% SE
    const SE95PP = zAlphaTwoSided * standardErrorDifference * 100;
    const valuePlusMinus95SEPP = [differencePercentagePoints, SE95PP];

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
        valuePlusMinus95SEPP,
        pValue,
        zScore,
        pValueOneSidedSignificance,
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