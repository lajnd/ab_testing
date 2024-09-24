
function fetchFormAndCalculate() {
    // Fetch input values from the form
    const controlVisitors = parseInt(document.getElementById('control-visitors').value);
    const controlConversions = parseInt(document.getElementById('control-conversions').value);
    const variantVisitors = parseInt(document.getElementById('variant-visitors').value);
    const variantConversions = parseInt(document.getElementById('variant-conversions').value);
    const confidenceLevelInput = parseFloat(document.getElementById('confidence-level').value);
    const statisticalPowerInput = parseFloat(document.getElementById('statistical-power').value);
    const desiredEffectSizePercentage = parseFloat(document.getElementById('desired-effect-size').value);

    // Call calculateAbTest with the fetched values
    const results = calculateAbTest(
        controlVisitors,
        controlConversions,
        variantVisitors,
        variantConversions,
        confidenceLevelInput,
        statisticalPowerInput
    );

    // Update the document with the results
    document.getElementById('conversion-rate-control').innerText = results.conversion_rate_control.toFixed(2) + '%';
    document.getElementById('conversion-rate-variant').innerText = results.conversion_rate_variant.toFixed(2) + '%';
    document.getElementById('lift').innerText = results.lift.toFixed(2) + '%';
    document.getElementById('difference-pp').innerText = results.difference_pp.toFixed(2) + ' pp';
    document.getElementById('confidence-interval-diff').innerText = 
        results.confidence_interval_difference_pp[0].toFixed(2) + '% to ' + results.confidence_interval_difference_pp[1].toFixed(2) + '%';
    document.getElementById('right-sided-interval').innerText = 
        results.right_sided_interval_pp[0].toFixed(2) + '% to ' + (results.right_sided_interval_pp[1] === Infinity ? '∞' : results.right_sided_interval_pp[1].toFixed(2) + '%');
    document.getElementById('left-sided-interval').innerText = 
        (results.left_sided_interval_pp[0] === -Infinity ? '-∞' : results.left_sided_interval_pp[0].toFixed(2)) + '% to ' + results.left_sided_interval_pp[1].toFixed(2) + '%';
    document.getElementById('value-plus-minus-SE').innerText = 
        results.value_plus_minus_95_SE_pp[0].toFixed(2) + ' ± ' + results.value_plus_minus_95_SE_pp[1].toFixed(2) + ' pp';
    document.getElementById('p-value').innerText = results.p_value.toFixed(4);
    document.getElementById('z-score').innerText = results.z_score.toFixed(2);
    document.getElementById('significance').innerText = results.p_value_1sided_significance;
    document.getElementById('sample-size').innerText = results.sample_size_per_group;
    document.getElementById('control-ci').innerText = 
        results.confidence_interval_control[0].toFixed(2) + '% to ' + results.confidence_interval_control[1].toFixed(2) + '%';
    document.getElementById('variant-ci').innerText = 
        results.confidence_interval_variant[0].toFixed(2) + '% to ' + results.confidence_interval_variant[1].toFixed(2) + '%';
    document.getElementById('relative-mde').innerText = results.relative_mde.toFixed(2) + '%';
    document.getElementById('bayesian-variant-wins').innerText = results.bayesian_variant_wins.toFixed(2) + '%';
    document.getElementById('bayesian-control-wins').innerText = results.bayesian_control_wins.toFixed(2) + '%';
    document.getElementById('bayes-factor').innerText = results.bayes_factor_H1_H0.toFixed(2);
}
function calculateAbTest(
    controlVisitors,
    controlConversions,
    variantVisitors,
    variantConversions,
    confidenceLevelInput,
    statisticalPowerInput
) {
    // Convert percentages to proportions
    let confidenceLevel = confidenceLevelInput / 100;
    let statisticalPower = statisticalPowerInput / 100;

    // Calculate conversion rates
    let conversionRateControl = controlConversions / controlVisitors;
    let conversionRateVariant = variantConversions / variantVisitors;

    // Handle edge cases - maybe this should be done elsewhere
    if (conversionRateControl <= 0 || conversionRateControl >= 1 || conversionRateVariant <= 0 || conversionRateVariant >= 1) {
        throw new Error('Conversion rates must be between 0% and 100% (exclusive).');
    }

    // Lets calculate p-value and z-score
    let pValueResult = calculatePValue(
        controlConversions,
        controlVisitors,
        variantConversions,
        variantVisitors,
        true  // One-sided test a default for AB testing
    );

    let zScore = pValueResult.z_score;
    let pValue = pValueResult.p_value;
    let pValue1SidedSignificance = pValue < (1 - confidenceLevel) ? 'Significant' : 'Not Significant';

    // Calculate Lift (relative difference)
    let absoluteDifference = conversionRateVariant - conversionRateControl;
    let lift = (absoluteDifference / conversionRateControl) * 100;

    // Difference in percentage points
    let differencePp = absoluteDifference * 100;

    // Calculate standard error for the difference
    let SE_diff = Math.sqrt(
        (conversionRateControl * (1 - conversionRateControl)) / controlVisitors +
        (conversionRateVariant * (1 - conversionRateVariant)) / variantVisitors
    );

    // Z-scores for confidence intervals
    let zAlphaTwoSided = jStat.normal.inv(1 - (1 - confidenceLevel) / 2, 0, 1);
    let zAlphaOneSided = jStat.normal.inv(confidenceLevel, 0, 1);

    // Two-sided confidence interval for the difference
    let lowerCiDiff = absoluteDifference - zAlphaTwoSided * SE_diff;
    let upperCiDiff = absoluteDifference + zAlphaTwoSided * SE_diff;
    let lowerCiDiffPp = lowerCiDiff * 100;
    let upperCiDiffPp = upperCiDiff * 100;

    // Right-sided (one-sided) confidence interval
    let lowerCiRightSided = absoluteDifference - zAlphaOneSided * SE_diff;
    let lowerCiRightSidedPp = lowerCiRightSided * 100;

    // Left-sided (one-sided) confidence interval
    let upperCiLeftSided = absoluteDifference + zAlphaOneSided * SE_diff;
    let upperCiLeftSidedPp = upperCiLeftSided * 100;

    // Value ± 95% SE
    let SE_95 = zAlphaTwoSided * SE_diff;
    let SE_95Pp = SE_95 * 100;
    let valuePlusMinusSE = [differencePp, SE_95Pp];

    // Sample size calculation
    let sampleSizePerGroup = calculateSampleSize(
        conversionRateControl,
        lift,
        confidenceLevelInput,
        statisticalPowerInput,
        2,
        true  // Should be one-sided tests
    );

    // Confidence intervals for control and variant
    let zAlpha = jStat.normal.inv(confidenceLevel, 0, 1);
    let controlCi = calculateConfidenceInterval(conversionRateControl, controlVisitors, zAlpha);
    let variantCi = calculateConfidenceInterval(conversionRateVariant, variantVisitors, zAlpha);

    // Relative MDE (Minimum Detectable Effect)
    // I am not a fan of this
    let relativeMde = calculateRelativeMde(
        controlVisitors,
        variantVisitors,
        conversionRateControl,
        confidenceLevel,
        statisticalPower
    );

    // Bayesian probability
    let bayesianResults = calculateBayesianProbability(
        controlConversions,
        controlVisitors,
        variantConversions,
        variantVisitors
    );

    // Prepare the results
    return {
        'conversion_rate_control': conversionRateControl * 100,  // As percentage
        'conversion_rate_variant': conversionRateVariant * 100,  // As percentage
        'lift': lift,  // In percentage
        'difference_pp': differencePp,  // Difference in percentage points
        'confidence_interval_difference_pp': [lowerCiDiffPp, upperCiDiffPp],  // 95% CI in percentage points
        'right_sided_interval_pp': [lowerCiRightSidedPp, Infinity],  // Right-sided interval
        'left_sided_interval_pp': [-Infinity, upperCiLeftSidedPp],  // Left-sided interval
        'value_plus_minus_95_SE_pp': valuePlusMinusSE,  // Value ± 95% SE
        'p_value': pValue,  // P-value for H0: B ≤ A
        'z_score': zScore,
        'p_value_1sided_significance': pValue1SidedSignificance,
        'sample_size_per_group': sampleSizePerGroup,
        'confidence_interval_control': [controlCi[0] * 100, controlCi[1] * 100],  // As percentage
        'confidence_interval_variant': [variantCi[0] * 100, variantCi[1] * 100],  // As percentage
        'relative_mde': relativeMde,  // In percentage
        'bayesian_variant_wins': bayesianResults.probability_variant_wins * 100,  // As percentage
        'bayesian_control_wins': bayesianResults.probability_control_wins * 100,  // As percentage
        'bayes_factor_H1_H0': bayesianResults.bayes_factor_H1_H0
    };
}

function calculateConfidenceInterval(conversionRate, visitors, zAlpha) {
    /*
    Calculate the confidence interval for a given conversion rate.

    Parameters:
        conversionRate (float): The conversion rate.
        visitors (int): The number of visitors.
        zAlpha (float): The z-score corresponding to the desired confidence level.

    Returns:
        Array[float, float]: The lower and upper bounds of the confidence interval.
    */

    let marginOfError = zAlpha * Math.sqrt((conversionRate * (1 - conversionRate)) / visitors);
    
    return [conversionRate - marginOfError, conversionRate + marginOfError];
}

function calculateRelativeMde(controlSampleSize, variantSampleSize, baselineConversionRate, confidenceLevel, statisticalPower) {
    // Check if confidence level and statistical power are in percentage
    if (confidenceLevel > 1) {
        confidenceLevel = confidenceLevel / 100;
    }
    if (statisticalPower > 1) {
        statisticalPower = statisticalPower / 100;
    }

    // Get z-scores for significance level and power (one-tailed test)
    let zAlpha = jStat.normal.inv(confidenceLevel, 0, 1); // One-tailed test
    let zBeta = jStat.normal.inv(statisticalPower, 0, 1);

    // Adjusted standard error for unequal sample sizes
    let se = Math.sqrt(
        (baselineConversionRate * (1 - baselineConversionRate) / controlSampleSize) +
        (baselineConversionRate * (1 - baselineConversionRate) / variantSampleSize)
    );

    // Calculate the absolute MDE (difference) using z-scores and standard error
    let absoluteMde = (zAlpha + zBeta) * se;

    // Calculate the relative MDE as a percentage of the baseline conversion rate
    let relativeMdePercentage = (absoluteMde / baselineConversionRate) * 100;

    return relativeMdePercentage;
}

function durationCalculator(
    baselineConversionRatePercentage,
    mdePercentage,
    numberOfVariants,
    dailyVisitorsTotal,
    confidenceLevel = 95,
    power = 80,
    isOneSided = true
) {
    // Convert percentages to proportions
    let p1 = baselineConversionRatePercentage / 100;
    let mde = mdePercentage / 100;
    confidenceLevel = confidenceLevel / 100;
    power = power / 100;

    // Calculate p2 based on MDE
    let p2 = p1 * (1 + mde);
    p2 = Math.min(p2, 0.9999);  // Ensure p2 is less than 1

    // Calculate z-scores
    let alpha = 1 - confidenceLevel;
    let zAlpha;
    if (isOneSided) {
        zAlpha = jStat.normal.inv(1 - alpha, 0, 1);
    } else {
        zAlpha = jStat.normal.inv(1 - alpha / 2, 0, 1);
    }
    let zBeta = jStat.normal.inv(power, 0, 1);

    // Adjust alpha for multiple variants if needed (Bonferroni correction)
    if (numberOfVariants > 2) {
        let alphaAdjusted = alpha / (numberOfVariants - 1);
        zAlpha = jStat.normal.inv(1 - alphaAdjusted, 0, 1);
    }

    // Calculate pooled standard deviation
    let pAvg = (p1 + p2) / 2;
    let sigmaPooled = Math.sqrt(2 * pAvg * (1 - pAvg));

    // Calculate numerator and denominator
    let se = Math.sqrt(p1 * (1 - p1) + p2 * (1 - p2));
    let numerator = Math.pow(zAlpha * sigmaPooled + zBeta * se, 2);
    let deltaP = p2 - p1;
    let denominator = Math.pow(deltaP, 2);

    // Calculate required sample size per group
    let sampleSizePerGroup = Math.ceil(numerator / denominator);

    // Calculate daily visitors per group
    let dailyVisitorsPerGroup = dailyVisitorsTotal / numberOfVariants;

    // Calculate number of days needed per group
    let daysNeeded = Math.ceil(sampleSizePerGroup / dailyVisitorsPerGroup);

    return daysNeeded;
}


function calculatePValue(
    controlConversions,
    controlVisitors,
    variantConversions,
    variantVisitors,
    isOneSided = true
) {
    /**
     * Calculate the p-value for an A/B test.
     *
     * Parameters:
     *   controlConversions (int): Number of conversions in the control group.
     *   controlVisitors (int): Number of visitors in the control group.
     *   variantConversions (int): Number of conversions in the variant group.
     *   variantVisitors (int): Number of visitors in the variant group.
     *   isOneSided (bool): Whether the test is one-sided (default is true).
     *
     * Returns:
     *   Object: An object containing the z-score, p-value, and isSignificant (boolean).
     */

    // Calculate conversion rates
    let conversionRateControl = controlConversions / controlVisitors;
    let conversionRateVariant = variantConversions / variantVisitors;

    // Pooled conversion rate across both groups
    let pooledConversionRate = (controlConversions + variantConversions) / (controlVisitors + variantVisitors);

    // Pooled standard error
    let pooledSE = Math.sqrt(
        pooledConversionRate * (1 - pooledConversionRate) * (1 / controlVisitors + 1 / variantVisitors)
    );

    // Handle cases where pooledSE is zero to avoid division by zero
    if (pooledSE === 0) {
        throw new Error("Pooled standard error is zero. Check your input values.");
    }

    // Calculate z-score
    let zScore = (conversionRateVariant - conversionRateControl) / pooledSE;

    // Calculate p-value based on one-sided or two-sided test
    let pValue;
    if (isOneSided) {
        pValue = 1 - jStat.normal.cdf(zScore, 0, 1);
    } else {
        pValue = 2 * (1 - jStat.normal.cdf(Math.abs(zScore), 0, 1));
    }

    // Ensure p-value is between 0 and 1
    pValue = Math.min(Math.max(pValue, 0), 1);

    return {
        z_score: zScore,
        p_value: pValue,
        is_significant: pValue < 0.05,
    };
}

function weeklyMdeAndSignificanceCalculator(
    baselineConversionRatePercentage,
    numberOfVariants,
    dailyVisitorsTotal,
    confidenceLevel = 95,
    power = 80,
    weeks = 7
) {
    // Convert baseline conversion rate to a proportion
    let baselineConversionRate = baselineConversionRatePercentage / 100;
    confidenceLevel = confidenceLevel / 100;
    power = power / 100;

    // Calculate daily visitors per group
    let dailyVisitorsPerGroup = dailyVisitorsTotal / numberOfVariants;

    let results = [];

    // Loop through each week and calculate MDE and significance
    for (let week = 1; week <= weeks; week++) {
        // Total visitors for the current week
        let weeklyVisitorsPerGroup = dailyVisitorsPerGroup * 7 * week;

        // Calculate the relative MDE for the current week
        let relativeMde = calculateRelativeMde(
            weeklyVisitorsPerGroup,
            weeklyVisitorsPerGroup,
            baselineConversionRate,
            confidenceLevel,
            power
        );

        // Calculate pooled standard error for significance testing
        let se = Math.sqrt(
            baselineConversionRate * (1 - baselineConversionRate) *
            (1 / weeklyVisitorsPerGroup + 1 / weeklyVisitorsPerGroup)
        );

        // Calculate absolute difference needed for significance (MDE)
        let zAlpha = jStat.normal.inv(1 - (1 - confidenceLevel), 0, 1); // One-tailed test
        let requiredDifference = zAlpha * se;

        let thisWeeksConversionRateControl = baselineConversionRate * weeklyVisitorsPerGroup;
        let thisWeeksConversionRateVariant = 
            (baselineConversionRate + (relativeMde / 100) * baselineConversionRate) * weeklyVisitorsPerGroup;

        // Calculate p-value for the current week
        let pValueResult = calculatePValue(
            Math.round(thisWeeksConversionRateControl),
            Math.round(weeklyVisitorsPerGroup),
            Math.round(thisWeeksConversionRateVariant),
            Math.round(weeklyVisitorsPerGroup),
            true
        );

        // Append result for the week
        results.push({
            week: week,
            relativeMde: relativeMde,
            requiredDifferenceForSignificance: requiredDifference * 100, // Convert to percentage
            pValue: pValueResult.pValue,
            significance: pValueResult.isSignificant
        });
    }

    return results;
}

function calculateSampleSize(p1, liftPercentage, confidenceLevel, power, numVariants = 2, isOneSided = true) {
    /*
    Calculate the required sample size for an A/B test.

    Parameters:
        p1 (float): The baseline conversion rate (proportion) for the control group.
        liftPercentage (float): The expected percentage increase in conversion rate for the treatment group.
        confidenceLevel (float): The desired confidence level (e.g., 95 for 95% confidence).
        power (float): The desired statistical power (e.g., 80 for 80% power).
        numVariants (int, optional): The number of variants being tested (default is 2 for control and one variant).
        isOneSided (bool, optional): Whether the test is one-sided (default is True).

    Returns:
        int: The required sample size per group for the A/B test.
    */

    // Calculate p2 based on lift percentage
    let p2 = p1 * (1 + liftPercentage / 100);

    // Ensure p2 does not exceed 1
    if (p2 >= 1) {
        p2 = 0.9999;
    }

    // Calculate alpha and z-scores
    let alpha = 1 - confidenceLevel / 100;
    let zAlpha;
    
    if (isOneSided) {
        zAlpha = jStat.normal.inv(1 - alpha, 0, 1);
    } else {
        zAlpha = jStat.normal.inv(1 - alpha / 2, 0, 1);
    }
    
    let zBeta = jStat.normal.inv(power / 100, 0, 1);

    // Adjust alpha for multiple variants (Bonferroni correction)
    if (numVariants > 2) {
        let alphaAdjusted = alpha / (numVariants - 1);
        zAlpha = jStat.normal.inv(1 - alphaAdjusted, 0, 1);
    }

    // Pooled proportion
    let pAvg = (p1 + p2) / 2;
    let sigmaPooled = Math.sqrt(2 * pAvg * (1 - pAvg));

    // Standard error of proportions
    let standardError = Math.sqrt(p1 * (1 - p1) + p2 * (1 - p2));

    // Calculate numerator and denominator for the sample size formula
    let numerator = Math.pow(zAlpha * sigmaPooled + zBeta * standardError, 2);
    let deltaP = p2 - p1;
    let denominator = Math.pow(deltaP, 2);

    // Calculate sample size per group
    let sampleSizePerGroup = Math.ceil(numerator / denominator);

    return sampleSizePerGroup;
}

function calculateBayesianProbability(controlConversions, controlVisitors, variantConversions, variantVisitors, simulations = 100000) {
    /*
    Calculate the Bayesian probability that the variant is better than the control.

    Parameters:
        controlConversions (int): The number of conversions in the control group.
        controlVisitors (int): The number of visitors in the control group.
        variantConversions (int): The number of conversions in the variant group.
        variantVisitors (int): The number of visitors in the variant group.
        simulations (int): The number of Monte Carlo simulations to run (default is 100,000).

    Returns:
        Object: A dictionary containing the probability that the variant is better than the control and the Bayes Factor.
    */

    // Beta prior parameters (assume non-informative priors)
    let alphaPrior = 1;
    let betaPrior = 1;

    // Posterior distributions
    let controlPosteriorAlpha = alphaPrior + controlConversions;
    let controlPosteriorBeta = betaPrior + (controlVisitors - controlConversions);
    let variantPosteriorAlpha = alphaPrior + variantConversions;
    let variantPosteriorBeta = betaPrior + (variantVisitors - variantConversions);

    // Monte Carlo simulation to estimate P(Variant > Control)
    let variantWins = 0;
    
    for (let i = 0; i < simulations; i++) {
        let controlSample = jStat.beta.sample(controlPosteriorAlpha, controlPosteriorBeta);
        let variantSample = jStat.beta.sample(variantPosteriorAlpha, variantPosteriorBeta);

        if (variantSample > controlSample) {
            variantWins++;
        }
    }

    let probabilityVariantWins = variantWins / simulations;
    let probabilityControlWins = 1 - probabilityVariantWins;

    // Log Marginal Likelihood under H1
    let ln_ml_H1 = (
        jStat.betaln(controlConversions + alphaPrior, controlVisitors - controlConversions + betaPrior) -
        jStat.betaln(alphaPrior, betaPrior) +
        jStat.betaln(variantConversions + alphaPrior, variantVisitors - variantConversions + betaPrior) -
        jStat.betaln(alphaPrior, betaPrior)
    );

    // Log Marginal Likelihood under H0
    let totalConversions = controlConversions + variantConversions;
    let totalVisitors = controlVisitors + variantVisitors;
    let ln_ml_H0 = (
        jStat.betaln(totalConversions + alphaPrior, totalVisitors - totalConversions + betaPrior) -
        jStat.betaln(alphaPrior, betaPrior)
    );

    // Compute the Bayes Factor (BF_10)
    let ln_BF_10 = ln_ml_H1 - ln_ml_H0;
    let BF_10 = Math.exp(ln_BF_10);

    // Handle potential issues with BF_10
    if (!isFinite(BF_10)) {
        BF_10 = ln_BF_10 > 0 ? Infinity : 0;
    }

    return {
        probability_variant_wins: probabilityVariantWins,
        probability_control_wins: probabilityControlWins,
        bayes_factor_H1_H0: BF_10
    };
}