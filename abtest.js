function calculateBayesianProbability(controlConversions, controlVisitors, variantConversions, variantVisitors) {
    // Beta prior parameters (assume non-informative priors)
    const alphaPrior = 1;
    const betaPrior = 1;

    // Posterior distributions
    const controlPosteriorAlpha = alphaPrior + controlConversions;
    const controlPosteriorBeta = betaPrior + (controlVisitors - controlConversions);
    const variantPosteriorAlpha = alphaPrior + variantConversions;
    const variantPosteriorBeta = betaPrior + (variantVisitors - variantConversions);

    // Monte Carlo simulation to estimate P(Variant > Control)
    const simulations = 10000;  // Reduced number of simulations for performance
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

    return {
        probabilityVariantWins: probabilityVariantWins,
        probabilityControlWins: probabilityControlWins
    };
}

function calculateABTest() {
    // Clear previous results
    clearResults();

    // Gather inputs
    const controlVisitors = parseInt(document.getElementById('control-visitors').value);
    const controlConversions = parseInt(document.getElementById('control-conversions').value);
    const variantVisitors = parseInt(document.getElementById('variant-visitors').value);
    const variantConversions = parseInt(document.getElementById('variant-conversions').value);
    const confidenceLevelInput = parseFloat(document.getElementById('confidence-level').value);
    const statisticalPowerInput = parseFloat(document.getElementById('statistical-power').value);
    const desiredEffectSizePercentage = parseFloat(document.getElementById('desired-effect-size').value); 

    const confidenceLevel = confidenceLevelInput / 100;
    const statisticalPower = statisticalPowerInput / 100;
    const desiredEffectSize = desiredEffectSizePercentage / 100;

    // Calculate conversion rates
    const conversionRateControl = controlConversions / controlVisitors;
    const conversionRateVariant = variantConversions / variantVisitors;

    // Handle edge cases
    if (conversionRateControl === 0 || conversionRateControl === 1 || conversionRateVariant === 0 || conversionRateVariant === 1) {
        alert('Conversion rates of 0% or 100% are not supported.');
        return;
    }

    // Display conversion rates
    document.getElementById('conversion-rate-control').innerText = `Conversion Rate (Control): ${(conversionRateControl * 100).toFixed(2)}%`;
    document.getElementById('conversion-rate-variant').innerText = `Conversion Rate (Variant): ${(conversionRateVariant * 100).toFixed(2)}%`;

    // Calculate absolute difference
    const absoluteDifference = conversionRateVariant - conversionRateControl;

    // Display Lift (relative difference)
    const lift = (absoluteDifference / conversionRateControl) * 100;
    document.getElementById('lift').innerText = `Lift: ${lift.toFixed(2)}%`;

    // Calculate pooled standard error for absolute difference
    const pooledProb = (controlConversions + variantConversions) / (controlVisitors + variantVisitors);
    const se = Math.sqrt(pooledProb * (1 - pooledProb) * (1 / controlVisitors + 1 / variantVisitors));

    // Check for small sample sizes
    if (controlVisitors + variantVisitors < 30) {
        alert('Warning: Small sample sizes may not provide accurate results using normal approximation.');
    }

    // Calculate z-score using absolute difference
    const z = absoluteDifference / se;

    // Calculate two-sided p-value using cumulative normal distribution function
    const pValue2Sided = 2 * (1 - jStat.normal.cdf(Math.abs(z), 0, 1));
    const pValue2SidedSignificance = pValue2Sided < (1 - confidenceLevel) ? 'Significant' : 'Not Significant';

    // Display p-value for two-sided test
    if (pValue2Sided < 0.001) {
        document.getElementById('p-value-2-side').innerText = `P-value (2-sided): <0.001`;
    } else {
        document.getElementById('p-value-2-side').innerText = `P-value (2-sided): ${pValue2Sided.toFixed(5)}`;
    }
    // Display p-value significance
    document.getElementById('significance-2-sided').innerText = `P-value (2-sided) Significance: ${pValue2SidedSignificance}`;

    // Calculate zAlpha and zBeta for sample size calculation (two-tailed test)
    const zAlpha = jStat.normal.inv(1 - (1 - confidenceLevel) / 2, 0, 1); // For two-tailed test
    const zBeta = jStat.normal.inv(statisticalPower, 0, 1);

    // Calculate required sample size using desired effect size
    const p = pooledProb;  // Use pooled probability for estimated conversion rate
    const sampleSizePerGroup = Math.ceil(((zAlpha + zBeta) ** 2 * 2 * p * (1 - p)) / (desiredEffectSize ** 2));

    document.getElementById('sample-size').innerText = `Required Sample Size (per group): ${sampleSizePerGroup}`;

    // Calculate confidence intervals for control and variant
    const controlCI = calculateConfidenceInterval(conversionRateControl, controlVisitors, zAlpha);
    const variantCI = calculateConfidenceInterval(conversionRateVariant, variantVisitors, zAlpha);

    // Display confidence intervals
    document.getElementById('confidence-interval-control').innerText = `Confidence Interval (Control): ${(controlCI[0] * 100).toFixed(2)}% - ${(controlCI[1] * 100).toFixed(2)}%`;
    document.getElementById('confidence-interval-variant').innerText = `Confidence Interval (Variant): ${(variantCI[0] * 100).toFixed(2)}% - ${(variantCI[1] * 100).toFixed(2)}%`;

    // Calculate and display relative MDE (Minimum Detectable Effect)
    const relativeMDE = calculateRelativeMDE(
        controlVisitors,                // Control sample size
        variantVisitors,                // Variant sample size
        conversionRateControl,          // Baseline conversion rate (from control group)
        confidenceLevel,                // Confidence level
        statisticalPower                // Statistical power
    );
    document.getElementById('mde').innerText = `Minimum Detectable Effect (MDE - Lift): ${relativeMDE.toFixed(2)}%`;

    // Calculate Bayesian probability
    const bayesianResults = calculateBayesianProbability(controlConversions, controlVisitors, variantConversions, variantVisitors);
    document.getElementById('bayesian-variant-wins').innerText = `Bayesian P(Variant > Control): ${(bayesianResults.probabilityVariantWins * 100).toFixed(2)}%`;
    document.getElementById('bayesian-control-wins').innerText = `Bayesian P(Control > Variant): ${(bayesianResults.probabilityControlWins * 100).toFixed(2)}%`;
}

function calculateConfidenceInterval(conversionRate, visitors, zAlpha) {
    const marginOfError = zAlpha * Math.sqrt((conversionRate * (1 - conversionRate)) / visitors);
    return [conversionRate - marginOfError, conversionRate + marginOfError];
}

function calculateRelativeMDE(controlSampleSize, variantSampleSize, baselineConversionRate, confidenceLevel, statisticalPower) {
    // Get z-scores for significance level and power (two-tailed test)
    const zAlpha = jStat.normal.inv(1 - (1 - confidenceLevel) / 2, 0, 1);
    const zBeta = jStat.normal.inv(statisticalPower, 0, 1);

    // Adjusted standard error for unequal sample sizes
    const se = Math.sqrt(
        (baselineConversionRate * (1 - baselineConversionRate) / controlSampleSize) +
        (baselineConversionRate * (1 - baselineConversionRate) / variantSampleSize)
    );

    // Calculate the absolute MDE (difference) using z-scores and standard error
    const absoluteMDE = (zAlpha + zBeta) * se;

    // Calculate the relative MDE as a percentage of the baseline conversion rate
    const relativeMDEPercentage = (absoluteMDE / baselineConversionRate) * 100;
    
    return relativeMDEPercentage;
}

function clearResults() {
    document.getElementById('conversion-rate-control').innerText = '';
    document.getElementById('conversion-rate-variant').innerText = '';
    document.getElementById('lift').innerText = '';  // Clear the lift value
    document.getElementById('p-value-2-side').innerText = '';  // Clear two-sided p-value
    document.getElementById('significance-2-sided').innerText = '';
    document.getElementById('sample-size').innerText = '';
    document.getElementById('confidence-interval-control').innerText = '';
    document.getElementById('confidence-interval-variant').innerText = '';
    document.getElementById('mde').innerText = '';
    document.getElementById('bayesian-variant-wins').innerText = '';
    document.getElementById('bayesian-control-wins').innerText = '';
}

// Event listener for calculating when the button is clicked
document.getElementById('calculate-button').addEventListener('click', calculateABTest);
