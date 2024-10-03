function fetchFormAndCalculate() {
    try {
        // Fetch and parse input values from the form, convert percentages to decimals
        const confidenceLevel = Number(document.getElementById('confidence-level').value) / 100;
        const powerLevel = Number(document.getElementById('power-level').value) / 100;
        const conversionRate = Number(document.getElementById('conversion-rate').value) / 100;
        const minEffect = Number(document.getElementById('min-effect').value) / 100;
        const weeklyTraffic = Number(document.getElementById('weekly-traffic').value);

        // Fetch number of variants
        const numVariants = Number(document.getElementById('num-variants').value);

        // Calculate alpha
        const alpha = 1.0 - confidenceLevel;

        // Fetch radio button values for isRelative and oneSided
        const isRelative = document.querySelector('input[name="is-relative"]:checked').value === "true";
        const oneSided = document.querySelector('input[name="one-sided"]:checked').value === "true";

        // Apply Bonferroni correction by dividing alpha by the number of comparisons
        const adjustedAlpha = alpha / numVariants;

        // Calculate the required sample size per variant
        const sampleSizePerVariant = calculateSampleSize(adjustedAlpha, powerLevel, conversionRate, minEffect, isRelative, oneSided);

        // Total sample size for all variants plus the control
        const totalSampleSize = sampleSizePerVariant * (numVariants + 1); // +1 for the control group

        // Calculate the duration in weeks based on the weekly traffic
        const durationInWeeks = totalSampleSize / weeklyTraffic;
        // take ceiling to ensure we have enough samples
        const durationInWeeksCeil = Math.ceil(durationInWeeks);

        // Update the result
        document.getElementById('sample-size-result-per-group').innerText = `${sampleSizePerVariant}`;
        document.getElementById('sample-size-result-total').innerText = `${totalSampleSize}`;
        document.getElementById('sample-size-result-weeks').innerText = `${durationInWeeksCeil}`;
    } catch (error) {
        console.error(error);
        alert(`An error occurred: ${error.message}`);
    }
}

/**
 * Calculate the required sample size based on a given significance level (alpha),
 * desired power (1 - beta), base conversion rate, and minimum detectable effect (MDE).
 * 
 * @param {number} alpha - The significance level (Type I error rate).
 * @param {number} powerLevel - The desired power level (1 - Type II error rate).
 * @param {number} conversionRate - The base conversion rate (p).
 * @param {number} minEffect - The minimum detectable effect (MDE).
 * @param {boolean} isRelative - Whether the MDE is relative to the base conversion rate.
 * @param {boolean} oneSided - Whether the test is one-sided.
 * 
 * @returns {number} - The required sample size.
 */
function calculateSampleSize(alpha, powerLevel, conversionRate, minEffect, isRelative = true, oneSided = true) {
    // Adjust for relative or absolute difference
    let delta = isRelative ? minEffect * conversionRate : minEffect;

    // Ensure the conversion rate is within the correct range
    if (conversionRate > 0.5) {
        conversionRate = 1.0 - conversionRate;
    }

    // Critical values for alpha and power level (one-sided or two-sided)
    let tAlpha = oneSided ? ppf(1.0 - alpha) : ppf(1.0 - alpha / 2);
    let tBeta = ppf(powerLevel);

    // Standard deviations for control and test groups
    // We could do pooling here if we assume equal variances, but lets not
    let sd1 = Math.sqrt(2 * conversionRate * (1.0 - conversionRate));
    let sd2 = Math.sqrt(conversionRate * (1.0 - conversionRate) + (conversionRate + delta) * (1.0 - (conversionRate + delta)));

    // Sample size calculation based on the standard t-test formula
    let sampleSize = Math.pow((tAlpha * sd1 + tBeta * sd2), 2) / Math.pow(delta, 2);

    return Math.ceil(sampleSize);
}

/**
 * Compute the inverse of the cumulative distribution function (CDF) for the normal distribution.
 * 
 * @param {number} p - The probability corresponding to the cumulative distribution function.
 * @returns {number} - The Z-value corresponding to the given probability.
 */
function ppf(p) {
    return jStat.normal.inv(p, 0, 1);
}