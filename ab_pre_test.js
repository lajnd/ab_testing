function fetchFormAndCalculateMDE() {
    // Fetch input values from the form
    const weeklyTraffic = parseInt(document.getElementById('weekly-traffic').value);
    const weeklyConversions = parseInt(document.getElementById('weekly-conversions').value);
    const confidenceLevel = parseFloat(document.getElementById('confidence-level').value);
    const power = parseFloat(document.getElementById('statistical-power').value);
    const numVariants = parseInt(document.getElementById('number-of-variants').value);

    const weeks = [1, 2, 3, 4, 5, 6]; // Default weeks

    // Call calculateMDEOneSided with the fetched values
    const results = calculateMDEOneSided(
        weeklyTraffic,
        weeklyConversions,
        confidenceLevel / 100,
        power / 100,
        numVariants,
        weeks,
        false // Use Dunnett's adjustment
    );

    // calculate baseline conversion rate in percenage
    const baselineCr = (weeklyConversions / weeklyTraffic) * 100;
    // update the baseline conversion rate in the UI
    document.getElementById('baseline-cr').innerText = baselineCr.toFixed(2) + '%';
    document.getElementById('mde_week_1').innerText = results[0].relativeMde.toFixed(2) + '%';
    document.getElementById('visitors_week_1').innerText = results[0].totalUsersPerVariant;
    document.getElementById('mde_week_2').innerText = results[1].relativeMde.toFixed(2) + '%';
    document.getElementById('visitors_week_2').innerText = results[1].totalUsersPerVariant;
    document.getElementById('mde_week_3').innerText = results[2].relativeMde.toFixed(2) + '%';
    document.getElementById('visitors_week_3').innerText = results[2].totalUsersPerVariant;
    document.getElementById('mde_week_4').innerText = results[3].relativeMde.toFixed(2) + '%';
    document.getElementById('visitors_week_4').innerText = results[3].totalUsersPerVariant;
    document.getElementById('mde_week_5').innerText = results[4].relativeMde.toFixed(2) + '%';
    document.getElementById('visitors_week_5').innerText = results[4].totalUsersPerVariant;
    document.getElementById('mde_week_6').innerText = results[5].relativeMde.toFixed(2) + '%';
    document.getElementById('visitors_week_6').innerText = results[5].totalUsersPerVariant;
}

function calculateMDEOneSided(
    weeklyTraffic,
    weeklyConversions,
    confidenceLevel,
    power,
    numVariants, // Includes control
    weeks = [1, 2, 3, 4, 5, 6],
    useDunnettsAdjustment = false
) {
    // Total number of groups including control
    var totalGroups = numVariants; // This includes the control group

    // Number of test variants (exclude control)
    var testVariants = numVariants - 1;

    // Adjust confidence level using Dunnett's adjustment if needed
    var adjustedConfidenceLevel = confidenceLevel;
    if (useDunnettsAdjustment && testVariants > 1) {
        adjustedConfidenceLevel = calculateDunnettCriticalValue(confidenceLevel, testVariants);
    }

    // Convert confidence level and power to Z-scores for one-sided test
    var zAlpha = jStat.normal.inv(1 - (1 - adjustedConfidenceLevel), 0, 1); // One-sided test
    var zBeta = jStat.normal.inv(power, 0, 1);

    var baselineCr = weeklyConversions / weeklyTraffic; // Baseline conversion rate

    var mdeList = [];

    for (var i = 0; i < weeks.length; i++) {
        var w = weeks[i];
        
        // Step 1: Calculate the number of visitors per variant for each week
        var visitorsPerVariant = (weeklyTraffic * w) / numVariants; // Divide by the total number of variants (including control)

        // Step 2: Calculate the variance of the outcome metric (OEC)
        var varianceOEC = baselineCr * (1 - baselineCr);

        // Step 3: Calculate the absolute MDE based on variance and Z-scores using visitors per variant
        var mdeAbsolute = (zAlpha + zBeta) * Math.sqrt(varianceOEC / visitorsPerVariant);

        // Step 4: Calculate the relative MDE as a percentage
        var relativeMde = (mdeAbsolute / baselineCr) * 100;

        // Step 5: Recalculate the sample size based on the calculated MDE
        var recalculatedSampleSize = calculateSampleSize(varianceOEC, mdeAbsolute, power * 100, testVariants, useDunnettsAdjustment);

        // Print the recalculated sample size and the original sample size for comparison
        console.log(`Week ${w} -> Original Sample Size: ${visitorsPerVariant}, Recalculated Sample Size: ${recalculatedSampleSize}`);

        // Append week, absolute MDE, relative MDE, visitors per test variant, and recalculated sample size
        mdeList.push({
            week: w,
            mde: mdeAbsolute,
            relativeMde: relativeMde,
            totalUsersPerVariant: visitorsPerVariant, // Per test variant
            recalculatedSampleSize: recalculatedSampleSize,
        });
    }

    return mdeList;
}

// Function to calculate the sample size based on variance, MDE, and power level
function calculateSampleSize(variance, mde, powerLevel = 80, numVariants = 1, conservative = false) {
    let coefficient;

    // Choose the coefficient based on the desired power level
    if (powerLevel === 80) {
        coefficient = 16;  // For 80% power
    } else if (powerLevel === 90) {
        coefficient = 21;  // For 90% power
    } else {
        throw new Error("Unsupported power level. Please choose either 80 or 90.");
    }

    // If using Wheeler's conservative formula for multiple comparisons, adjust the coefficient
    if (conservative && numVariants > 1) {
        coefficient = Math.pow(4 * numVariants, 2);  // Wheeler's conservative formula
    }

    // Calculate the sample size
    const sampleSize = (coefficient * variance) / Math.pow(mde, 2);

    return Math.ceil(sampleSize); // Round up to the nearest integer
}

// Function to calculate Dunnett's critical value for multiple comparisons
function calculateDunnettCriticalValue(confidenceLevel, numComparisons) {
    // For Dunnett's test with one-sided adjustment, use the student t-distribution
    const df = numComparisons;  // Degrees of freedom approximated to number of comparisons
    const criticalValue = jStat.studentt.inv(confidenceLevel, df); // One-sided t-distribution

    return criticalValue;
}

function calculateMDE(sampleSize, variance, powerLevel = 80, numVariants = 1, conservative = false) {
    let coefficient;

    // Choose the coefficient based on the desired power level
    if (powerLevel === 80) {
        coefficient = 16;  // For 80% power
    } else if (powerLevel === 90) {
        coefficient = 21;  // For 90% power
    } else {
        throw new Error("Unsupported power level. Please choose either 80 or 90.");
    }

    // If using Wheeler's conservative formula for multiple comparisons, adjust the coefficient
    if (conservative && numVariants > 1) {
        coefficient = Math.pow(4 * numVariants, 2);  // Wheeler's conservative formula
    }

    // Reverse the formula to calculate MDE
    const mde = Math.sqrt((coefficient * variance) / sampleSize);

    return mde; // Return the MDE
}


function testCalculateSampleSizeAndMDE() {
    const variance = 0.02; // Example variance
    const powerLevel = 80; // Desired power level (80%)
    const numVariants = 2; // Number of variants (excluding control)
    const conservative = false; // Whether to use conservative formula

    // Step 1: Calculate sample size based on a known MDE
    const mde = 0.05; // Example MDE (5%)
    const calculatedSampleSize = calculateSampleSize(variance, mde, powerLevel, numVariants, conservative);
    
    // Step 2: Use the calculated sample size to reverse-calculate the MDE
    const recalculatedMDE = calculateMDE(calculatedSampleSize, variance, powerLevel, numVariants, conservative);

    // Step 3: Test equivalence (allowing for small floating-point differences)
    const epsilon = 1e-6; // Small tolerance for floating-point comparisons
    const mdeEquivalent = Math.abs(mde - recalculatedMDE) < epsilon;

    if (mdeEquivalent) {
        console.log("Test passed: calculateSampleSize and calculateMDE are equivalent.");
    } else {
        console.log(`Test failed: Original MDE = ${mde}, Recalculated MDE = ${recalculatedMDE}`);
    }

    // Repeat the process but starting with sample size
    // Step 1: Calculate MDE based on a known sample size
    const sampleSize = 10000; // Example sample size
    const calculatedMDE = calculateMDE(sampleSize, variance, powerLevel, numVariants, conservative);

    // Step 2: Use the calculated MDE to reverse-calculate the sample size
    const recalculatedSampleSize = calculateSampleSize(variance, calculatedMDE, powerLevel, numVariants, conservative);

    // Step 3: Test equivalence (allowing for small floating-point differences)
    const sampleSizeEquivalent = Math.abs(sampleSize - recalculatedSampleSize) < epsilon;

    if (sampleSizeEquivalent) {
        console.log("Test passed: calculateMDE and calculateSampleSize are equivalent.");
    } else {
        console.log(`Test failed: Original Sample Size = ${sampleSize}, Recalculated Sample Size = ${recalculatedSampleSize}`);
    }
}

// Call the test function
testCalculateSampleSizeAndMDE();