function fetchAndCalculate() {
    try {
        // Fetch data from the HTML input fields
        const weeklyTraffic = Number(document.getElementById('weekly-traffic').value);
        const weeklyConversions = Number(document.getElementById('weekly-conversions').value);
        const confidenceLevel = Number(document.getElementById('confidence-level').value) / 100;  // Convert to decimal
        const powerLevel = Number(document.getElementById('power-level').value) / 100;  // Convert to decimal
        const oneSided = document.querySelector('input[name="test-type"]:checked').value === "one-sided";  // Get test type
        const numConditions = 2;  // Assuming 2 conditions for now, will fix later

        // Call the calculateMDEFor6Weeks function to calculate the MDE for each week
        const mdeList = calculateMDEFor6Weeks(weeklyTraffic, weeklyConversions, numConditions, confidenceLevel, powerLevel, true, oneSided);

        // Update the table with the results for each week
        for (let week = 1; week <= mdeList.length; week++) {
            const mde = mdeList[week - 1] * 100;  // Convert MDE to percentage
            const visitorsPerVariant = (weeklyTraffic * week) / numConditions;

            // Get the table elements and check if they exist
            const weekCell = document.getElementById(`week-${week}`);
            const mdeCell = document.getElementById(`mde-week-${week}`);
            const visitorsCell = document.getElementById(`visitors-week-${week}`);

            if (weekCell && mdeCell && visitorsCell) {
                // Update the respective table cells if they exist
                weekCell.innerText = week;
                mdeCell.innerText = `${mde.toFixed(2)}%`;  // Express MDE as percentage
                visitorsCell.innerText = visitorsPerVariant.toFixed(0);  // Round visitors to whole number
            } else {
                console.error(`Table cell for week ${week} is missing.`);
            }
        }
    } catch (error) {
        console.error(error);
        alert(`An error occurred: ${error.message}`);
    }
}

/**
 * Calculate the Minimum Detectable Effect (MDE) for 6 weeks based on weekly traffic and conversions.
 * 
 * @param {number} weeklyTraffic - The total traffic per week.
 * @param {number} weeklyConversions - The total number of conversions per week.
 * @param {number} numConditions - The number of conditions (default is 2).
 * @param {number} alpha - Significance level (Type I error rate, default is 0.05).
 * @param {number} powerLevel - Power level (1 - Type II error rate, default is 0.8).
 * @param {boolean} isRelative - Whether the MDE is relative or absolute (default is true).
 * @param {boolean} oneSided - Whether the test is one-sided or two-sided (default is true).
 * 
 * @returns {Array} - A list of MDE values for each of the 6 weeks.
 */
function calculateMDEFor6Weeks(weeklyTraffic, weeklyConversions, numConditions = 2, alpha = 0.05, powerLevel = 0.8, isRelative = true, oneSided = true) {
    // Calculate the conversion rate per week
    const conversionRate = weeklyConversions / weeklyTraffic;
    
    // Initialize variables to store results
    let mdeResults = [];
    
    // Loop through each week and calculate the MDE
    for (let week = 1; week <= 6; week++) {
        // Calculate the total traffic for the current week, per variant (split between conditions)
        const visitorsPerVariant = (weeklyTraffic * week) / numConditions;

        // Calculate the MDE for the current week 
        const mde = calculateMDE(alpha, powerLevel, conversionRate, visitorsPerVariant, isRelative, oneSided);

        // Add the result to the list
        mdeResults.push(mde);
    }
    
    return mdeResults;
}

/**
 * Helper function to compute the inverse of the cumulative distribution function (CDF) for the normal distribution.
 * 
 * @param {number} p - The probability corresponding to the cumulative distribution function.
 * @returns {number} - The Z-value corresponding to the given probability.
 */
function ppnd(p) {
    return jStat.normal.inv(p, 0, 1);
}

/**
 * Calculate the minimum detectable effect (MDE) based on the input parameters.
 * 
 * @param {number} alpha - Significance level (Type I error rate).
 * @param {number} powerLevel - Power level (1 - Type II error rate).
 * @param {number} conversionRate - The base conversion rate (p).
 * @param {number} numSubjects - The number of subjects for each group.
 * @param {boolean} [isRelative=true] - Whether the MDE is relative or absolute.
 * @param {boolean} [oneSided=true] - Whether the test is one-sided or two-sided.
 * 
 * @returns {number} - The minimum detectable effect (MDE) for the given sample size.
 */
function calculateMDE(alpha, powerLevel, conversionRate, numSubjects, isRelative = true, oneSided = true) {
    // Ensure the conversion rate is within the correct range
    if (conversionRate > 0.5) {
        conversionRate = 1.0 - conversionRate;
    }

    // Critical values for alpha and power level (one-sided or two-sided)
    const tAlpha = oneSided ? ppnd(1.0 - alpha) : ppnd(1.0 - alpha / 2); // Z value for significance level
    const tBeta = ppnd(powerLevel); // Z value for power level

    // Standard deviation for the control group
    const sd1 = Math.sqrt(2 * conversionRate * (1.0 - conversionRate));

    // Calculate the minimum detectable effect (MDE)
    const delta = Math.sqrt((tAlpha * sd1 + tBeta * sd1) ** 2 / numSubjects);

    // Return the MDE as relative if requested (i.e., percentage change)
    return isRelative ? delta / conversionRate : delta;
}