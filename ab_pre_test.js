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
        confidenceLevel/100,
        power/100,
        numVariants,
        weeks
    );

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
    numVariants,
    weeks = [1, 2, 3, 4, 5, 6]
) {
    // Convert confidence level and power to Z-scores for one-sided test
    var zAlpha = jStat.normal.inv(confidenceLevel, 0, 1); // One-sided test
    var zBeta = jStat.normal.inv(power, 0, 1);

    var baselineCr = weeklyConversions / weeklyTraffic; // Baseline conversion rate

    var mdeList = [];

    for (var i = 0; i < weeks.length; i++) {
        var w = weeks[i];
        // Number of users per variant after w weeks
        var totalUsersPerVariant = (weeklyTraffic / numVariants) * w;

        // Calculate standard error
        var se = Math.sqrt((2 * baselineCr * (1 - baselineCr)) / totalUsersPerVariant);

        // Calculate MDE for w weeks (absolute)
        var mde = (zAlpha + zBeta) * se;

        // Calculate relative MDE as a percentage of the baseline conversion rate
        var relativeMde = (mde / baselineCr) * 100;

        // Append week, absolute MDE, relative MDE, and visitors per variant
        mdeList.push({
            week: w,
            mde: mde,
            relativeMde: relativeMde,
            totalUsersPerVariant: totalUsersPerVariant,
        });
    }

    return mdeList;
}