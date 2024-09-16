const jStat = require('jstat');

function calculateBayesianProbability(controlConversions, controlVisitors, variantConversions, variantVisitors) {
    const alphaPrior = 1;
    const betaPrior = 1;

    const controlPosteriorAlpha = alphaPrior + controlConversions;
    const controlPosteriorBeta = betaPrior + (controlVisitors - controlConversions);
    const variantPosteriorAlpha = alphaPrior + variantConversions;
    const variantPosteriorBeta = betaPrior + (variantVisitors - variantConversions);

    const simulations = 100000;
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

// Helper function to compute log of binomial coefficient
function logBinomialCoefficient(n, k) {
    if (k > n) return -Infinity;
    let logCoeff = 0;
    for (let i = 1; i <= k; i++) {
        logCoeff += Math.log(n - (k - i)) - Math.log(i);
    }
    return logCoeff;
}

// Helper function to compute log-sum-exp in a numerically stable way
function logsumexp(a, b) {
    if (a === -Infinity) return b; // log(0) + b = b
    if (b === -Infinity) return a; // log(0) + a = a
    if (a > b) {
        return a + Math.log(1 + Math.exp(b - a));
    } else {
        return b + Math.log(1 + Math.exp(a - b));
    }
}

// Function to compute the probability that Variant > Control using logarithms
function calculateBayesianProbabilityAnalytically(controlConversions, controlVisitors, variantConversions, variantVisitors) {
    const alphaPrior = 1;
    const betaPrior = 1;

    // Directly assign to alpha1, beta1 for control group and alpha2, beta2 for variant group
    const alpha1 = alphaPrior + controlConversions;
    const beta1 = betaPrior + (controlVisitors - controlConversions);
    const alpha2 = alphaPrior + variantConversions;
    const beta2 = betaPrior + (variantVisitors - variantConversions);

    // Closed-form solution to compute log(P(Variant > Control))
    let logProbabilityVariantWins = -Infinity; // Initialize as log(0)

    for (let i = 0; i < alpha2; i++) {
        const logBinomCoeff = logBinomialCoefficient(alpha2 + beta2 - 1, i);
        const logBetaTerm = Math.log(jStat.betafn(alpha1 + i, beta1 + beta2));
        logProbabilityVariantWins = logsumexp(logProbabilityVariantWins, logBinomCoeff + logBetaTerm);
    }

    const logDenom = Math.log(jStat.betafn(alpha1, beta1)) + Math.log(jStat.betafn(alpha2, beta2));
    const logProbabilityControlWins = Math.log(1 - Math.exp(logProbabilityVariantWins - logDenom));

    return {
        probabilityVariantWins: Math.exp(logProbabilityVariantWins - logDenom),
        probabilityControlWins: Math.exp(logProbabilityControlWins)
    };
}


function testCalculateBayesianProbability() {
    // Test case 1: Equal conversions and visitors
    let result = calculateBayesianProbability(50, 1000, 50, 1000);
    console.log("Test Case 1:", result);
    console.assert(result.probabilityVariantWins > 0.45 && result.probabilityVariantWins < 0.55, 
                   "Test case 1 failed: Expected variant win probability to be around 50%.");
    
    // Test case 2: Variant performs better
    result = calculateBayesianProbability(100, 1000, 150, 1000);
    console.log("Test Case 2:", result);
    console.assert(result.probabilityVariantWins > 0.95, 
                   "Test case 2 failed: Expected variant win probability to be greater than 95%.");

    // Test case 3: Control performs better
    result = calculateBayesianProbability(150, 1000, 100, 1000);
    console.log("Test Case 3:", result);
    console.assert(result.probabilityControlWins > 0.95, 
                   "Test case 3 failed: Expected control win probability to be greater than 95%.");
    
    // Test case 4: Very small data
    result = calculateBayesianProbability(1, 10, 9, 10);
    console.log("Test Case 4:", result);
    console.assert(result.probabilityVariantWins > 0.9, 
                   "Test case 4 failed: Expected variant win probability to be greater than 90%.");

    // Run the same tests for calculateBayesianProbabilityAnalytically
    // Test case 1: Equal conversions and visitors
    result = calculateBayesianProbabilityAnalytically(50, 1000, 50, 1000);
    console.log("Test Case 1 (Analytical):", result);
    console.assert(result.probabilityVariantWins > 0.45 && result.probabilityVariantWins < 0.55, 
                   "Test case 1 (Analytical) failed: Expected variant win probability to be around 50%.");
            
    // Test case 2: Variant performs better
    result = calculateBayesianProbabilityAnalytically(100, 1000, 150, 1000);
    console.log("Test Case 2 (Analytical):", result);
    console.assert(result.probabilityVariantWins > 0.95, 
                   "Test case 2 (Analytical) failed: Expected variant win probability to be greater than 95%.");
    
    // Test case 3: Control performs better
    result = calculateBayesianProbabilityAnalytically(150, 1000, 100, 1000);
    console.log("Test Case 3 (Analytical):", result);
    console.assert(result.probabilityControlWins > 0.95, 
                   "Test case 3 (Analytical) failed: Expected control win probability to be greater than 95%.");

    // Test case 4: Very small data
    result = calculateBayesianProbabilityAnalytically(1, 10, 9, 10);
    console.log("Test Case 4 (Analytical):", result);
    console.assert(result.probabilityVariantWins > 0.9, 
                   "Test case 4 (Analytical) failed: Expected variant win probability to be greater than 90%.");
    
    console.log("All test cases completed.");
}

// Run test cases
testCalculateBayesianProbability();