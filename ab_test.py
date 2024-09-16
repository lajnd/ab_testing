import math
from scipy.stats import norm, beta
import numpy as np

def duration_calculator(
    baseline_conversion_rate_percentage: float,
    mde_percentage: float,
    number_of_variants: int,
    daily_visitors_total: int,
    confidence_level=95,
    power=80,
    is_one_sided=True
    ) -> int:
    """
    Calculate the duration (in days) required for an A/B test.

    Parameters:
    - baseline_conversion_rate_percentage (float): The baseline conversion rate as a percentage.
    - mde_percentage (float): The minimum detectable effect as a percentage.
    - number_of_variants (int): The number of variants being tested.
    - daily_visitors_total (int): The total number of daily visitors across all groups.
    - confidence_level (float, optional): The desired confidence level for the test (default is 95).
    - power (float, optional): The desired statistical power of the test (default is 80).
    - is_one_sided (bool, optional): Whether the test is one-sided (default is True).

    Returns:
    - int: The number of days needed to run the A/B test.
    """

    # Convert percentages to proportions
    p1 = baseline_conversion_rate_percentage / 100
    mde = mde_percentage / 100
    confidence_level = confidence_level / 100
    power = power / 100
    
    # Calculate p2 based on MDE 
    # relative difference MDE
    p2 = p1 * (1 + mde)
    if p2 >= 1:
        p2 = 0.9999  # Ensure p2 is less than 1
    
    # Calculate z-scores
    alpha = 1 - confidence_level
    if is_one_sided:
        z_alpha = norm.ppf(1 - alpha)
    else:
        z_alpha = norm.ppf(1 - alpha / 2)
    z_beta = norm.ppf(power)
    
    # Adjust alpha for multiple variants if needed
    # Doing Bonferroni correction here, which is a bit conservative
    if number_of_variants > 2:
        alpha_adjusted = alpha / (number_of_variants - 1)
        z_alpha = norm.ppf(1 - alpha_adjusted)
    
    # Calculate pooled standard deviation
    # TODO: figure out if this is the right way to calculate pooled standard deviation
    p_avg = (p1 + p2) / 2
    sigma_pooled = math.sqrt(2 * p_avg * (1 - p_avg))
    
    # Calculate numerator and denominator
    se = math.sqrt(p1 * (1 - p1) + p2 * (1 - p2))
    numerator = (z_alpha * sigma_pooled + z_beta * se) ** 2
    delta_p = p2 - p1
    denominator = delta_p ** 2
    
    # Calculate required sample size per group
    sample_size_per_group = math.ceil(numerator / denominator)
    
    # Calculate daily visitors per group
    daily_visitors_per_group = daily_visitors_total / number_of_variants
    
    # Calculate number of days needed per group
    days_needed = math.ceil(sample_size_per_group / daily_visitors_per_group)
    
    return days_needed

def calculate_p_value(
    control_conversions: int,
    control_visitors: int,
    variant_conversions: int,
    variant_visitors: int,
    is_one_sided: bool = True
) -> dict:
    """
    Calculate the p-value for an A/B test.
    
    Parameters:
        control_conversions: int - Number of conversions in the control group.
        control_visitors: int - Number of visitors in the control group.
        variant_conversions: int - Number of conversions in the variant group.
        variant_visitors: int - Number of visitors in the variant group.
        is_one_sided: bool - Whether the test is one-sided or two-sided.
        
    Returns:
        dict: A dictionary containing the z-score and p-value.
    """
    
    # Calculate conversion rates
    conversion_rate_control = control_conversions / control_visitors
    conversion_rate_variant = variant_conversions / variant_visitors
    
    # Pooled conversion rate across both groups
    pooled_conversion_rate = (control_conversions + variant_conversions) / (control_visitors + variant_visitors)
    
    # Pooled standard error
    pooled_se = math.sqrt(pooled_conversion_rate * (1 - pooled_conversion_rate) * (1 / control_visitors + 1 / variant_visitors))
    
    # Calculate z-score
    z_score = (conversion_rate_variant - conversion_rate_control) / pooled_se
    
    # Calculate p-value based on one-sided or two-sided test
    if is_one_sided:
        p_value = 1 - norm.cdf(z_score)
    else:
        p_value = 2 * (1 - norm.cdf(abs(z_score)))
    
    return {
        'z_score': z_score,
        'p_value': p_value,
        'is_significant': p_value < 0.05,
    }

def weekly_mde_and_significance_calculator(
    baseline_conversion_rate_percentage: float,
    number_of_variants: int,
    daily_visitors_total: int,
    confidence_level=95,
    power=80,
    weeks=7
    ) -> list:

    """
    Calculates the Minimum Detectable Effect (MDE) and significance for an A/B test over a specified number of weeks.
    Parameters:
    - baseline_conversion_rate_percentage (float): The baseline conversion rate as a percentage.
    - number_of_variants (int): The number of variants being tested.
    - daily_visitors_total (int): The total number of daily visitors across all groups.
    - confidence_level (float, optional): The desired confidence level for the test (default is 95).
    - power (float, optional): The desired statistical power of the test (default is 80).
    - weeks (int, optional): The number of weeks over which to conduct the test (default is 7).
    Returns:
    - List[Dict[str, Union[int, float]]]: A list of dictionaries containing the results for each week, including:
        - week (int): The week number.
        - relative_mde (float): The relative minimum detectable effect for the week.
        - required_difference_for_significance (float): The absolute difference required for significance, expressed as a percentage.
        - p_value (float): The p-value calculated for the week.
        - significance (bool): Indicates whether the result is statistically significant.
    """
    # Convert baseline conversion rate to a proportion
    baseline_conversion_rate = baseline_conversion_rate_percentage / 100
    confidence_level = confidence_level / 100
    power = power / 100

    # Calculate daily visitors per group
    daily_visitors_per_group = daily_visitors_total / number_of_variants

    results = []

    # Loop through each week and calculate MDE and significance
    for week in range(1, weeks + 1):
        # Total visitors for the current week
        weekly_visitors_per_group = daily_visitors_per_group * 7 * week

        # Calculate the relative MDE for the current week
        relative_mde = calculate_relative_mde(
            control_sample_size=int(weekly_visitors_per_group),
            variant_sample_size=int(weekly_visitors_per_group),
            baseline_conversion_rate=baseline_conversion_rate,
            confidence_level=confidence_level,
            statistical_power=power
        )

        # Calculate pooled standard error for significance testing
        se = math.sqrt(baseline_conversion_rate * (1 - baseline_conversion_rate) * (1 / weekly_visitors_per_group + 1 / weekly_visitors_per_group))

        # Calculate absolute difference needed for significance (MDE)
        z_alpha = norm.ppf(1 - (1 - confidence_level))  # One-tailed test
        required_difference = z_alpha * se

        this_weeks_conversion_rate_control = baseline_conversion_rate * weekly_visitors_per_group
        this_weeks_conversion_rate_variant = (baseline_conversion_rate + (relative_mde / 100) * baseline_conversion_rate) * weekly_visitors_per_group
        

        # Calculate p-value for the current week
        p_value_result = calculate_p_value(
            control_conversions=int(this_weeks_conversion_rate_control),
            control_visitors=int(weekly_visitors_per_group),
            variant_conversions=int(this_weeks_conversion_rate_variant),
            variant_visitors=int(weekly_visitors_per_group),
            is_one_sided=True
        )

        # Append result for the week
        results.append({
            'week': week,
            'relative_mde': relative_mde,
            'required_difference_for_significance': required_difference * 100,  # Convert to percentage
            'p_value': p_value_result['p_value'],
            'significance': p_value_result['is_significant']  # Use significance result from calculate_p_value
        })

    return results

def calculate_confidence_interval(conversion_rate: float, 
                                  visitors: int, 
                                  z_alpha: float) -> tuple:
    """
    Calculate the confidence interval for a given conversion rate.
    Parameters:
        conversion_rate (float): The conversion rate.
        visitors (int): The number of visitors.
        z_alpha (float): The z-score corresponding to the desired confidence level.
    Returns:
        Tuple[float, float]: The lower and upper bounds of the confidence interval.
    """
    margin_of_error: float = z_alpha * math.sqrt((conversion_rate * (1 - conversion_rate)) / visitors)
    return (conversion_rate - margin_of_error, conversion_rate + margin_of_error)


def calculate_relative_mde(control_sample_size, variant_sample_size, baseline_conversion_rate, confidence_level, statistical_power):
    # check if confidence level and statistical power are in percentage
    if confidence_level > 1:
        confidence_level = confidence_level / 100
    if statistical_power > 1:
        statistical_power = statistical_power / 100
    
    # Get z-scores for significance level and power (one-tailed test)
    z_alpha = norm.ppf(confidence_level)  # One-tailed test
    z_beta = norm.ppf(statistical_power)

    # Adjusted standard error for unequal sample sizes
    se = math.sqrt(
        (baseline_conversion_rate * (1 - baseline_conversion_rate) / control_sample_size) +
        (baseline_conversion_rate * (1 - baseline_conversion_rate) / variant_sample_size)
    )

    # Calculate the absolute MDE (difference) using z-scores and standard error
    absolute_mde = (z_alpha + z_beta) * se

    # Calculate the relative MDE as a percentage of the baseline conversion rate
    relative_mde_percentage = (absolute_mde / baseline_conversion_rate) * 100

    return relative_mde_percentage

def calculate_bayesian_probability(control_conversions: int, 
                                   control_visitors: int, 
                                   variant_conversions: int, 
                                   variant_visitors: int,
                                   simulations=100_000) -> dict:
    # Beta prior parameters (assume non-informative priors)
    alpha_prior = 1
    beta_prior = 1

    # Posterior distributions
    control_posterior_alpha = alpha_prior + control_conversions
    control_posterior_beta = beta_prior + (control_visitors - control_conversions)
    variant_posterior_alpha = alpha_prior + variant_conversions
    variant_posterior_beta = beta_prior + (variant_visitors - variant_conversions)

    # Monte Carlo simulation to estimate P(Variant > Control)
    variant_wins = 0
    for _ in range(simulations):
        control_sample = beta.rvs(control_posterior_alpha, control_posterior_beta)
        variant_sample = beta.rvs(variant_posterior_alpha, variant_posterior_beta)

        if variant_sample > control_sample:
            variant_wins += 1

    probability_variant_wins = variant_wins / simulations
    probability_control_wins = 1 - probability_variant_wins

    return {
        'probability_variant_wins': probability_variant_wins,
        'probability_control_wins': probability_control_wins
    }

def calculate_sample_size(
    p1: float,
    lift_percentage: float,
    confidence_level: float,
    power: float,
    num_variants=2,
    is_one_sided=True
) -> int:
    """
    Calculate the required sample size for an A/B test.
    I have compared this method with the ones described here:
    https://towardsdatascience.com/required-sample-size-for-a-b-testing-6f6608dd330a

    And I get the same results.
    
    Parameters:
        p1 (float): The baseline conversion rate (proportion) for the control group.
        lift_percentage (float): The expected percentage increase in conversion rate for the treatment group.
        confidence_level (float): The desired confidence level (e.g., 95 for 95% confidence).
        power (float): The desired statistical power (e.g., 80 for 80% power).
        num_variants (int, optional): The number of variants being tested (default is 2 for control and one variant).
        is_one_sided (bool, optional): Whether the test is one-sided (default is True).
    Returns:
        int: The required sample size per group for the A/B test.
    """
    # Calculate p2 based on lift percentage
    p2 = p1 * (1 + lift_percentage / 100)

    # Ensure p2 does not exceed 1
    if p2 >= 1:
        p2 = 0.9999

    # Calculate z-scores
    alpha = 1 - confidence_level / 100
    if is_one_sided:
        z_alpha = norm.ppf(1 - alpha)
    else:
        z_alpha = norm.ppf(1 - alpha / 2)
    z_beta = norm.ppf(power / 100)

    # Adjust alpha for multiple variants if needed
    # I am doing a Bonferroni correction here
    if num_variants > 2:
        alpha_adjusted = alpha / (num_variants - 1)
        z_alpha = norm.ppf(1 - alpha_adjusted)

    # Pooled proportion
    p_avg = (p1 + p2) / 2
    sigma_pooled = math.sqrt(2 * p_avg * (1 - p_avg))

    # Numerator and denominator
    # standard error of proportions 
    standard_error = math.sqrt(p1 * (1 - p1) + p2 * (1 - p2))
    
    numerator = (z_alpha * sigma_pooled + z_beta * standard_error) ** 2
    delta_p = p2 - p1
    denominator = delta_p ** 2
    sample_size_per_group = math.ceil(numerator / denominator)

    return sample_size_per_group

def calculate_ab_test(
    control_visitors: int,
    control_conversions: int,
    variant_visitors: int,
    variant_conversions: int,
    confidence_level_input: float,
    statistical_power_input: float,
    ) -> dict:
    

    # Convert percentages to proportions
    confidence_level = confidence_level_input / 100
    statistical_power = statistical_power_input / 100

    # Calculate conversion rates
    conversion_rate_control = control_conversions / control_visitors
    conversion_rate_variant = variant_conversions / variant_visitors

    # Handle edge cases
    if conversion_rate_control <= 0 or conversion_rate_control >= 1 or conversion_rate_variant <= 0 or conversion_rate_variant >= 1:
        raise ValueError('Conversion rates must be between 0% and 100% (exclusive).')

    # Calculate p-value and z-score
    p_value_result = calculate_p_value(
        control_conversions=control_conversions,
        control_visitors=control_visitors,
        variant_conversions=variant_conversions,
        variant_visitors=variant_visitors,
        is_one_sided=True
    )
    
    z_score = p_value_result['z_score']
    p_value_1sided = p_value_result['p_value']
    p_value_1sided_significance = 'Significant' if p_value_1sided < (1 - confidence_level) else 'Not Significant'

    # Calculate Lift (relative difference)
    absolute_difference = conversion_rate_variant - conversion_rate_control
    lift = (absolute_difference / conversion_rate_control) * 100

    # Update sample size calculation
    sample_size_per_group = calculate_sample_size(
        p1=conversion_rate_control,
        lift_percentage=lift,
        confidence_level=confidence_level_input,
        power=statistical_power_input,
        num_variants=len([control_visitors, variant_visitors]),
        is_one_sided=True
    )

    # Calculate confidence intervals for control and variant
    z_alpha = norm.ppf(confidence_level)
    control_ci = calculate_confidence_interval(conversion_rate_control, control_visitors, z_alpha)
    variant_ci = calculate_confidence_interval(conversion_rate_variant, variant_visitors, z_alpha)

    # Calculate relative MDE (Minimum Detectable Effect)
    relative_mde = calculate_relative_mde(
        control_visitors,
        variant_visitors,
        conversion_rate_control,
        confidence_level,
        statistical_power
    )

    # Calculate Bayesian probability
    bayesian_results = calculate_bayesian_probability(
        control_conversions,
        control_visitors,
        variant_conversions,
        variant_visitors
    )

    # Prepare the results
    results = {
        'conversion_rate_control': conversion_rate_control * 100,  # As percentage
        'conversion_rate_variant': conversion_rate_variant * 100,  # As percentage
        'lift': lift,  # In percentage
        'p_value_1sided': p_value_1sided,
        'p_value_1sided_significance': p_value_1sided_significance,
        'z_score': z_score,
        'sample_size_per_group': sample_size_per_group,
        'confidence_interval_control': (control_ci[0] * 100, control_ci[1] * 100),  # As percentage
        'confidence_interval_variant': (variant_ci[0] * 100, variant_ci[1] * 100),  # As percentage
        'relative_mde': relative_mde,  # In percentage
        'bayesian_variant_wins': bayesian_results['probability_variant_wins'] * 100,  # As percentage
        'bayesian_control_wins': bayesian_results['probability_control_wins'] * 100,  # As percentage
    }

    return results

# Example usage
if __name__ == "__main__":
    # test duration calculator
    number_of_days = duration_calculator(
                    baseline_conversion_rate_percentage=9.09,
                    mde_percentage=30,
                    number_of_variants=2,
                    daily_visitors_total=150,
                    )
    print("Number of days required for A/B test:", number_of_days)
    
    # Input parameters
    control_visitors = 1100
    control_conversions = 100
    variant_visitors = 1100
    variant_conversions = 140
    confidence_level_input = 95  # in percentage
    statistical_power_input = 80  # in percentage

    try:
        # Calculate A/B test statistics
        results = calculate_ab_test(
            control_visitors,
            control_conversions,
            variant_visitors,
            variant_conversions,
            confidence_level_input,
            statistical_power_input,
        )

        # Display the results
        print("Conversion Rate (Control): {:.2f}%".format(results['conversion_rate_control']))
        print("Conversion Rate (Variant): {:.2f}%".format(results['conversion_rate_variant']))
        print("Lift: {:.2f}%".format(results['lift']))
        print("P-value (1-sided): {:.5f}".format(results['p_value_1sided']))
        print("Significance: {}".format(results['p_value_1sided_significance']))
        print("Required Sample Size per Group: {}".format(results['sample_size_per_group']))
        print("Confidence Interval (Control): {:.2f}% - {:.2f}%".format(*results['confidence_interval_control']))
        print("Confidence Interval (Variant): {:.2f}% - {:.2f}%".format(*results['confidence_interval_variant']))
        print("Minimum Detectable Effect (MDE): {:.2f}%".format(results['relative_mde']))
        print("Bayesian P(Variant > Control): {:.2f}%".format(results['bayesian_variant_wins']))
        print("Bayesian P(Control > Variant): {:.2f}%".format(results['bayesian_control_wins']))

    except ValueError as e:
        print("Error:", e)