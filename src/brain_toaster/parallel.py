from joblib import Parallel, delayed
import numpy as np
from scipy.optimize import brentq
from scipy.stats import nct

def ncpci(t_stat, df, confLevel):  
    """  
    Compute confidence interval for the non-centrality parameter (NCP)  
    for a t-distribution.  

    Parameters:  
      t_stat   : observed t-statistic value  
      df       : degrees of freedom  
      confLevel: confidence level (e.g., 0.95)  

    Returns:  
      A NumPy array [ncp_lower, ncp_upper]  
    """  
    alpha = 1 - confLevel  

    def lower_bound(ncp):  
        return nct.cdf(t_stat, df, ncp) - alpha / 2  

    def upper_bound(ncp):  
        return nct.cdf(t_stat, df, ncp) - (1 - alpha / 2)  

    # Find lower bound  
    try:  
        if lower_bound(-100) * lower_bound(t_stat) > 0:  
            # print(f"No sign change for lower_bound: f(-100)={lower_bound(-100)}, f(t_stat)={lower_bound(t_stat)}")  
            ncp_lower = np.nan  # Assign NaN if no root exists  
        else:  
            ncp_lower = brentq(lower_bound, -100, t_stat)  
    except ValueError as e:  
        print(f"Error finding lower bound: {e}")  
        ncp_lower = np.nan  

    # Find upper bound  
    try:  
        if upper_bound(t_stat) * upper_bound(100) > 0:  
            # print(f"No sign change for upper_bound: f(t_stat)={upper_bound(t_stat)}, f(100)={upper_bound(100)}")  
            ncp_upper = np.nan  # Assign NaN if no root exists  
        else:  
            ncp_upper = brentq(upper_bound, t_stat, 100)  
    except ValueError as e:  
        print(f"Error finding upper bound: {e}")  
        ncp_upper = np.nan  

    return np.array([ncp_lower, ncp_upper])

def compute_single_ci(t_val, td_fac, DoF, confLevel):  
    """  
    Compute confidence interval for a single t-value.  

    Parameters:  
      t_val     : t-statistic value  
      td_fac    : t-to-d conversion factor  
      DoF       : degrees of freedom  
      confLevel : confidence level (e.g., 0.95)  

    Returns:  
      A tuple (ci_lower, ci_upper)  
    """  
    # If t value is zero (or nearly so), assign zero CIs.  
    if np.abs(t_val) < np.finfo(float).eps:  
        return (0.0, 0.0)  
    else:  
        ci_tmp = ncpci(t_val, DoF, confLevel)  
        if np.isnan(ci_tmp).any():  
            # print(f"Warning: NaN confidence interval for t_val={t_val}")  
            return (np.nan, np.nan)  
        # Multiply by td_fac to obtain effect size confidence interval  
        return (ci_tmp[0] * td_fac, ci_tmp[1] * td_fac)

def compute_ci_for_all(ts, td_fac, DoF, confLevel):
    # Parallel computation for each t-value in the vector ts.
    results = Parallel(n_jobs=-1)(
        delayed(compute_single_ci)(t_val, td_fac, DoF, confLevel)
        for t_val in ts
    )
    # Convert list of tuples to a (2, n) NumPy array.
    results = np.array(results)  # shape: (num_voxels, 2)
    return results.T  # Transpose to shape (2, num_voxels)
