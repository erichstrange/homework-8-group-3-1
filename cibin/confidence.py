"""
Implementation of 2 sided CI of treatment effect.
"""

import math
import numpy as np
import scipy as sp
from scipy.special import comb
from itertools import filterfalse, combinations


def tau_twosided_ci(n11, n10, n01, n00, alpha, exact=True,
                    max_combinations=10**5, reps=10**3):
    '''
    Find 2-sided 1−alpha confidence bounds for the average treatment effect.

    Assuming a randomized experiment with binary outcomes and two treatments,
    (active and control), finds the confidence bounds for tau, the average
    treatment effect.

    Parameters
    ----------
    n11 : int
        number of subjects assigned to treatment 1 who had outcome 1
    n10 : int
        number of subjects assigned to treatment 1 who had outcome 0
    n01 : int
        number of subjects assigned to treatment 0 who had outcome 1
    n00 : int
        number of subjects assigned to treatment 0 who had outcome 0
    alpha : float in (0, 1)
        The desired type 1 error level.
    exact: boolean
        Whether to compute exact result. If true, sampling distributions are
        computed by enumerating all possible samples. If false, simulations
        are used to approximate the sampling distribution.
    max_combinations: int
        Maximum amount of combinations to sample per table if exact=True
    reps: int
        Number of samples to generate if exact=False

    Returns
    -------
    Tuple of 3 items:
        [lb, ub]: lower/upper bound of the confidence interval,
        [allocation that gives lb, allocation that gives ub],
        [# tables examined, total reps across simulations]
    '''
    N = n11 + n10 + n01 + n00
    n = n11 + n10                     # size of treatment sample
    t_star = n11/n - n01/(N-n)        # unbiased estimate of tau

    n_combs = int(comb(N, n))         # total number of samples for exact ans
    if exact and n_combs > max_combinations:
        raise ValueError(f"Please raise max_combinations to {n_combs} for \
                          exact solution.")
        
    if ((alpha <= 0) or (alpha >= 1)):
        raise ValueError("Invalid value for alpha!")
    
    if (n11 < 0) or (n10 < 0) or (n01 < 0) or (n00 < 0):
        raise ValueError("subject count cannot be negative!")

    conf_set = {}
    n_tables, n_reps = 0, 0

    for Nt in N_generator(N, n00, n01, n10, n11):
        table = potential_outcomes(Nt)    # generate potential outcomes table
        tau = (Nt[1]-Nt[2])/N             # average treatment effect for table
        t = abs(t_star - tau)             # test statistic

        if exact:
            # generate samples     
            mask = np.zeros((n_combs,N), dtype=bool)
            for i, sample in enumerate(combinations(range(N), n)):
                mask[i,sample] = True

            # calculate test statistic for each sample
            tau_hat = mask.dot(table[:, 1])/n - (~mask).dot(table[:, 0])/(N-n)
            dist = abs(tau_hat-tau)
            n_reps += n_combs

            # include in confidence set if t is within confidence bounds
            if t <= np.percentile(dist, (1-alpha)*100):
                conf_set[tau] = Nt

        else:
            mask = np.zeros((reps,N), dtype=bool)
            for i in range(reps):
                sample = np.random.choice(range(N),n, replace=False)
                mask[i,sample] = True

            tau_hat = mask.dot(table[:,1])/n - (~mask).dot(table[:,0])/(N-n)
            dist = abs(tau_hat-tau)
            n_reps += reps

            if t <= np.percentile(dist, (1-alpha)*100):
                conf_set[tau] = Nt

        n_tables += 1

    lower, upper = min(conf_set.keys()), max(conf_set.keys())
    lower_alloc, upper_alloc = conf_set[lower], conf_set[upper]
    return ([lower, upper], [lower_alloc, upper_alloc], [n_tables, n_reps])


def N_generator(N, n00, n01, n10, n11):
    '''
    Generate tables algebraically consistent with data from an experiment
    with binary outcomes.

    Parameters
    ----------
    N : int
        number of subjects
    n00 : int
        number of subjects assigned to treatment 0 who had outcome 0
    n01 : int
        number of subjects assigned to treatment 0 who had outcome 1
    n10 : int
        number of subjects assigned to treatment 1 who had outcome 0
    n11 : int
        number of subjects assigned to treatment 1 who had outcome 1

    Returns
    -------
    Nt : list of 4 ints
    N00, subjects with potential outcome 0 under control and treatment
    N01, subjects with potential outcome 0 under control and 1 under treatment
    N10, subjects with potential outcome 1 under control and 0 under treatment
    N11, subjects with potential outcome 1 under control and treatment
    '''
    
    if N < (n00 + n01 + n10 + n11):
        raise ValueError("Number of subjects do not match!")
    if (n11 < 0) or (n10 < 0) or (n01 < 0) or (n00 < 0):
        raise ValueError("subject count cannot be negative!")
    for i in range(N):
        N00 = i
        for j in range(N-i):
            N01 = j
            for k in range(N-i-j):
                N10 = k
                N11 = N-i-j-k
                if filterTable([N00, N01, N10, N11], n00, n01, n10, n11):
                    yield [N00, N01, N10, N11]
                else:
                    pass


def filterTable(Nt, n00, n01, n10, n11):
    '''
    Check whether summary table Nt of binary outcomes is consistent with
    observed counts.

    Implements the test in Theorem 1 of Li and Ding (2016).

    Parameters:
    ----------
    Nt : list of four ints
        the table of counts of subjects with each combination of potential
        outcomes, in the order N_00, N_01, N_10, N_11
    n00 : int
        number of subjects assigned to control whose observed response was 0
    n01 : int
        number of subjects assigned to control whose observed response was 1
    n10 : int
        number of subjects assigned to treatment whose observed response was 0
    n11 : int
        number of subjects assigned to treatment whose observed response was 1

    Returns:
    --------
    ok : boolean
        True if table is consistent with the data
    '''
    
    if N < (n00 + n01 + n10 + n11):
        raise ValueError("Number of subjects do not match!")
    if (n11 < 0) or (n10 < 0) or (n01 < 0) or (n00 < 0):
        raise ValueError("subject count cannot be negative!")
    N = np.sum(Nt)   # total subjects
    return max(0, n11-Nt[1], Nt[3]-n01, Nt[2]+Nt[3]-n10-n01) <= \
        min(Nt[3], n11, Nt[2]+Nt[3]-n01, N-Nt[1]-n01-n10)


def potential_outcomes(Nt):
    '''
    make a 2xN table of potential outcomes from the 2x2 summary table Nt.

    Parameters
    ----------
    Nt : list of 4 ints
        N00, N01, N10, N11

    Returns
    -------
    po : Nx2 table of potential outcomes consistent with Nt
    '''
    if len(Nt) != 4:
        raise ValueError("Table size must be 4: N00, N01, N10, N11 ")
        
    for i in range(len(Nt)):
        if Nt[i] < 0:
            raise ValueError("Cannot have a negative number as a potential outcome")
    return np.reshape(np.array([0, 0] * Nt[0] + [0, 1] * Nt[1] +
                               [1, 0] * Nt[2] + [1, 1] * Nt[3]), [-1, 2])
