"""
Implementation of 2 sided CI of treatment effect.
"""

import math
import numpy as np
import scipy as sp
from scipy.special import comb
from itertools import filterfalse, combinations
from scipy.stats import hypergeom


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

    n_combs = int(comb(N, n,exact=True))         # total number of samples for exact ans
    if exact and n_combs > max_combinations:
        raise ValueError(f"Please raise max_combinations to {n_combs} for \
                          #exact solution.")
        
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
    
    if sum(Nt) < (n00 + n01 + n10 + n11):
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


def hypergeom_accept(N, G, n, cl=0.975, randomized=False):
    """
    Acceptance region for a randomized hypergeometric test for sterne method
    If randomized==True, find the acceptance region for a randomized, exact 
    level-alpha test of the null hypothesis X~Hypergeom(N,G,n). 
    The acceptance region is the smallest possible. (And not, for instance, symmetric.)
    If randomized==False, find the smallest conservative acceptance region.
    Parameters
    ----------
    N : integer
        number of independent trials
    G : int
        number of "good" in population
    n : integer
        number of sample
    cl : float
        confidence level  
    ramndomized : Boolean
        return randomized exact test or conservative non-randomized test?
    Returns
    --------
    If randomized:
    I : list
        values for which the test never rejects
    If not randomized:
    I : list
        values for which the test does not reject
    """
    assert n <= N, 'impossible sample size'
    assert G <= N, 'impossible number of "good"'
    assert 0 < cl < 1, 'silly confidence level'

    alpha = 1 - cl
    prob = np.arange(0, n+1)
    I = list(prob)
    pmf = hypergeom.pmf(prob, N, G, n)
    bottom = 0
    top = n
    J = []
    prob_J = 0
    prob_tail = 0
    while prob_tail < alpha:
        p_bot = pmf[bottom]
        p_top = pmf[top]
        if p_bot < p_top:
            J = [bottom]
            prob_J = p_bot
            bottom += 1
        elif p_bot > p_top:
            J = [top]
            prob_J = p_top
            top -= 1
        else:
            if bottom < top:
                J = [bottom, top]
                prob_J = p_bot + p_top
                bottom += 1
                top -= 1
            else:
                J = [bottom]
                prob_J = p_bot
                bottom += 1
        prob_tail += prob_J
        for j in J:
            I.remove(j)
    if randomized == False:
        while prob_tail > alpha:
            j = J.pop()
            prob_tail -= pmf[j]
            I.append(j)
    return I


def hypergeom_conf_interval(n, x, N, cl=0.975, alternative="two-sided", G=None):
    """
    Confidence interval for a hypergeometric distribution parameter G, the number of good
    objects in a population in size N, based on the number x of good objects in a simple
    random sample of size n.
    Parameters
    ----------
    n : int
        The number of draws without replacement.
    x : int
        The number of "good" objects in the sample.
    N : int
        The number of objects in the population.
    cl : float in (0, 1)
        The desired confidence level.
    alternative : {"two-sided", "lower", "upper"}
        Indicates the alternative hypothesis. Defalut is two-sided
    G : int in [0, N]
        Starting point in search for confidence bounds for the hypergeometric parameter G.

    Returns
    -------
    ci_low, ci_upp : tuple
        lower and upper confidence level with coverage (at least)
        1-alpha.
    Notes
    -----
    xtol : float
        Tolerance
    rtol : float
        Tolerance
    maxiter : int
        Maximum number of iterations.
    """
    assert alternative in ("two-sided", "lower", "upper")
    assert 0 <= x <= n, 'impossible arguments'
    assert n <= N, 'impossible sample size'
    assert 0 < cl < 1, 'silly confidence level'

    if G is None:
        G = (x / n) * N
    ci_low = 0
    ci_upp = N

    if alternative == "lower" and x > 0:
        ci_low = x
        tail = hypergeom.sf(x - 1, N, ci_low, n)
        while tail < (1 - cl):
            ci_low += 1
            tail = hypergeom.sf(x - 1, N, ci_low, n)

    if alternative == "upper" and x < n:
        tail = hypergeom.sf(x, N, ci_upp, n)
        while tail > cl:
            ci_upp -= 1
            tail = hypergeom.sf(x, N, ci_upp, n)

    if alternative == 'two-sided':
        cl = 1 - (1 - cl) / 2
        eps = 0.1
        if x > 0:
            while x not in hypergeom_accept(N, ci_low, n, cl, randomized=False):
                ci_low += eps
            ci_low -= eps
            ci_low = round(ci_low)
        if x < n:
            while x not in hypergeom_accept(N, ci_upp, n, cl, randomized=False):
                ci_upp -= eps
            ci_upp += eps
            ci_upp = round(ci_upp)

    return ci_low, ci_upp