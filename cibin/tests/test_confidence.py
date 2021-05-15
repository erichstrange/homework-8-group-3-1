"""Test Suite for Confidence Intervals."""

from cibin import *
import pytest


def test_tau_twosided_ci_1():
    """tau_twosided_ci test.

    Tests legal inputs to test_tau_twosided_ci
    function and checks for proper output.
    Uses examples from Method 3 in Li & Ding's paper
    """
    bounds, a, b = tau_twosided_ci(2, 6, 8, 0, 0.05)
    new_bounds = [i*(2+6+8) for i in bounds]
    assert((new_bounds == [-14.0, -5.0]))
    bounds_1, a_1, b_1 = tau_twosided_ci(1, 1, 1, 13, 0.05)
    new_bounds_1 = [i*(1+1+1+13) for i in bounds_1]
    assert(new_bounds_1 == [-1.0, 14.0])
    bounds_2, a_2, b_2 = tau_twosided_ci(6, 0, 11, 3, 0.05)
    new_bounds_2 = [i*(6+11+3) for i in bounds_2]
    assert((new_bounds_2 == [-4.0, 8.0]))


def test_tau_twosided_ci_2():
    """tau_twosided_ci test.

    Tests that with small numbers, the lower bound
    and upper bound are the same whether exact ==True
    or exact== False
    """
    bounds1, a, b = tau_twosided_ci(1, 1, 2, 0, 0.05, True)
    bounds2, c, d = tau_twosided_ci(1, 1, 2, 0, 0.05, False)
    assert(bounds1 == bounds2)
    assert(a == c)


def test_filterTable_1():
    """filterTable test.

    Tests legal inputs to filterTable
    function and checks for proper output.
    """
    assert(!filterTable([5, 10, 10, 5], 6, 11, 10, 3))


def test_potential_outcomes_1():
    """Test potential_outcomes.

    Tests legal inputs to potential_outcomes
    function and checks for proper output.
    """
    assert(potential_outcomes([1, 0, 1, 0]).tolist() == [[0, 0], [1, 0]])


def test_filterTable_badinput_2():
    """filterTable test.

    Checks for negative subjects
    """
    pytest.raises(ValueError, filterTable, [5, 10, 2, 1], 10, -1, 5, 4)


def test_tau_twosided_ci_badinput_1():
    """tau_twosided_ci test.

    Asserts that a value error is raised
    if alpha is not in range (0,1)
    """
    pytest.raises(ValueError, tau_twosided_ci, 5, 10, 10, 5, 1.1)


def test_tau_twosided_ci_badinput_2():
    """tau_twosided_ci test.

    Checks for negative subjects
    """
    pytest.raises(ValueError, tau_twosided_ci, -5, 10, 10, 5, 0.05)


def test_potential_outcomes_badinput_1():
    """Test for potential_outcomes.

    Checks for input list containing too many numbers
    """
    pytest.raises(ValueError, potential_outcomes, [1, 2, 3, 4, 5])


def test_potential_outcomes_badinput_2():
    """Test for potential_outcomes.

    Checks for input list containing negative numbers
    """
    pytest.raises(ValueError, potential_outcomes, [1, 2, 3, -4])
