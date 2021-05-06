from confidence import *



def test_tau_twosided_ci_1():
    """tau_twosided_ci test
    
    Tests legal inputs to test_tau_twosided_ci
    function and checks for proper output.
    Uses examples from Method 3 in Li & Ding's paper
    """
    
    pass


def test_tau_twosided_ci_2():
    """tau_twosided_ci test
    
    Tests legal inputs to test_tau_twosided_ci
    function and checks for proper output.
    """
    
    pass

def test_filterTable_1():
    """filterTable test
    
    Tests legal inputs to filterTable
    function and checks for proper output.
    """
    
    pass

def test_potential_outcomes_1():
    """potential_outcomes test
    
    Tests legal inputs to potential_outcomes
    function and checks for proper output.
    """
    
    pass
    
    
def test_N_generator_badinput_1():
    """N generator test
    
    Checks for invalid number of subjects"""
    pytest.raises(ValueError, N_generator, 5, 10, 10, 5, 5)
    
def test_N_generator_badinput_2():
    """N generator test
    
    Checks for negative subjects"""
    pytest.raises(ValueError, N_generator, -5, 10, 10, 5, 5)
    
    
def test_filterTable_badinput_1():
    """filterTable test
    
    Checks for invalid number of subjects"""
    pytest.raises(ValueError, filterTable, 5, 10, 10, 5, 5)
    
def test_filterTable_badinput_2():
    """filterTable test
    
    Checks for negative subjects"""
    pytest.raises(ValueError, filterTable, -5, 10, 10, 5, 5)
    
    
def test_tau_twosided_ci_badinput_1():
    """tau_twosided_ci test
    
    Asserts that a value error is raised
    if alpha is not in range (0,1)"""
    pytest.raises(ValueError, tau_twosided_ci, 5, 10, 10, 5, 1.1)

def test_tau_twosided_ci_badinput_2():
    """tau_twosided_ci test
    
    Checks for negative subjects"""
    pytest.raises(ValueError, tau_twosided_ci, -5, 10, 10, 5, 0.05)

def test_potential_outcomes_badinput_1():
    """potential_outcomes test
    
    Checks for input list containing too many numbers
    """
    
    pytest.raises(ValueError, potential_outcomes, [1,2,3,4,5])
    
def test_potential_outcomes_badinput_2():
    """potential_outcomes test
    
    Checks for input list containing negative numbers
    """
    
    pytest.raises(ValueError, potential_outcomes, [1,2,3,-4])