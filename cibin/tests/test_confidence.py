from confidence import *



def test_tau_twosided_ci_1():
    """Tests legal inputs to def test_tau_twosided_ci
    function and checks for proper output """
    
    

def test_tau_twosided_ci_badinput_1():
    """Asserts that a value error is raised
    if alpha is not in range (0,1)"""
    pytest.raises(ValueError, tau_twosided_ci, 5, 10, 10, 5, 1.1)

def test_tau_twosided_ci_badinput_2():
    pass

def test_tau_twosided_ci_badinput_3():
    pass

def test_tau_twosided_ci_badinput_4():
    pass