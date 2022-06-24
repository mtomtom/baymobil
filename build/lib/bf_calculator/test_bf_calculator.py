import pytest
import bf_calculator as bf
import numpy as np
## Function to test our stirling binom function
## If one of the values is 0, then this function should return np.log(1)
def test_stirling_binom1():
    assert bf.stirling_binom(0,0) == np.log(1)

def test_stirling_binom2():
    assert bf.stirling_binom(100,0) == np.log(1)

def test_stirling_binom3():
    assert bf.stirling_binom(0,100) == np.log(1)

## If N = n, then the function returns np.log(1)
def test_stirling_binom3_equal_values():
    assert bf.stirling_binom(100,100) == np.log(1)

## Check that the two methods give similar results
def test_fasterpostN2_stirling_approx():
    test_output = bf.fasterpostN2_stirling(100,10,100,10,100,10,10)
    assert test_output[2] == pytest.approx(-0.2988, 0.01)

def test_fasterpostN2_approx():
    test_output = bf.fasterpostN2(100,10,100,10,100,10,10)
    assert test_output[2] == pytest.approx(-0.2988, 0.01)