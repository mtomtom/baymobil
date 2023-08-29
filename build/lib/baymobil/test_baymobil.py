import pytest
import baymobil as baymob
import numpy as np

import hypothesis.strategies as st
from hypothesis import given


@given(fields=st.lists(st.text(), min_size=1, max_size=10))
@example([","])
def test_read_write_csv_hypothesis(fields):
    formatted_row = naive_write_csv_row(fields)
    parsed_row = naive_read_csv_row(formatted_row)
    assert fields == parsed_row


## Write our tests for the functions

## Main functions

## Simulations

## Plotting functions
