import pytest
from numpy.testing import assert_allclose

from ..fit_gaussian_with_class import (values_list_from_result_lines, 
                                       n_lines_from_values_list)


def test_values_list_from_result_lines():

    # arrange
    raw_strings = ['  0.6056423       3.998480       12.70451       0.000000       0.000000    ', 
                   '   0.000000       0.000000       0.000000       0.000000       0.000000    ', 
                   '   0.000000       0.000000       0.000000       0.000000       0.000000    ']

    expected_values = [0.6056423, 3.998480, 12.70451, 0., 0.,
                       0., 0., 0., 0., 0., 
                       0., 0., 0., 0., 0.] 

    # act
    actual_values = values_list_from_result_lines(raw_strings)

    # assert
    assert_allclose(actual_values, expected_values)


def test_n_lines_from_values_list():

    # arrange
    values_list = [0.6056423, 3.998480, 12.70451, 0., 0.,
                   0., 0., 0., 0., 0., 
                   0., 0., 0., 0., 0.]
    expected_n_lines = 1

    # act
    n_lines = n_lines_from_values_list(values_list)

    # assert
    assert n_lines == expected_n_lines