#!/usr/bin/env python3
"""Tests for mdinterface.utils.auxiliary utility functions."""

import numpy as np
import pytest
import ase

from mdinterface.utils.auxiliary import (
    as_list,
    mass2symbol,
    same_rev_check,
    find_smallest_missing,
    remove_inverted_tuples,
    round_list_to_sum,
    atoms_to_indexes,
)


class TestAsList:
    def test_none_returns_empty(self):
        assert as_list(None) == []

    def test_int_returns_single_item_list(self):
        assert as_list(5) == [5]

    def test_numpy_int_returns_single_item_list(self):
        assert as_list(np.int64(3)) == [3]

    def test_list_passthrough(self):
        assert as_list([1, 2, 3]) == [1, 2, 3]

    def test_tuple_converted_to_list(self):
        assert as_list((1, 2)) == [1, 2]

    def test_numpy_array_converted_to_list(self):
        result = as_list(np.array([4, 5, 6]))
        assert result == [4, 5, 6]

    def test_string_wrapped_in_list(self):
        # Strings are not iterable for our purposes
        assert as_list("H") == ["H"]

    def test_float_wrapped_in_list(self):
        assert as_list(1.5) == [1.5]


class TestMass2Symbol:
    def test_exact_match_oxygen(self):
        result = mass2symbol(15.999, ["O", "N", "C"])
        assert result == "O"

    def test_exact_match_carbon(self):
        result = mass2symbol(12.011, ["O", "N", "C"])
        assert result == "C"

    def test_within_tolerance(self):
        result = mass2symbol(15.95, ["O", "N", "C"], tolerance=0.1)
        assert result == "O"

    def test_no_match_raises_valueerror(self):
        # Mass 999 is far from any real element in the list
        with pytest.raises(ValueError):
            mass2symbol(999.0, ["H", "O"], tolerance=0.1)


class TestSameRevCheck:
    def test_identical_lists(self):
        assert same_rev_check([1, 2, 3], [1, 2, 3]) is True

    def test_reversed_lists(self):
        assert same_rev_check([1, 2, 3], [3, 2, 1]) is True

    def test_different_lists(self):
        assert same_rev_check([1, 2, 3], [1, 3, 2]) is False

    def test_single_element(self):
        assert same_rev_check([42], [42]) is True


class TestFindSmallestMissing:
    def test_missing_zero(self):
        assert find_smallest_missing([1, 2, 3]) == 0

    def test_missing_in_middle(self):
        assert find_smallest_missing([0, 1, 3, 4]) == 2

    def test_empty_list(self):
        assert find_smallest_missing([]) == 0

    def test_with_start_offset(self):
        assert find_smallest_missing([1, 2, 3], start=1) == 4

    def test_consecutive_from_zero(self):
        assert find_smallest_missing([0, 1, 2, 3]) == 4


class TestRemoveInvertedTuples:
    def test_removes_reverse(self):
        # Iterates backwards: keeps the last-seen of each pair and removes the earlier one
        data = [(1, 2), (2, 1), (3, 4)]
        remove_inverted_tuples(data)
        assert len(data) == 2
        assert (3, 4) in data
        # Exactly one of (1, 2) / (2, 1) remains
        assert (1, 2) in data or (2, 1) in data

    def test_no_duplicates(self):
        data = [(1, 2), (3, 4)]
        remove_inverted_tuples(data)
        assert len(data) == 2

    def test_empty_list(self):
        data = []
        remove_inverted_tuples(data)
        assert data == []


class TestRoundListToSum:
    def test_sums_to_target(self):
        lst = [0.1, 0.2, 0.3]
        result = round_list_to_sum(lst, 1.0)
        assert abs(sum(result) - 1.0) < 1e-9

    def test_already_correct(self):
        lst = [0.333, 0.333, 0.334]
        result = round_list_to_sum(lst, 1.0)
        assert abs(sum(result) - 1.0) < 1e-9

    def test_length_preserved(self):
        lst = [1.0, 2.0, 3.0]
        result = round_list_to_sum(lst, 6.0)
        assert len(result) == 3


class TestAtomsToIndexes:
    def setup_method(self):
        self.atoms = ase.Atoms("H2O", positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0)])

    def test_by_symbol(self):
        idxs = atoms_to_indexes(self.atoms, "O")
        assert idxs == [2]

    def test_by_symbol_multiple(self):
        idxs = atoms_to_indexes(self.atoms, "H")
        assert idxs == [0, 1]

    def test_by_index(self):
        idxs = atoms_to_indexes(self.atoms, 0)
        assert idxs == [0]

    def test_all(self):
        idxs = atoms_to_indexes(self.atoms, "all")
        assert idxs == [0, 1, 2]

    def test_mixed_list(self):
        idxs = atoms_to_indexes(self.atoms, [0, "O"])
        assert set(idxs) == {0, 2}
