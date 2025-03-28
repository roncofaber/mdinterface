#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 14:04:43 2024

@author: roncofaber
"""

import re
import collections
import numpy as np

import ase
from ase.data import atomic_masses, chemical_symbols
#%%

def mass2symbol(mass, possible_symbols, tolerance=0.1):
    """
    Guess the atomic symbol from the given mass, restricted to a set of possible symbols.

    Parameters:
    mass (float): The atomic mass to guess the symbol for.
    possible_symbols (list): A list of possible atomic symbols to consider.
    tolerance (float): The tolerance within which to consider a match (default is 0.1).

    Returns:
    str: The guessed atomic symbol, or None if no match is found within the tolerance.
    """
    closest_mass = None
    closest_symbol = None
    min_diff = float('inf')

    for symbol in possible_symbols:
        index = chemical_symbols.index(symbol)
        atomic_mass = atomic_masses[index]
        diff = abs(mass - atomic_mass)
        if diff < min_diff:
            min_diff = diff
            closest_mass = atomic_mass
            closest_symbol = symbol

    if min_diff <= tolerance:
        return closest_symbol
    else:
        raise "Could NOT guess atom type."
    
def label_to_element(atostr, atomss):
    """
    Attempts to determine the chemical element symbol corresponding to a given 
    string and atomic mass.

    Args:
        atostr (str): A string potentially representing a chemical element.
        atomss (float): The approximate atomic mass of the element.

    Returns:
        str: The inferred chemical element symbol.

    Raises:
        ValueError: If the function cannot determine a valid element from the input.
    """


    new_label = re.sub(r'[^A-Za-z]', '', atostr).capitalize()    # Clean up input
    
    # load ase info
    atomic_masses = ase.data.atomic_masses
    elements = ase.data.chemical_symbols
    
    # initialize variables
    is_ready = False  # Flag to track if the element is found
    tried_last_resort = False 

    while not is_ready:
        try:
            try_atom = ase.Atom(new_label)  # Attempt to create an Atom object
            existent = True
        except:
            existent = False

        if existent and np.abs(try_atom.mass - atomss) < 1:  # Check mass match
            is_ready = True

        if not is_ready:
            new_label = new_label[:-1]  # Shorten the label for the next attempt

            if not new_label:  # If the label is empty, try a last-resort approach
                
                new_label = elements[np.argmin(np.abs(atomss - atomic_masses))]

                if tried_last_resort:
                    raise ValueError("{} is not a valid element".format(new_label))

                tried_last_resort = True

    return new_label


# return copy of input as list if not one
def as_list(inp):
    if inp is None:
        return []
    elif isinstance(inp, int) or isinstance(inp, np.int64):
        return [inp]
    elif isinstance(inp, collections.abc.Iterable) and not isinstance(inp, str): 
        # Handles lists, tuples, NumPy arrays, etc. (Excludes strings)
        return list(inp)  
    else:
        return [inp] # prone to error?
        # raise TypeError(f"Cannot convert type {type(inp)} to list")
        

def find_smallest_missing(data, start=0):
    """Finds the next smallest integer that is not in the list.
  
    This function efficiently finds the next smallest integer that is not present in the input list.
    It leverages sets for fast membership checks.
    
    Args:
        data: A list of integers.
  
    Returns:
        The next smallest integer that is not in the list.
      """

    data_set = set(data)  # Convert the list to a set for efficient membership checks
    smallest = start          # Start with the smallest possible positive integer
    while smallest in data_set:  # Check if 'smallest' is in the set
        smallest += 1              # If found, increment 'smallest'
    return smallest            # Return the first integer not found in the set 


def remove_inverted_tuples(list_of_tuples):
    seen = set()
    for i in range(len(list_of_tuples) - 1, -1, -1):  # Iterate backwards
        tup = list_of_tuples[i]
        reversed_tup = tup[::-1]
        if reversed_tup in seen:
            del list_of_tuples[i]  # Remove the tuple
        else:
            seen.add(tup)
            
# return list of indexes from mixed input of indexes and string (elements)
def atoms_to_indexes(system, symbols):

    # check if symbols is a list of strings
    if isinstance(symbols, str):
        if symbols == 'all':
            return list(range(len(system.get_chemical_symbols())))

    symbols = as_list(symbols)

    indexes = []
    for symbol in symbols:
        if not isinstance(symbol, str):
            indexes.append(symbol)
        else:
            for cc, atom in enumerate(system.get_chemical_symbols()):
                if atom == symbol:
                    indexes.append(cc)
    return np.unique(indexes).tolist()

# chunk a sequence in bits of approx. same size
def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

# check if two lists are the same even if reversed
def same_rev_check(list1, list2):
    if list(list1) == list(list2):
        return True
    elif list(list1) == list(reversed(list2)):
        return True
    return False


def round_list_to_sum(lst, target_sum, decimals=3):
    rounded_list = [round(x, decimals) for x in lst]
    current_sum = sum(rounded_list)
    difference = target_sum - current_sum
    
    # Sort the list by the decimal part to minimize the adjustment impact
    sorted_indices = sorted(range(len(lst)), key=lambda i: lst[i] - rounded_list[i], reverse=True)
    
    for i in sorted_indices:
        if abs(difference) < 10**(-(decimals+1)):
            break
        adjustment = round(difference, decimals)
        new_value = round(rounded_list[i] + adjustment, decimals)
        if new_value != rounded_list[i]:
            rounded_list[i] = new_value
            difference -= adjustment
    
    return rounded_list
