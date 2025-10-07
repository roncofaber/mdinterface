#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 17:17:30 2024

@author: roncofaber
"""

# numpy or nothing
import numpy as np

# ase stuff
import ase
import ase.build
from ase.data import atomic_numbers, covalent_radii

#%%

def calculate_overlap_score(oligomer, monomer, cutoff=5.0):
    """
    Calculate overlap score between oligomer and monomer based on interatomic distances.
    Lower score indicates less overlap (better positioning).
    """
    oli_pos = oligomer.get_positions()
    mon_pos = monomer.get_positions()

    # Calculate all pairwise distances
    distances = np.linalg.norm(oli_pos[:, np.newaxis, :] - mon_pos[np.newaxis, :, :], axis=2)

    # Calculate weighted overlap score (closer atoms contribute more)
    overlap_score = np.sum(1.0 / (distances + 1e-8))

    return overlap_score

def calculate_overlap_gradient(oligomer, monomer, rotation_center):
    """
    Calculate the gradient of overlap score with respect to rotations around x,y,z axes.
    Returns torque vector that indicates direction to rotate for minimal overlap.
    """
    oli_pos = oligomer.get_positions()
    mon_pos = monomer.get_positions()

    # Center monomer positions relative to rotation center
    mon_centered = mon_pos - rotation_center

    # Calculate pairwise distance vectors and distances
    diff_vectors = oli_pos[:, np.newaxis, :] - mon_pos[np.newaxis, :, :]  # (n_oli, n_mon, 3)
    distances = np.linalg.norm(diff_vectors, axis=2) + 1e-8  # (n_oli, n_mon)

    # Gradient of 1/r with respect to monomer positions: -1/r^2 * (r_vec/|r_vec|)
    force_magnitude = -1.0 / (distances**2)  # (n_oli, n_mon)
    unit_vectors = diff_vectors / distances[:, :, np.newaxis]  # (n_oli, n_mon, 3)

    # Force on each monomer atom from all oligomer atoms
    forces_on_monomer = np.sum(force_magnitude[:, :, np.newaxis] * unit_vectors, axis=0)  # (n_mon, 3)

    # Convert forces to torques around rotation center: τ = r × F
    torques = np.cross(mon_centered, forces_on_monomer)  # (n_mon, 3)

    # Sum all torques to get total gradient
    total_gradient = np.sum(torques, axis=0)  # (3,)

    return total_gradient

def generate_random_normalized_vector():
    # Generate random x, y, z coordinates
    x, y, z = np.random.rand(3) * 2 - 1  # Generate values from -1 to 1

    # Create a vector with these components
    vector = np.array([x, y, z])

    # Normalize the vector
    norm = np.linalg.norm(vector)  # Compute the magnitude
    if norm == 0:  # To handle the edge case where vector length is 0
        return np.array([1, 0, 0])  # Return a default unit vector
    normalized_vector = vector / norm  # Divide by norm to normalize

    return normalized_vector

def generate_spherical_points(n_samples):
    """
    Generate uniformly distributed points on the surface of a sphere.

    Parameters
    ----------
    n_samples : int
        The number of points to generate.

    Returns
    -------
    np.array
        An array of shape (n_samples, 3) representing points on the sphere.
    """
    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians
    
    for i in range(n_samples):
        y = 1 - (i / float(n_samples - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y
        
        theta = phi * i  # golden angle increment
        
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        
        points.append([x, y, z])
    
    return np.array(points)

def optimize_monomer_rotation_gradient(oligomer, monomer, x0, v1, max_iterations=50, learning_rate=0.1):
    """
    Optimize monomer rotation using gradient descent on overlap function.
    """
    best_monomer = monomer.copy()
    best_score = calculate_overlap_score(oligomer, best_monomer)

    current_monomer = monomer.copy()

    for iteration in range(max_iterations):
        # Calculate gradient (torque direction)
        gradient = calculate_overlap_gradient(oligomer, current_monomer, x0)
        gradient_norm = np.linalg.norm(gradient)

        if gradient_norm < 1e-6:  # Converged
            break

        # Normalize and apply small rotation in gradient direction
        rotation_axis = gradient / gradient_norm
        rotation_angle = learning_rate * gradient_norm

        # Limit rotation angle to prevent overshooting
        rotation_angle = min(rotation_angle, 0.2)  # Max ~11 degrees per step

        test_monomer = current_monomer.copy()
        test_monomer.rotate(rotation_axis, rotation_angle, center=x0)

        test_score = calculate_overlap_score(oligomer, test_monomer)

        if test_score < best_score:
            best_score = test_score
            best_monomer = test_monomer.copy()
            current_monomer = test_monomer.copy()
        else:
            # Reduce learning rate if no improvement
            learning_rate *= 0.8
            if learning_rate < 1e-3:
                break

    return best_monomer

def optimize_monomer_rotation(oligomer, monomer, x0, v1, n_samples=3600):
    """
    Rotate monomer around the attachment axis to minimize overlap with oligomer.

    Parameters:
    -----------
    oligomer : ase.Atoms
        The existing oligomer structure
    monomer : ase.Atoms
        The monomer to be attached (already translated to attachment position)
    x0 : np.array
        Position where attachment occurs (the center of rotation)
    v1 : np.array
        The vector indicating the preferred rotation direction
    n_samples : int
        Number of rotation angles to sample (default is 1000)

    Returns:
    --------
    ase.Atoms : Optimally rotated monomer
    """

    best_score = float('inf')
    best_monomer = monomer.copy()

    # Generate points on the sphere
    spherical_points = generate_spherical_points(n_samples)

    for v_rand in spherical_points:
        # Copy the monomer to apply the rotation
        test_monomer = monomer.copy()

        # Apply rotation around attachment point
        test_monomer.rotate(v1, v_rand, center=x0)

        # Calculate overlap score
        score = calculate_overlap_score(oligomer, test_monomer)

        if score < best_score:
            best_score = score
            best_monomer = test_monomer.copy()

    return best_monomer

def optimize_monomer_rotation_mc(oligomer, monomer, x0, v1, n_samples=3600):
    """
    Rotate monomer around the attachment axis to minimize overlap with oligomer.

    Parameters:
    -----------
    oligomer : ase.Atoms
        The existing oligomer structure
    monomer : ase.Atoms
        The monomer to be attached (already translated to attachment position)
    attachment_point : np.array
        Position where attachment occurs
    monomer_attachment_point : int
        Index of the attachment atom in the monomer
    n_samples : int
        Number of rotation angles to sample (default 36 = 10° increments)

    Returns:
    --------
    ase.Atoms : Optimally rotated monomer
    """

    best_score = float('inf')
    best_monomer = monomer.copy()

    for trial in range(n_samples):
       
        v_rand = generate_random_normalized_vector()

        # Apply rotation around attachment point
        test_monomer = monomer.copy()

        test_monomer.rotate(v1, v_rand, center=x0)

        # Calculate overlap score
        score = calculate_overlap_score(oligomer, test_monomer)

        if score < best_score:
            best_score = score
            best_monomer = test_monomer.copy()

    return best_monomer

def start_chain(monomer):

    monomer = monomer.copy()

    # find starting index
    mon_idx = int(np.where(monomer.arrays["polymerize"] == 1)[0])

    # Find atom connected to X (the bonding atom)
    end_idx = ase.build.connected_indices(monomer, mon_idx)[1]

    # remember connecting points (mark the bonding atom, not the X atom)
    is_connected = np.array(len(monomer)*[False])
    # is_connected[end_idx] = True
    monomer.new_array("is_connected", is_connected)

    # mark monomer index
    monomer_label = 0
    monomer.new_array("mon_id", np.array(len(monomer)*[monomer_label]))

    return monomer

def attach_to_chain(oligomer, monomer):

    # copy to prevent changing
    oligomer = oligomer.copy()
    monomer  = monomer.copy()

    # Find indexes of atoms that will be polymerized (H atoms)
    oli_idx = int(np.where(oligomer.arrays["polymerize"] == 2)[0])
    mon_idx = int(np.where(monomer.arrays["polymerize"] == 1)[0])

    # Find atoms connected to X (most likely C atoms)
    ini_idx = ase.build.connected_indices(oligomer, oli_idx)[1]
    end_idx = ase.build.connected_indices(monomer, mon_idx)[1]

    # Find bond distance #TODO here we can change the distance if needed
    d1 = covalent_radii[atomic_numbers[oligomer.get_chemical_symbols()[ini_idx]]]
    d2 = covalent_radii[atomic_numbers[monomer.get_chemical_symbols()[end_idx]]]
    target_distance = (d1 + d2)/2
    
    oligomer.set_distance(ini_idx, oli_idx, target_distance, fix=0)
    monomer.set_distance(end_idx, mon_idx, target_distance, fix=0)
    
    # remember connecting points
    is_connected = np.array(len(monomer)*[False])
    is_connected[end_idx] = True
    monomer.new_array("is_connected", is_connected)
    # oligomer.arrays["is_connected"][ini_idx] = True

    # mark monomer index
    monomer_label = int(np.max(oligomer.arrays["mon_id"]+1))
    monomer.new_array("mon_id", np.array(len(monomer)*[monomer_label]))

    # Translate monomer so its bonding atom reaches the target position
    x0 = oligomer.get_positions()[oli_idx]
    x1 = monomer.get_positions()[mon_idx]
    monomer.translate(x0-x1)

    # Optimize monomer rotation to minimize overlap
    v1 = monomer.get_distance(mon_idx, end_idx, vector=True)
    monomer = optimize_monomer_rotation(oligomer, monomer, x0, v1, n_samples=3600)
    monomer = optimize_monomer_rotation_gradient(oligomer, monomer, x0, v1)

    # remove fluff
    del oligomer[oli_idx]
    del monomer[mon_idx]

    # return monomer attached to bigger object
    return oligomer + monomer


def build_polymer(monomers, sequence=None, nrep=None):
    
    if not isinstance(monomers, list):
        monomers = [monomers]
    
    if sequence is None:
        sequence = [0]*nrep
    
    for cc, seq in enumerate(sequence):
        
        if cc == 0:
            oligomer = start_chain(monomers[seq])
        else:
            oligomer = attach_to_chain(oligomer, monomers[seq])
    
    return oligomer
