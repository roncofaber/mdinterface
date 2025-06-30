#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 12:01:53 2024

@author: roncofaber
"""

# scientific libraries
import numpy as np
import scipy.integrate
import scipy.interpolate as sint
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.integrate import simpson
from scipy.integrate import solve_bvp

# internal modules
from mdinterface.utils.units import eps00

import copy

#%%

def integrate_poisson_1D(r, charge_dens_profile, periodic=False,
                         return_efield=False, rule="trapezoid"):
    """Integrate Poisson's equation twice to get electrostatic potential from charge density.

    Inputs:
    r : 1D numpy array. Values of coordinates, in Angstroms.
    charge_dens_profile : 1D numpy array. Values of charge density, 
        in e*Angstrom^-3.
    periodic : Boolean. Optional keyword. If True, adds linear function to 
        ensure values at boundaries are periodic.
    Outputs:
    phi : Numpy array of same length as r. Electrostatic potential in V. 
        The value of phi at the first r value is set to be 0.
    """
    
    if rule == "trapezoid":
        integrator = scipy.integrate.cumulative_trapezoid
    elif rule == "simpson":
        integrator = scipy.integrate.cumulative_simpson
    
    Efield =  integrator(charge_dens_profile, x=r, initial=0)
    Vfield = -integrator(Efield, x=r, initial=0)
    
    if periodic:
        # Ensure periodic boundary conditions by adding a linear function such
        # that the values at the boundaries are equal
        Vfield = Vfield - ((Vfield[-1]-Vfield[0])/(r[-1]-r[0])*(r-r[0]) + Vfield[0])
    
    if return_efield:
        return Vfield/eps00, Efield/eps00
    else:
        return Vfield/eps00

def integrate_poisson_ND(coords, charge_dens_profile, periodic=True):
    """Integrate Poisson's equation in N dimensions to get electrostatic potential.

    Inputs:
        coords: List of 1D NumPy arrays representing coordinates along each dimension.
        charge_dens_profile: N-dimensional NumPy array of charge density.
        periodic: Boolean (optional, default True). Enforce periodic boundary conditions.

    Outputs:
        phi: N-dimensional NumPy array of electrostatic potential (V).
    """

    phi = charge_dens_profile.copy()  # Start with charge density

    for axis, coord in enumerate(coords):
        phi = -simpson(phi, coord, axis=axis) / eps00  # Integrate along each axis

    if periodic:
        phi -= np.mean(phi)  # Center potential for periodic systems

    return phi


class Poisson1D_bvp:
    """
    A class to solve the 1D Poisson equation using boundary value problem (BVP) approach.
    
    Attributes:
    ----------
    _solvent : object
        An object representing the solvent with methods eps(z) and epsp(z).
    _z_grid : numpy.ndarray
        The grid points along the z-axis.
    """
    
    def __init__(self, z_grid, solvent):
        """
        Initializes the Poisson1D_bvp solver with the given grid and solvent properties.
        
        Parameters:
        ----------
        z_grid : numpy.ndarray
            The grid points along the z-axis.
        solvent : object
            An object representing the solvent with methods eps(z) and epsp(z).
        """
        self._solvent = solvent
        self._z_grid  = z_grid 
        
    def solve(self, rho_at_z, phi, pd_iter):
        """
        Solves the 1D Poisson equation for the given charge density.
        
        Parameters:
        ----------
        rho_at_z : numpy.ndarray
            Charge density at each grid point.
        phi : numpy.ndarray
            Initial guess for the potential.
        pd_iter : float
            Boundary condition at the first grid point.
        
        Returns:
        -------
        numpy.ndarray
            The potential and its gradient at each grid point.
        """
        # Interpolate rho as a function
        rho = sint.interp1d(self._z_grid, rho_at_z, kind='cubic', fill_value="extrapolate")
            
        def poisson_equation(z, potential):
            """
            Computes the derivative for the Poisson equation in 1D.
            
            Parameters:
            z (float): Position variable.
            potential (array): Array containing [phi, phi'] where phi is the potential and phi' is its first derivative.
            
            Returns:
            array: The derivatives [phi', phi''].
            """
            eps_z  = self._solvent.eps(z)
            phi_p  = potential[1]
            phi_pp = -(self._solvent.epsp(z) / eps_z) * potential[1] - rho(z) / eps_z
            
            return np.vstack((phi_p, phi_pp))
      
        def boundary_conditions(phi_at_start, phi_at_end):
            """
            Defines the boundary conditions for the BVP solver.
            
            Parameters:
            phi_at_start (array): Potential and its derivative at the start of the grid.
            phi_at_end (array): Potential and its derivative at the end of the grid.
            
            Returns:
            array: The boundary conditions.
            """
            return np.array([phi_at_start[0] - pd_iter, phi_at_end[0]])
      
        # Solve for phi using solve_bvp
        sol_ode = solve_bvp(poisson_equation, boundary_conditions, self._z_grid, phi)
    
        if sol_ode.status != 0:
            raise RuntimeError("BVP solver did not converge.")
        
        return sol_ode.sol(self._z_grid)
    
    
class Poisson1D_fd:
    """
    Class to solve the 1D Poisson equation using finite differences.
    """
    
    def __init__(self, z_grid, eps, epsp):
        """
        Initialize the Poisson1D_fd solver.

        Parameters:
        z_grid : array-like
            Grid points along the z-axis.
        eps : callable, dielectric function
        epsp: callable, derivative of dielectric function
            Solvent object with methods to get permittivity and its derivative.
        """
        
        h = (z_grid[1] - z_grid[0])
        N = len(z_grid)
        
        # Pre-calculate terms on grid
        eps_d  = eps(z_grid)
        epsp_d = epsp(z_grid)

        # Main diagonal
        main_diag     = np.full(N, -4 * eps_d)
        main_diag[0]  = 1
        main_diag[-1] = 1

        # Off-diagonal terms
        upper_diag = np.zeros(N-1)
        lower_diag = np.zeros(N-1)
        upper_diag[1:N-1] = 2 * eps_d[1:N-1] + h * epsp_d[1:N-1]
        lower_diag[0:N-2] = 2 * eps_d[1:N-1] - h * epsp_d[1:N-1]

        # Sparse matrix creation
        diagonals = [lower_diag, main_diag, upper_diag]
        A = diags(diagonals, offsets=[-1, 0, 1], shape=(N, N), format='csr')
        
        self._A = A
        self._N = N
        self._h = h
        self._z_grid = z_grid
    
    def solve(self, rho_at_z, phi, Vleft, Vright=0):
        """
        Solve the Poisson equation for the given charge density.

        Parameters:
        rho_at_z : array-like
            Charge density profile.
        phi : array-like
            Initial guess for the electrostatic potential.
        pd_iter : float
            Parameter for the boundary condition at the first grid point.

        Returns:
        array
            Electrostatic potential and its gradient.
        """
        
        # Define the right-hand side vector
        b     = -(2 * self._h**2) * rho_at_z
        b[0]  = Vleft
        b[-1] = Vright
        
        # Solve the system
        phi = spsolve(self._A, b)
        
        # Calculate the gradient of the potential
        phi_gradient = np.gradient(phi, self._z_grid, edge_order=2)
    
        return np.array([phi, phi_gradient])
