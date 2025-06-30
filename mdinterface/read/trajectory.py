#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 30 14:21:18 2025

@author: roncofaber
"""

# ASE stuff
import ase
import ase.io
from ase import units
import ase.visualize

import numpy as np

# package stuff
from mdinterface.utils.auxiliary import atoms_to_indexes, as_list
from mdinterface.utils.poisson import integrate_poisson_1D
from mdinterface.utils.units import au_per_ang3_to_g_per_cm3

#%%

common_metals = [
    "Fe", "Al", "Cu", "Au", "Ag", "Zn", "Ni", "Pb", "Sn", "Ti", "Mg", "Hg", "Pt",
    "Na", "K", "Ca", "Cr", "Mn", "Co", "Mo", "W", "V", "Zr", "Pd", "Cd", "In", 
    "Li", "Rb", "Cs", "Be", "Sr", "Ba", "Sc", "Y", "Hf", "Nb", "Ta", "Re", "Os",
    "Ir", "Rh", "Ru"
]

class Trajectory(object):
    
    def __init__(self, filename, datafile=None, index=":", parallel=False,
                 electrode=None, electrode_atoms=None, zaxis=2, every=1,
                 is_molecule=False):
        
        # setup internal variables
        self._filename  = filename
        self._datafile  = datafile
        self._zaxis     = zaxis
        self._electrode = electrode
        
        self._traj, self._data = self._read_trajectory(filename, datafile,
                                                       index, parallel=parallel,
                                                       every=every)
        
        if is_molecule:
            self._make_molecule(self._traj)
        
        # setup vars
        if electrode:
            self._setup_electrode(electrode_atoms)
        
        return
    
    # this is format specific
    @staticmethod
    def _read_trajectory(filename, datafile, index, parallel=False, every=1):
        pass

    def _setup_electrode(self, electrode_atoms):

        # guess atoms
        if electrode_atoms is None:
            first_atom = np.array(self[0].get_chemical_symbols())[0]
            if first_atom in common_metals:
                electrode_atoms = first_atom
            else:
                raise "Please suggest electrode atoms"
        
        # find electrode indexes
        idxs = np.array(atoms_to_indexes(self[0], electrode_atoms))
        
        if len(idxs) == 0:
            self._int_idxs = []
            self._end_idxs = []
            return
        
        self._int_idxs = idxs[np.where(self[0].get_positions(wrap=True)[idxs]
                                       [:,self._zaxis] < self.depth/2)[0]]
        self._end_idxs = idxs[np.where(self[0].get_positions(wrap=True)[idxs]
                                       [:,self._zaxis] > self.depth/2)[0]]
        
        return
    
    def get_density(self, return_mean=False):
        
        tot_mass = self[0].get_masses().sum()
        volumes  = np.array([frame.get_volume() for frame in self])
        
        if not return_mean:
            return au_per_ang3_to_g_per_cm3*tot_mass/volumes
        else:
            return au_per_ang3_to_g_per_cm3*tot_mass/volumes.mean()
    
    # get 1D density along specific axis. Works best if cell is constant
    def get_1D_density(self, axis=None, nbins=1000, symbols="all", mass=None,
                       molar=None, charge=None, return_mean=True,
                       zmin=None, zmax=None, background=None):
        
        if axis is None:
            axis = self._zaxis
        
        weights = None
        
        # check not trying to do too much
        assert [ii is not None for ii in [molar, mass, charge]].count(True) == 1
        
        # make bins
        if zmin is None:
            zmin  = 0  # Assuming your simulation cell starts at z = 0
        if zmax is None:
            zmax  = self[0].get_cell()[axis, axis]  # Get the z-dimension of the cell
        z_bins = np.linspace(zmin, zmax, nbins + 1)  
        
        # calculate bin cell
        # bc_volume = ((zmax-zmin)/self.cell[2][2])*self.cell.volume/(nbins + 1)
        bc_volume = self.cell[0][0]*self.cell[1][1]*(z_bins[1]-z_bins[0])
        # print(bc_volume)
        # print(v2)
        
        # mask elements
        mask = atoms_to_indexes(self.traj[0], symbols)
        
        # get relevant frames positions
        positions = [frame[mask].get_positions(wrap=True)[:,axis] for frame in self]
        
        if charge:
            weights = [frame[mask].get_initial_charges() for frame in self]
        
        # calculate histogram
        if return_mean:
            density, _ = np.histogram(positions, bins=z_bins, weights=weights)
            density = density/len(positions)
        else:
            if not charge:
                weights = len(positions)*[None]
            density = []
            for pos, chg in zip(positions, weights):
                density.append(np.histogram(pos, bins=z_bins, weights=chg)[0])
        
        # return in q/Ang^3
        if charge is not None:
            
            density = density/bc_volume
            
            if background is not None:
                bk_chg = background/self.cell.volume
                density -= bk_chg
            
        # return in mol/L
        if molar is not None:
            density = density/(molar*units.mol*1e-27*bc_volume)    
            
        # return mass density in g/cm3
        if mass is not None:
            density = mass*au_per_ang3_to_g_per_cm3*(density/bc_volume)
        
        return (z_bins[:-1] + z_bins[1:])/2, density
    
    # get 2D density
    def get_2D_density(self, axes=None, nbins=[100, 100], symbols="all", mass=None,
                       molar=None, charge=None, return_mean=True,
                       zmin=None, zmax=None, background=None, zrange=None):
        
        if axes is None: # if none, assume it's for the non z-axis (xy)
            axes = [ii for ii in [0, 1, 2] if ii != self._zaxis]
        
        zaxis = [ii for ii in [0, 1, 2] if ii not in axes]
        
        if len(as_list(nbins)) == 1:
            nbins = 2*as_list(nbins)
        
        # check not trying to do too much
        assert [ii is not None for ii in [molar, mass, charge]].count(True) == 1
        
        # make bins
        if zmin is None:
            zmin  = [0, 0]  # Assuming your simulation cell starts at z = 0
        if zmax is None:
            zmax  = [self[0].get_cell()[ax, ax] for ax in axes] # Get the z-dimension of the cell
            
        z_bins = [np.linspace(zmin[cc], zmax[cc], nbins[cc] + 1)  for cc in range(2)]
        
        # calculate bin cell
        if zrange is None:
            zlength = self[0].get_cell()[zaxis, zaxis]
        else:
            zlength = np.subtract(*zrange)
        
        bc_volume = abs(zlength*np.prod(np.subtract(zmax, zmin)))/np.prod(nbins)
        
        # mask elements
        mask = atoms_to_indexes(self.traj[0], symbols)
        
        # get relevant frames positions
        density = []
        for frame in self:
            
            pos = frame[mask].get_positions(wrap=True)
            
            if zrange is not None:
                zmask = np.where(np.logical_and(pos[:, zaxis]>=zrange[0],
                                                pos[:, zaxis]<=zrange[1]))[0]
            else:
                zmask = slice(None)
                
            if charge is None:
                density.append(np.histogram2d(pos[zmask, axes[0]], pos[zmask, axes[1]],
                                              bins=z_bins)[0])
            else:
                chg = frame[mask].get_initial_charges()
                density.append(np.histogram2d(pos[zmask, axes[0]], pos[zmask, axes[1]],
                                              bins=z_bins)[0], weights=chg[zmask])
                
        if return_mean:
            density = np.mean(density, axis=0)

        # return in q/Ang^3
        if charge is not None:
            
            density = density/bc_volume
            
            if background is not None:
                bk_chg = background/self.cell.volume
                density -= bk_chg
            
        # return in mol/L
        if molar is not None:
            density = density/(molar*units.mol*1e-27*bc_volume)    
            
        # return mass density in g/cm3
        if mass is not None:
            density = mass*density / (units.mol*1e-24*bc_volume)
        
        return (z_bins[0][:-1] + z_bins[0][1:])/2, (z_bins[1][:-1] + z_bins[1][1:])/2, density
    
    
    def get_1D_potential(self, axis=None, nbins=100, symbols="all",
                         return_mean=True, periodic=False, zmin=None, zmax=None,
                         correct_background=False, rule="trapezoid",
                         return_efield=False):
        
        zpos, chgden = self.get_1D_density(axis=axis, nbins=nbins, symbols=symbols,
                                           charge=True, return_mean=return_mean,
                                           zmin=zmin, zmax=zmax)
        
        if correct_background and not return_mean:
            bkg_chg = np.array([frame.get_initial_charges().sum() for frame in self])/self.cell.volume
            chgden  = (chgden.T-bkg_chg).T
        elif correct_background:
            bkg_chg = np.array([frame.get_initial_charges().sum() for frame in self])/self.cell.volume
            chgden = chgden - np.mean(bkg_chg)
            
        pot = integrate_poisson_1D(zpos, chgden, periodic=periodic, rule=rule,
                                   return_efield=return_efield)
        
        return zpos, pot
    
    def get_charges_on_electrodes(self, electrode_atoms=None):
        
        assert self._electrode
        
        int_chgs = []
        end_chgs = []
        for frame in self:
            charges = frame.get_initial_charges()
            int_chgs.append(charges[self._int_idxs].sum())
            end_chgs.append(charges[self._end_idxs].sum())
               
        return np.array(int_chgs), np.array(end_chgs)
    
    def get_potential_difference(self, start_coord, end_coord, nbins=2000,
                                 axis=None, periodic=False, correct_background=False):
        
        # get 1D potential
        zpos, pot = self.get_1D_potential(nbins=2000, symbols="all",
                                          return_mean=False, periodic=periodic,
                                          correct_background=correct_background,
                                          axis=axis)
        
        # find closest indexes
        st_idx = (np.abs(zpos - start_coord)).argmin()
        en_idx = (np.abs(zpos - end_coord)).argmin()

        dV = pot[:,en_idx] - pot[:,st_idx]
        
        return dV
    
    def get_dV_on_electrodes(self, electrode_atoms=None, nbins=2000,
                                 axis=None, periodic=False, correct_background=False):
        
        assert self._electrode
        
        if axis is None:
            axis = self._zaxis
 
        st_crd = np.mean(self[0][self._int_idxs].get_positions()[:,axis])
        en_crd = np.mean(self[0][self._end_idxs].get_positions()[:,axis])
        
        return self.get_potential_difference(st_crd, en_crd, nbins=nbins,
                                     axis=axis, periodic=periodic,
                                     correct_background=correct_background)
    
    def get_coordination_number(self, ref_atom, tar_atom, r0, frames=slice(None),
                                n=6, m=12, a=None):
        
        if a is not None:
            def fun_CN(dist, r0, n, m, a):
                return 1 - 1/(1 + np.exp(-a*(dist-r0)))
        else:
            def fun_CN(dist, r0, n, m, a):
                return ( 1 - (dist/r0)**n ) / (  1 - (dist/r0)**m )
           
        ref_idxs = atoms_to_indexes(self[0], ref_atom)
        tar_idxs = atoms_to_indexes(self[0], tar_atom)
        
        CN = []
        for frame in self[frames]:
            
            # distance = frame.get_all_distances(mic=True)
            
            frame_CN = []
            for ref_idx in ref_idxs:
                dist = frame.get_distances(ref_idx, tar_idxs, mic=True)
                # dist = distance[ref_idx][tar_idxs]
                
                cn = fun_CN(dist, r0, n, m, a)
                frame_CN.append(cn)
                
            CN.append(frame_CN)
        
        return np.squeeze(CN)
    
    def get_electrode_COM(self, average=True):
        
        COMs = [frame[self._int_idxs].get_center_of_mass()[2] for frame in self]
        
        if average:
            return np.mean(COMs)
        
        return np.array(COMs)
    
    @staticmethod
    def _make_molecule(traj):
        
        for frame in traj:
            frame.center()
            frame.set_pbc(False)
            frame.set_cell(None)
        
        return
    
    @property
    def depth(self):
        return self._traj[0].cell[self._zaxis, self._zaxis]
    
    @property
    def area(self):
        # Define the available axes
        all_axes = [0, 1, 2]
        
        # Filter out the _zaxis
        area_axes = [axis for axis in all_axes if axis != self._zaxis]
        
        # Assuming the area is calculated as the product of the two axes dimensions
        axis1, axis2 = area_axes
        return self._traj[0].cell[axis1, axis1] * self._traj[0].cell[axis2, axis2]
    
    @property
    def cell(self):
        return self._traj[0].cell
    
    @property
    def traj(self):
        return self._traj
    
    @property
    def data(self):
        return self._data
    
    def view(self):
        ase.visualize.view(self.traj)
        return
    
    # make class iterable
    def __getitem__(self, index):
        return self._traj[index]
    
    def __iter__(self):
        return iter(self._traj)
    
    def __len__(self):
        return len(self._traj)

    def __add__(self, other):
        if isinstance(other, Trajectory):
            new_traj = self.traj + other.traj
            self._traj = new_traj
            return self
        else:
            raise TypeError("Unsupported operand type for +")