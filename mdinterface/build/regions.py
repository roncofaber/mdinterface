#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Geometric regions for spatially constrained PACKMOL placement.

A Region describes a 3D shape usable to constrain where PACKMOL places
molecules, either "inside" or "outside" it. Regions are the building block
for spatially heterogeneous solvent layers (e.g. a gas pocket embedded in a
bulk solvent) - see SimCell.add_solvent's `regions` parameter.
"""

import math
from dataclasses import dataclass
from typing import Any, List, Optional, Tuple


class Region:
    """
    Base class for a 3D geometric region usable as a PACKMOL constraint.

    Subclasses implement volume(), packmol_line(), and bounding_box(). This
    interface is deliberately minimal so combinators (union/intersection/
    invert) can be added later as new subclasses without changing anything
    that already exists.
    """

    def volume(self) -> float:
        raise NotImplementedError

    def packmol_line(self, mode: str) -> str:
        raise NotImplementedError

    def bounding_box(self) -> Tuple[float, float, float, float, float, float]:
        raise NotImplementedError

    def fill(self, solvent: Optional[Any] = None, solute: Optional[List[Any]] = None,
              nsolute=None, density: Optional[float] = None, nsolvent=None,
              concentration: Optional[float] = None, conmodel: Optional[dict] = None,
              ratio: Optional[List[float]] = None) -> "FilledRegion":
        """
        Bundle this region with the content that should fill it.

        Parameters mirror SimCell.add_solvent's content parameters, minus
        zdim/xysize since the region defines its own extent.
        """
        return FilledRegion(
            region=self, solvent=solvent, solute=solute, nsolute=nsolute,
            density=density, nsolvent=nsolvent, concentration=concentration,
            conmodel=conmodel, ratio=ratio,
        )

    @staticmethod
    def _check_mode(mode: str) -> None:
        if mode not in ("inside", "outside"):
            raise ValueError(f"mode must be 'inside' or 'outside', got {mode!r}")


@dataclass
class FilledRegion:
    """A Region paired with the species/density spec that should fill it."""
    region: Region
    solvent: Optional[Any] = None
    solute: Optional[List[Any]] = None
    nsolute: Optional[Any] = None
    density: Optional[float] = None
    nsolvent: Optional[Any] = None
    concentration: Optional[float] = None
    conmodel: Optional[dict] = None
    ratio: Optional[List[float]] = None


class Sphere(Region):
    """A spherical region centered at `center` with the given `radius`."""

    def __init__(self, center: Tuple[float, float, float], radius: float):
        if radius <= 0:
            raise ValueError(f"radius must be positive, got {radius}")
        self.center = tuple(float(c) for c in center)
        self.radius = float(radius)

    def volume(self) -> float:
        return (4.0 / 3.0) * math.pi * self.radius ** 3

    def packmol_line(self, mode: str) -> str:
        self._check_mode(mode)
        x, y, z = self.center
        return f"{mode} sphere {x} {y} {z} {self.radius}"

    def bounding_box(self) -> Tuple[float, float, float, float, float, float]:
        x, y, z = self.center
        r = self.radius
        return (x - r, y - r, z - r, x + r, y + r, z + r)

    def __repr__(self):
        return f"Sphere(center={self.center}, radius={self.radius})"


class Box(Region):
    """An axis-aligned box region centered at `center` with the given `size`."""

    def __init__(self, center: Tuple[float, float, float], size: Tuple[float, float, float]):
        if any(s <= 0 for s in size):
            raise ValueError(f"size components must be positive, got {size}")
        self.center = tuple(float(c) for c in center)
        self.size = tuple(float(s) for s in size)

    @classmethod
    def from_bounds(cls, xmin: float, ymin: float, zmin: float,
                      xmax: float, ymax: float, zmax: float) -> "Box":
        center = ((xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2)
        size = (xmax - xmin, ymax - ymin, zmax - zmin)
        return cls(center=center, size=size)

    def volume(self) -> float:
        return self.size[0] * self.size[1] * self.size[2]

    def packmol_line(self, mode: str) -> str:
        self._check_mode(mode)
        xmin, ymin, zmin, xmax, ymax, zmax = self.bounding_box()
        return f"{mode} box {xmin} {ymin} {zmin} {xmax} {ymax} {zmax}"

    def bounding_box(self) -> Tuple[float, float, float, float, float, float]:
        cx, cy, cz = self.center
        sx, sy, sz = self.size
        return (cx - sx / 2, cy - sy / 2, cz - sz / 2, cx + sx / 2, cy + sy / 2, cz + sz / 2)

    def __repr__(self):
        return f"Box(center={self.center}, size={self.size})"


class Cylinder(Region):
    """
    A cylindrical region centered at `center`, with the given `radius` and
    `height`, extending along `axis`. `height` is centered on `center` along
    `axis` (the cylinder spans center[axis] - height/2 to center[axis] + height/2).
    """

    _AXES = {"x": 0, "y": 1, "z": 2}

    def __init__(self, center: Tuple[float, float, float], radius: float,
                  height: float, axis: str = "z"):
        if radius <= 0:
            raise ValueError(f"radius must be positive, got {radius}")
        if height <= 0:
            raise ValueError(f"height must be positive, got {height}")
        if axis not in self._AXES:
            raise ValueError(f"axis must be 'x', 'y', or 'z', got {axis!r}")
        self.center = tuple(float(c) for c in center)
        self.radius = float(radius)
        self.height = float(height)
        self.axis = axis

    def volume(self) -> float:
        return math.pi * self.radius ** 2 * self.height

    def packmol_line(self, mode: str) -> str:
        self._check_mode(mode)
        axis_idx = self._AXES[self.axis]
        base = list(self.center)
        base[axis_idx] -= self.height / 2
        direction = [0.0, 0.0, 0.0]
        direction[axis_idx] = 1.0
        x0, y0, z0 = base
        dx, dy, dz = direction
        return f"{mode} cylinder {x0} {y0} {z0} {dx} {dy} {dz} {self.radius} {self.height}"

    def bounding_box(self) -> Tuple[float, float, float, float, float, float]:
        axis_idx = self._AXES[self.axis]
        cx, cy, cz = self.center
        lo = [cx - self.radius, cy - self.radius, cz - self.radius]
        hi = [cx + self.radius, cy + self.radius, cz + self.radius]
        lo[axis_idx] = self.center[axis_idx] - self.height / 2
        hi[axis_idx] = self.center[axis_idx] + self.height / 2
        return tuple(lo + hi)

    def __repr__(self):
        return (f"Cylinder(center={self.center}, radius={self.radius}, "
                f"height={self.height}, axis={self.axis!r})")
