#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# simulationFFT.py
# Created: April 6th, 2018
# Python implementation of GAIA lab's MGSimulFFT.m available here: http://wp.unil.ch/gaia/downloads/
# Author: Noah Athens
import numpy as np
from numpy.fft import fftn, fftshift, ifftn
from numpy.random import uniform as rand

def transform_distribution(grid, new_distribution):
    """ Transforms grid to new distribution."""
    old_distribution = np.sort(grid.flatten())
    new_distribution = np.sort(np.random.choice(new_distribution, size = grid.size))
    d = dict(zip(old_distribution, new_distribution))
    return np.vectorize(d.get)(grid)

def simulFFT(nx, ny, nz, mu, sill, m, lx , ly, lz):
    """ Performs unconditional simulation with specified mean, variance,
    and correlation length.
    """
    if nz == 0: nz = 1 # 2D case
    xx, yy, zz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz))
    points = np.stack((xx.ravel(), yy.ravel(), zz.ravel())).T
    centroid = points.mean(axis=0)
    length = np.array([lx, ly, lz])
    h = np.linalg.norm((points - centroid) / length, axis = 1).reshape((ny, nx, nz))

    if m == 'Exponential':
        c = np.exp(-3*h) * sill
    elif m == 'Gaussian':
        c = np.exp(-3*h**2) * sill
    else:
        raise(Exception('For m enter either "Exponential" or "Gaussian"'))

    grid = fftn(fftshift(c)) / (nx * ny * nz)
    grid = np.abs(grid)
    grid[0, 0, 0] = 0 # reference level
    ran = np.sqrt(grid) * np.exp(1j * np.angle(fftn(rand(size=(ny, nx, nz)))))
    grid = np.real(ifftn(ran * nx * ny * nz))
    std = np.std(grid)
    if nx == 1 or ny == 1 or nz == 1: grid = np.squeeze(grid)
    return grid / std * np.sqrt(sill) + mu

def simulFFT_rotate(nx, ny, mu, sill, m, lx , ly, angle):
    """ Separate version of simulFFT for 2D rotation. """
    ox = nx
    oy = ny
    if angle >= 0:
        multiplier = max(nx, ny) * 1.45
    else:
        multiplier = min(nx, ny) * 1.45
    nx = int(np.ceil(multiplier))
    ny = int(np.ceil(multiplier))
    xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
    points = np.stack((xx.ravel(), yy.ravel())).T
    centroid = points.mean(axis=0)
    length = np.array([lx, ly])
    h = np.linalg.norm((points - centroid) / length, axis = 1).reshape((ny, nx))

    if m == 'Exponential':
        c = np.exp(-3*h) * sill
    elif m == 'Gaussian':
        c = np.exp(-3*h**2) * sill
    else:
        raise(Exception('For m enter either "Exponential" or "Gaussian"'))

    grid = fftn(fftshift(c)) / (nx * ny)
    grid = np.abs(grid)
    grid[0, 0] = 0 # reference level
    ran = np.sqrt(grid) * np.exp(1j * np.angle(fftn(rand(size=(ny, nx)))))
    grid = np.real(ifftn(ran * nx * ny))
    std = np.std(grid)
    grid = grid / std * np.sqrt(sill) + mu
    grid = rotate(grid, angle, order=3, cval=0)
    lowerleft = int(grid.shape[0] / 2) - ox / 2 , int(grid.shape[1] / 2 - oy / 2)
    grid = grid[lowerleft[0]: lowerleft[0] + ox, lowerleft[1]: lowerleft[1] + oy]
    return grid

#def truncate_grf(grid,fac_proportions):
