# 2D Acoustic Field FEM Demo

This repository contains a MATLAB finite element implementation for solving a simple 2D acoustic pressure field problem.

The code was originally written as a final project for EN234 in December 2017. It demonstrates the basic workflow of a finite element solver, including mesh generation, element stiffness assembly, boundary condition enforcement, solution of the global linear system, and comparison with an analytical 1D acoustic solution.

## Problem Description

The program solves a frequency-domain acoustic wave problem in a rectangular domain.

The governing equation is a Helmholtz-type equation for acoustic pressure:

\[
\nabla^2 p + k^2 p = 0
\]

where \(p\) is the acoustic pressure and \(k = 2\pi f / c_0\) is the acoustic wavenumber.

A prescribed background pressure is applied at the left boundary. The numerical solution is compared with the analytical 1D solution:

\[
p(x) = \frac{\cos(kx)}{\cos(kL)}
\]

## Features

- Generates a 2D rectangular grid
- Uses Delaunay triangulation to create triangular elements
- Assembles the global FEM stiffness matrix
- Applies Dirichlet boundary conditions
- Solves the acoustic pressure field
- Visualizes:
  - triangular mesh
  - sparsity pattern of the stiffness matrix
  - pressure field
  - numerical vs. analytical pressure profile

## Main File

```text
FEM_Final.m
