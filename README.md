# 2D Acoustic Field FEM Demo

This repository contains a MATLAB finite element implementation for solving a simple 2D acoustic pressure field problem.

The code was originally written as a final project for EN234 in December 2017. It demonstrates the basic workflow of a finite element solver, including mesh generation, element stiffness assembly, boundary condition enforcement, solution of the global linear system, and comparison with an analytical 1D acoustic solution.

## Problem Description

The program solves a frequency-domain acoustic wave problem in a rectangular domain.

The governing equation is a Helmholtz-type equation for acoustic pressure:

$$
\nabla^2 p + k^2 p = 0
$$

where $p$ is the acoustic pressure and

$$
k = \frac{2\pi f}{c_0}
$$

is the acoustic wavenumber.

A prescribed background pressure is applied at the left boundary. The numerical solution is compared with the analytical 1D solution:

$$
\nabla^2 p + k^2 p = 0
$$

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
# Additional Examples

This repository also contains several extended examples demonstrating different acoustic FEM configurations and porous media effects.

---

# Example 2 — Acoustic Scattering Around a Circular Void

This example investigates acoustic wave propagation in a rectangular domain containing a circular obstacle (void region).

## Description

A structured rectangular grid is first generated and triangulated using Delaunay triangulation. Elements whose centroids fall inside the circular region are removed to create the internal cavity.

The simulation demonstrates:

- unstructured mesh generation
- element filtering using geometric criteria
- acoustic scattering around an internal obstacle
- pressure field visualization with irregular geometry

## Geometry

- Rectangular domain
- Circular void at the center

The circular region satisfies:

$$
(x-L/2)^2 + (y-H/2)^2 < R^2
$$

where:

- $L$ = domain length
- $H$ = domain height
- $R$ = circular void radius

## Features

- Automatic mesh generation
- Delaunay triangulation
- Internal obstacle removal
- 2D pressure field visualization
- NaN masking for void regions

## Output

The code produces:

1. Finite element mesh with internal cavity
2. Sparsity pattern of the global stiffness matrix
3. Reordered sparse matrix using `symrcm`
4. 2D pressure distribution around the obstacle

---

# Example 3 — Thermoviscous Porous Acoustic Medium

This example extends the acoustic FEM formulation to include thermoviscous effects inside narrow porous channels.

## Description

Instead of using the standard acoustic wavenumber, the model computes an effective complex density and effective bulk modulus using a slit-pore thermoviscous model.

The resulting acoustic coefficient becomes complex-valued:

$$
k = \omega \sqrt{\frac{\rho_{\mathrm{eff}}}{K_{\mathrm{eff}}}}
$$

where:

- $\rho_{\mathrm{eff}}$ is the effective density
- $K_{\mathrm{eff}}$ is the effective bulk modulus
- $\omega = 2\pi f$

## Thermoviscous Effects

The model accounts for:

- viscous boundary layer effects
- thermal boundary layer effects
- frequency-dependent attenuation
- phase lag in porous acoustic propagation

The slit radius controls the strength of these effects.

## Features

- Complex-valued FEM system
- Frequency-dependent porous material properties
- Real and imaginary pressure field comparison
- Validation against analytical 1D solution

## Output

The code produces:

1. FEM mesh
2. Sparse stiffness matrix visualization
3. Real part of pressure field
4. Imaginary part of pressure field
5. Comparison with analytical solution

---

# Example 4 — Multi-Material Acoustic Domain

This example demonstrates acoustic wave propagation through multiple regions with different acoustic properties.

## Description

The computational domain is divided into several subregions:

- free air region
- porous material region 1
- porous material region 2

Different effective acoustic coefficients are assigned to different regions based on element centroid locations.

## Material Assignment

The code assigns material properties according to spatial position:

- Left region:
  - standard air acoustic coefficient
- Right upper/lower regions:
  - thermoviscous porous material
- Right middle region:
  - different porous material with larger slit radius

This creates a heterogeneous acoustic medium.

## Features

- Spatially varying material properties
- Multiple porous acoustic domains
- Complex-valued acoustic propagation
- Region-dependent element stiffness assembly

## Output

The code produces:

1. Multi-material FEM mesh
2. Sparse stiffness matrix
3. Reordered matrix visualization
4. 2D acoustic pressure distribution

---

# Numerical Method

All examples use:

- Linear triangular finite elements
- Gaussian integration
- Helmholtz-type acoustic formulation
- Direct linear system solution

The element matrix is assembled from:

$$
K_e = \int_{\Omega}
\left(
\nabla N \cdot \nabla N
-
k^2 N N
\right)
d\Omega
$$

where:

- $N$ are the shape functions
- $k$ is the acoustic wavenumber

---

# Educational Purpose

These examples were originally developed as educational FEM projects exploring:

- acoustic wave propagation
- porous acoustic materials
- thermoviscous effects
- heterogeneous media
- unstructured mesh FEM implementation in MATLAB

The implementation prioritizes clarity and accessibility over computational efficiency.
