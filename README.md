# Finite Element Simulation in Porous Media

A MATLAB-based finite element implementation for simulating acoustic wave propagation in porous and heterogeneous media.

This repository contains several educational FEM examples developed for EN234 (Brown University, 2017). The codes demonstrate how to construct a finite element solver for acoustic Helmholtz problems, including mesh generation, stiffness matrix assembly, porous acoustic modeling, and multi-material domain simulations.

---

# Repository Structure

```text
FEM_Final.m
```

Basic 2D acoustic FEM simulation in a homogeneous medium with comparison against the analytical 1D solution.

```text
FEM_Final_example2.m
```

Acoustic propagation around a circular void using an unstructured triangular mesh.

```text
FEM_Final_example3.m
```

Thermoviscous porous acoustic model with complex-valued effective material properties.

```text
FEM_Final_example4_low.m
```

Multi-material acoustic simulation with spatially varying porous media properties.

---

# Governing Equation

The simulations solve the frequency-domain acoustic Helmholtz equation:

$$
\nabla^2 p + k^2 p = 0
$$

where:

- $p$ is the acoustic pressure
- $k$ is the acoustic wavenumber

The acoustic wavenumber is defined as:

$$
k = \frac{2\pi f}{c_0}
$$

where:

- $f$ is the frequency
- $c_0$ is the sound speed in air

For the 1D validation case, the analytical solution is:

$$
p(x) = \frac{\cos(kx)}{\cos(kL)}
$$

---

# Numerical Method

The implementation uses:

- Linear triangular finite elements
- Delaunay triangulation
- Gaussian integration
- Direct stiffness matrix assembly
- Sparse matrix visualization

The elemental FEM formulation is:

$$
K_e =
\int_{\Omega}
\left(
\nabla N \cdot \nabla N
-
k^2 N N
\right)
d\Omega
$$

where:

- $N$ denotes the shape functions
- $k$ is the local acoustic coefficient

---

# Example 1 — Homogeneous Acoustic Medium

This example demonstrates the basic FEM solution of acoustic wave propagation in a rectangular domain.

## Features

- Structured rectangular mesh
- Triangular element generation
- Global stiffness assembly
- Dirichlet boundary condition
- Comparison with analytical solution

## Output

The code generates:

- FEM mesh
- Sparsity pattern of stiffness matrix
- Reordered sparse matrix using `symrcm`
- 2D pressure field
- Numerical vs analytical comparison

---

# Example 2 — Acoustic Scattering Around a Circular Void

This example introduces an internal circular cavity inside the computational domain.

Elements located inside the circular region are removed after Delaunay triangulation.

## Features

- Unstructured mesh generation
- Internal obstacle removal
- Acoustic scattering behavior
- Pressure field masking using `NaN`

## Geometry

The circular void satisfies:

$$
(x-L/2)^2 + (y-H/2)^2 < R^2
$$

---

# Example 3 — Thermoviscous Porous Acoustic Medium

This example incorporates thermoviscous effects in narrow slit-like porous channels.

The porous material is modeled using complex-valued effective density and effective bulk modulus.

## Effective Acoustic Model

The effective acoustic coefficient is:

$$
k =
\omega
\sqrt{
\frac{\rho_{\mathrm{eff}}}
{K_{\mathrm{eff}}}
}
$$

where:

- $\rho_{\mathrm{eff}}$ is the effective density
- $K_{\mathrm{eff}}$ is the effective bulk modulus
- $\omega = 2\pi f$

## Features

- Complex-valued FEM system
- Frequency-dependent attenuation
- Real and imaginary pressure comparison
- Thermoviscous porous media modeling

---

# Example 4 — Multi-Material Acoustic Domain

This example demonstrates wave propagation through multiple acoustic regions with different porous properties.

Different element coefficients are assigned according to element centroid positions.

## Features

- Spatially varying acoustic properties
- Multiple porous media regions
- Heterogeneous FEM assembly
- Complex-valued acoustic propagation

---

# Material Parameters

Typical parameters used in the simulations:

```matlab
density = 1.23;              % air density [kg/m^3]
viscosity = 1.95e-5;         % dynamic viscosity [Pa·s]
c0 = 343;                    % speed of sound [m/s]
p0 = 1.013e5;                % ambient pressure [Pa]
gamma = 1.4;                 % heat capacity ratio
```

---

# Running the Code

Run any example directly in MATLAB:

```matlab
FEM_Final
```

or

```matlab
FEM_Final_example2
```

etc.

---

# Educational Purpose

These examples were primarily developed as educational finite element projects exploring:

- acoustic wave propagation
- porous acoustic materials
- thermoviscous losses
- heterogeneous media
- unstructured triangular FEM implementation

The implementation prioritizes readability and transparency over computational efficiency.

---

# Possible Future Improvements

Potential extensions include:

- absorbing boundary conditions
- impedance boundary conditions
- sparse matrix optimization
- higher-order elements
- frequency sweep automation
- PML (Perfectly Matched Layer)
- 3D acoustic FEM
- coupled structural-acoustic simulations

---

# Author

Siyuan Song  
Brown University  
EN234 Final Project  
December 2017
