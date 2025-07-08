# 2D Diffusion Equation Solver

This project implements a numerical solver for the **two-dimensional diffusion (heat) equation** using a **finite difference method** (FDM) and explicit time-stepping.

## Mathematical Model

We solve the unsteady diffusion equation:

\[
\frac{\partial u}{\partial t} = \alpha \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right), \quad (x, y) \in \Omega \subset \mathbb{R}^2,\ t > 0
\]

with:
- \( u(x, y, t) \): temperature or scalar concentration field
- \( \alpha \): diffusion coefficient
- Appropriate **Dirichlet** or **Neumann** boundary conditions
- Initial condition: \( u(x, y, 0) = u_0(x, y) \)

## Numerical Method

- **Spatial discretization**:
  - Second-order centered finite differences for the Laplacian
  - Uniform Cartesian mesh
- **Time discretization**:
  - Explicit Forward Euler method
  - Time step limited by CFL condition for stability

The update scheme at each grid point \((i, j)\) is:

\[
u_{i,j}^{n+1} = u_{i,j}^n + \Delta t \cdot \alpha \left(
\frac{u_{i+1,j}^n - 2u_{i,j}^n + u_{i-1,j}^n}{\Delta x^2} +
\frac{u_{i,j+1}^n - 2u_{i,j}^n + u_{i,j-1}^n}{\Delta y^2}
\right)
\]

## Features

- Handles both **Dirichlet** and **Neumann** boundary conditions
- Written in modular Python
- Visualization of the diffusion process (animated and static)
- Sparse matrix construction of the Laplacian operator
- Optional use of **LU decomposition** via `scipy.sparse.linalg` for validation or implicit extensions

## Files

- `main.py`: solver and visualization
- `utils.py`: setup routines and boundary condition handling
- `laplacian.py`: constructs the 2D Laplacian operator using sparse matrices

## How to Run

```bash
python main.py
