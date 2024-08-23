---
sidebar_position: 1
---

# Introduction

Pytop is a powerful tool designed to extend the capabilities of FEniCS for topology optimization. It provides robust tools for advanced computational design and supports multiphysics simulations. Pytop is fully compatible with MPI, allowing for efficient monitoring of parallel applications.

## General Optimization Form

This software scopes following general optimization form:
$$
\min_{\bm d} F(\bm d, \bm u) = \int_\Omega f(\bm d, \bm u) \ dV
$$
$$
G_j(\bm d, \bm u) \leq 0
$$
where $\Omega$ is the design domain and $\bm u$ represents specific fields that satisfy governing equations for any physics defined in $\Omega$. For example, $\bm u$ is the displacement vector in an elastic deformation problem. $F$ is the functional that represents the objective of the optimization. $\bm d$ is the design vector, and $G_j$ are the functionals representing the constraints.

## NLOpt ```opt``` class

This software is designed to construct the ```nlopt.opt``` class in [NLOpt](https://nlopt.readthedocs.io/en/latest/NLopt_Python_Reference/) Python interface.
You define design variables as ```DesignVariables``` class and optimization problem with ```PloblemStatement``` class. Then, extended ```nlopt.opt``` class will be defined using ```NloptOptimizer``` class with
```DesignVariables``` and ```PloblemStatement```. Some local gradient-based optimizers in Nlopt are available, including **GCMMA**, **SLSQP**, and **CCSA**. Please refer detail of each algorithms in [NLopt algorithms](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/). 

## Visualization
Currently, this software does not include unique vilualization methods. However, The XDMF and H5 files can be exported. XDMF file can be visualized the VTK-based visualization software, like [paraview](https://www.paraview.org/).


## Container

We provide a container for this repository. The container includes python 3.11, FEniCS bundles, and NLOpt with python interface.
The container is avaiable in [dockerhub](https://hub.docker.com/repository/docker/ichiharanaruki/pytop/general).

## Installation

To install Pytop, you need to have Python and FEniCS installed on your system. You can install Pytop via pip:

```bash
pip install git+https://github.com/Naruki-Ichihara/pytop.git@main
```

More detailed installation can be found in [this page](/docs/Getting-started/Installing-pytop).

## MPI parallelization
This software is designed as parallelization-ready. Problem is automatically divided into partial problems and computed in each cpus. We provide useful detaclass for MPI, `MPI_Communicator`. You can construct a mpi group use this class,
```python
import pytop as pt
comm_group = pt.MPI_Communicator.comm_world
```
Single problem must be in the same group. Then, `mpirun` executes the solving with cpus you mentioned:
```bash
mpirun -n 36 python your_problem.py
```
In the container, you need to give option to allow run as a root user:
```bash
mpirun -n 36 --allow-run-as-root python your_problem.py
```
Above example use 36 cpus. 

## Working with pygmsh
You can use ```pygmsh``` as mesh generator.
Documentation of ```pygmsh``` is available [here](https://pygmsh.readthedocs.io/en/latest/).
```python
import pytop as pt
import pygmsh

with pygmsh.geo.Geometry() as geom:
    lcar = 0.1
    p1 = geom.add_point([0.0, 0.0], lcar)
    p2 = geom.add_point([1.0, 0.0], lcar)
    p3 = geom.add_point([1.0, 0.5], lcar)
    p4 = geom.add_point([1.0, 1.0], lcar)
    s1 = geom.add_bspline([p1, p2, p3, p4])

    p2 = geom.add_point([0.0, 1.0], lcar)
    p3 = geom.add_point([0.5, 1.0], lcar)
    s2 = geom.add_spline([p4, p3, p2, p1])

    ll = geom.add_curve_loop([s1, s2])
    pl = geom.add_plane_surface(ll)

    mesh = geom.generate_mesh()

mesh = pt.from_pygmsh(mesh)
```

Use ```planation``` option is used, all ```z``` coordinates will be neglected.
```python
mesh = pt.from_pygmsh(mesh, planation=True)
```