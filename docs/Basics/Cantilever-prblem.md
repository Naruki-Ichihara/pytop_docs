---
sidebar_position: 1
---

# Cantilever problem
In this example, a linear (2D plane stress) elasticity problem is considered. First, we show forward definition of the problem using FEniCS in pytop. Then, we describe the
how to build the optimization problem using pytop's features.

## Forward Problem
### Governing Equation
The equation for linear elastic deformation problem can be written as

$$
-\nabla \cdot \sigma = {\bm f} \ {\rm in}\  \Omega
$$

$$
\sigma = \lambda {\rm tr} (\varepsilon)I+2\mu \varepsilon
$$

$$
\varepsilon = \frac{1}{2}\left(\nabla {\bm u}+(\nabla {\bm u})^{\rm T} \right)
$$

where $\sigma$ is the stress tensor, ${\bm f}$ is a body force, and ${\bm u}$ is the deformation vector. The strain $\varepsilon$ can be calculated from the ${\bm u}$. 
The material property is described as Lame's constants, $\lambda$ and $\mu$. All variables are defined in closed body $\Omega$.

### Variational Form
The variational (or weak) form of abovementioned equation is obtained by the inner product of the equation and *test function* ${\bm v}$. We get

$$
\int_\Omega\sigma : \nabla {\bm v}\  dV = \int_{\partial \Omega}{\bm f}\cdot{\bm v}\ dV + \int_{\partial \Omega}\left( \sigma \cdot {\bm n}\right)\cdot {\bm v}\ dS.
$$

Here, $\left( \sigma \cdot {\bm n}\right)={\bm t}$ is a traction force on the Dirichlet bondary ${\partial \Omega_D}$. This form can summarize as: find ${\bm u}\in V$ that satified

$$
a({\bm u},{\bm v})=L({\bm v})
$$

$$
a({\bm u},{\bm v})=\int_\Omega\sigma : \nabla {\bm v}\  dV
$$

$$
L({\bm v})=\int_{\partial \Omega}{\bm f}\cdot{\bm v}\ dV + \int_{\partial \Omega}\left( \sigma \cdot {\bm n}\right)\cdot {\bm v}\ dS.
$$

Because of symmetrity of the $\sigma$, the bilinear form $a({\bm u},{\bm v})$ can be replaced into elastic potential:

$$
a({\bm u},{\bm v})=\int_\Omega\sigma({\bm u}) : \varepsilon({\bm v})\  dV
$$

## Optimization problem
### Design variable
The density variable $\rho({\bm x})$ is defined to represent the material distribution. We parametarize the modulus $E_\rho$ with $\rho({\bm x})$ as following:

$$
E_\rho=((1-\epsilon)\rho^p + \epsilon)E_s={\rm simp}(\rho; p)E_s,
$$

where $E_s$ and $p$ are bulk modulus and penalty parameter, and $\epsilon$ is small value to avoid the computational problem. In the isotropic elasticity problem, the bilinear form $a({\bm u},{\bm v})$
can be rewritten as:

$$
a({\bm u},{\bm v},\rho)=\int_\Omega\sigma({\bm u}){\rm simp}(\rho; p) : \varepsilon({\bm v})\  dV.
$$


### Objective
Here, compliance energy minimization problem is considered. The objective fuction $F$ is

$$
F = a({\bm u_s},{\bm u_s},\rho)=\int_\Omega\sigma({\bm u_s}){\rm simp}(\rho; p) : \varepsilon({\bm u_s})\  dV
$$

where ${\bm u_s}$ is a solution deformation of elastic problem.

## Implimentation
### import pytop
First, you import pytop and pygmsh packages as:
```python
import pytop as pt
import pygmsh
```

:::warning
Do not import FEniCS like `from fenics import *` because all method in FEniCS are overloaded by *pyadjoint* and imported in `__init__.py` of pytop.
:::

### Parameters

```python
TARGET_DENSITY = 0.25
FILTER_RADIUS = 1.0
RECTANGLE_SHAPE = (20, 10)
MESH_SIZE = 0.2
FORCE = -1e3
E = 1e6  # Youngs modulus
nu = 0.3 # Poisson's ratio
```

### Mesh

```python
with pygmsh.geo.Geometry() as geom:
    # generate polygon
    L = RECTANGLE_SHAPE(0)
    H = RECTANGLE_SHAPE(1)
    geom.add_polygon([
        [0, 0],
        [L, 0],
        [L, H],
        [0, H]
    ], mesh_size=MESH_SIZE)
    gmesh = geom.generate_mesh()

mesh = pt.from_pygmsh(mesh, planation=True)
```
`from_pygmsh` method convert `pygmsh` to `dolfin-xml` mesh format. `planation` option neglects all z coordinates and convert to 2D mesh.

### Function spaces and functions
Two function spaces, for displacement vector and density field, are defined.

```python
U = pt.VectorFunctionSpace(mesh, "CG", 1) # Displacement
X = pt.FunctionSpace(mesh, "CG", 1)       # Density
```
Both finction spaces are continuous Garelkin space (1st order Lagrange). The functions are defined as:
```python
us = pt.Function(U, name="displacement")
u = pt.TrialFunction(U)
du = pt.TestFunction(U)
```
### Boundary condition

```python
class Left(pt.SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < 1e-6 and on_boundary
bc = pt.DirichletBC(U, pt.Constant((0, 0)), Left())

class Loading(pt.SubDomain):
    def inside(self, x, on_boundary):
        CENTER = RECTANGLE_SHAPE(1)/2
        return x[0] > RECTANGLE_SHAPE(1)-1e-6 and CENTER-0.5 < x[1] < CENTER+0.5 and on_boundary

bc = pt.DirichletBC(U, pt.Constant((0, 0)), Left())
ds = pt.make_noiman_boundary_domains(mesh, [Loading()], True)
```
`ds` is `Measure` in FEniCS, please refer [FEniCS docs](https://fenicsproject.org/olddocs/dolfin/latest/python/demos/subdomains-poisson/documentation.html?highlight=measure). `make_noiman_boundary_domains` method sets all boundary as 0, then, set 1, 2, ... with given subdomain's list. In this example, the subdomain `Loading` is 1, and other bondaries are 0.

### Design Variables

```python
design_variables = pt.DesignVariables()
design_variables.register(X,
                          "density",
                          [TARGET_DENSITY],
                          [(0, 1)],
                          lambda x: pt.helmholtz_filter(x, R=FILTER_RADIUS),
                          recording_path="PATH FOR RECORDING",
                          recording_interval=1)
```
The desing variable name is "density". The initial value is `TARGET_DENSITY` and range is from 0 to 1.
The helmholtz filter is used for regularization. Also you can provide directory path to save history.

### Problem statement

```python
from pytop.physics import elasticity as el
from pytop.physics.utils import penalized_weight

class Problem(pt.ProblemStatement):
    def objective(self, design_variables):
        self.rho = design_variables["density"]
        a = el.linear_2D_elasticity_bilinear_form(u, du, E, nu, penalized_weight(self.rho, eps=1e-4))
        L = pt.inner(f, du) * ds(1)
        pt.solve(a == L, uh, bc)
        return pt.assemble(pt.inner(f, uh) * ds(1))
    def constraint_volume(self, design_variables):
        unitary = pt.project(pt.Constant(1), X)
        return pt.assemble(self.rho*pt.dx)/pt.assemble(unitary*pt.dx) - TARGET_DENSITY
```
`python.physics` has some bilinear form of built-in physics. Here, `linear_2D_elasticity_bilinear_form` is used to define 2D elasticity problem.

### Optimizer
```python
opt = pt.NloptOptimizer(design_variables, Problem(), "LD_MMA")
opt.set_maxeval(300)
opt.set_ftol_rel(1e-4)
opt.set_param("verbosity", 1)
opt.run("PATH FOR GROWTH CURVE")
```
`NloptOptimizer` class constructs `nlopt.opt` class with design variables and statements. Some local gradient-based optimizers in Nlopt are available, including **GCMMA**, **SLSQP**, and **CCSA**. Please refer detail of each algorithms in [NLopt algorithms](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/). Here, we select **GCMMA** with keyword "LD_MMA".

## Complete code
```python
import pytop as pt
import pygmsh

TARGET_DENSITY = 0.25
FILTER_RADIUS = 1.0
RECTANGLE_SHAPE = (20, 10)
MESH_SIZE = 0.2
FORCE = -1e3
E = 1e6  # Youngs modulus
nu = 0.3 # Poisson's ratio

with pygmsh.geo.Geometry() as geom:
    # generate polygon
    L = RECTANGLE_SHAPE(0)
    H = RECTANGLE_SHAPE(1)
    geom.add_polygon([
        [0, 0],
        [L, 0],
        [L, H],
        [0, H]
    ], mesh_size=MESH_SIZE)
    gmesh = geom.generate_mesh()

mesh = pt.from_pygmsh(mesh, planation=True)

U = pt.VectorFunctionSpace(mesh, "CG", 1) # Displacement
X = pt.FunctionSpace(mesh, "CG", 1)       # Density

us = pt.Function(U, name="displacement")
u = pt.TrialFunction(U)
du = pt.TestFunction(U)

class Left(pt.SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < 1e-6 and on_boundary
bc = pt.DirichletBC(U, pt.Constant((0, 0)), Left())

class Loading(pt.SubDomain):
    def inside(self, x, on_boundary):
        CENTER = RECTANGLE_SHAPE(1)/2
        return x[0] > RECTANGLE_SHAPE(1)-1e-6 and CENTER-0.5 < x[1] < CENTER+0.5 and on_boundary

bc = pt.DirichletBC(U, pt.Constant((0, 0)), Left())
ds = pt.make_noiman_boundary_domains(mesh, [Loading()], True)

design_variables = pt.DesignVariables()
design_variables.register(X,
                          "density",
                          [TARGET_DENSITY],
                          [(0, 1)],
                          lambda x: pt.helmholtz_filter(x, R=FILTER_RADIUS),
                          recording_path="PATH FOR RECORDING",
                          recording_interval=1)

from pytop.physics import elasticity as el
from pytop.physics.utils import penalized_weight

class Problem(pt.ProblemStatement):
    def objective(self, design_variables):
        self.rho = design_variables["density"]
        a = el.linear_2D_elasticity_bilinear_form(u, du, E, nu, penalized_weight(self.rho, eps=1e-4))
        L = pt.inner(f, du) * ds(1)
        pt.solve(a == L, uh, bc)
        return pt.assemble(pt.inner(f, uh) * ds(1))
    def constraint_volume(self, design_variables):
        unitary = pt.project(pt.Constant(1), X)
        return pt.assemble(self.rho*pt.dx)/pt.assemble(unitary*pt.dx) - TARGET_DENSITY

opt = pt.NloptOptimizer(design_variables, Problem(), "LD_MMA")
opt.set_maxeval(300)
opt.set_ftol_rel(1e-4)
opt.set_param("verbosity", 1)
opt.run("PATH FOR GROWTH CURVE")
```