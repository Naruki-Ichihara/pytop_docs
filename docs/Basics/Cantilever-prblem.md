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

### Implimentation
First, we import pytop package as:
```python
import pytop as pt
```

:::warning
Do not import FEniCS like `from fenics import *` because all method in FEniCS are overloaded by *pyadjoint* and imported in `__init__.py` of pytop.
:::

Second 
