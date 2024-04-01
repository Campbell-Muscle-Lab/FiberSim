---
layout: default
title: KX=F
grand_parent: FiberCpp
parent: Half-sarcomeres
has_children: false
nav_order: 3
---

# KX=F

## Overview

[Lattice equations](../lattice_equations/lattice_equations.html) showed how force-balance can be applied to determine the position of each node. This process yields a system of equations that can be written as

$$ K X = F $$

where:
 + $ K $ is a square matrix of size $ n $ x $ n $
 + $X$ is a vector of size $ n $ x 1
 + $F$ is a vector of size $ n $ x 1

 $K$ and $F$ can be calculated as described at [lattice equations](../lattice_equations/lattice_equations.html). The goal is to calculate $X$. This process can be written as:
 + $X = K \backslash F$

## Size of the problem

The main challenge is the size of the problem. Consider the smallest possible lattice which has just 4 thick filaments.

+  $n$ = (no of thick filaments * no of nodes per thick filament) + (no of thin filaments * no of nodes per thin filament)
+ $n$ = (4 * 54) + (8 * 189) = 1728

Doing anything 1728 times will require some resources. In practice, simulations often need 100 or more thick filaments to produce smooth records at low levels of activation. With 100 filaments, $n$ has increased to 43200.

Standard approaches to solve $X = K \backslash F$ scale as $n^3$ implying that calculating the node positions for 100 thick filaments will require 73 million operations. And that is for each time-step! The computational cost becomes prohibitive.

## Strategy

As explained at [lattice equations](../lattice_equations/lattice_equations.html), the $K$ matrix is sparse, meaning that most of the entries are 0. FiberSim takes advantage of this and uses special techniques that are designed for such matrices. In fact, the main code-base never stores the full $K$ matrix in so-called dense form.

### Principles

+ If a sparse matrix $A$ is stored in [Compressed Sparse Column (CSC)](https://www.gnu.org/software/gsl/doc/html/spmatrix.html#compressed-sparse-column-csc) form, and $B$ is a vector, the operation $A.B$ can be performed efficiently.

+ In the special case where a matrix $T$ is tri-diagonal (values only on the leading diagonal, and diagonals immediately above and immediately below), there is an efficient algorithm to solve $X = T \backslash F$

+ As explained at [lattice equations](../lattice_equations/lattice_equations.html), the $K$ matrix for FiberSim is dominated by the tridiagonals.

### Approach

The goal is to solve $X$ where $K$ and $F$ are known and 
$$ K X = F $$

Define $K_0$ as the portion of $K$ that correspond to the thin and thick filament backbones, noting that $K_0$'s non-zero elements all lie on the tridiagonals.

Define $G(X)$ as a vector that depends on X and represents perturbations induced by cross-links (in FiberSim's case, titin, bound myosins, and bound cMyBP-C molecules).

$$ K_0 X + G(X) = F $$

Further define

$$ X_i = X_{g,i} + \Delta X $$

where

+ $X_i$ is the 'true' $X$
+ $X_{g,i}$ is the i'th 'guess' at $X_i$
+ $\Delta X_i$ is the difference between them.

Substituting in, leads to

$ K_0 (X_{g,i} + \Delta X_i) + G(X_{g,i} + \Delta X_i) = F $

This can be re-arranged to yield

$ \Delta X  = K_0 ^{-1} (F - K_0 X_{g,i} - G(X_{g,i})) $

An improved guess for $X$ can then be calculated from

$ X_{g,i+1} = X_{g,i} + \Delta X_i $

Note that

+ $ K_0 X_{g,i} $ can be evaluated efficently if $K_0$ is stored in compressed sparse column format

+ $ (F - K_0 X_{g,i} - G(X_{g,i})) $ simplifies to $ some \ vector $

+ $ \Delta X  = K_0 ^{-1} (some\ vector) $ can be calculated efficiently because $K_0$ is tridiagonal

This leads to an iterative approach that runs as follows:

+ Start with an initial guess for $X_{g,i}$
  + FiberSim uses the tridiagonal solution for the first time-step and $X$ from the last time-step thereafter.
+ In a loop
  + calculate $G(X_{g,i})$
  + calculate $\Delta X_i$
  + form a new $X_{g,i+1}$ = $X_{g,i} + \alpha \Delta X_i$ 
  + repeat until $X_g$ converges

Since FiberSim Ver 2.3, $ \alpha$ is initialized at 1 but can be reduced to improve convergence. Specifically, if $X$ starts to diverge, $ \alpha $ is multiplied by 0.5 for the next step. This adds a smaller correction.

This approach

+ maintains fast convergence for situations where the filaments are relatively stiff
+  improves the likelihood of finding the correct $X$ vector when the filaments are so compliant that the perturbations induced by cross-links are comparable to the inter-node spacing.


#### An aside that didn't work

FiberCpp uses the [GSL Scientific Library](https://www.gnu.org/software/gsl/doc/html/index.html) for most of the scientific computing procedues.

GSL includes code that can solve $X = K \backslash X$ when $K$ is stored in a sparse format.

GSL's [approach](https://www.gnu.org/software/gsl/doc/html/splinalg.html) uses a Generalized Minimum Residual Method (GMRES) algorithm. Early tests showed that the approach correctly calcuated $X$ but that the algorithm was so slow that the calculations became impractical.

Alternative methods may be more efficent.

