# Adaptive Finite Element Method (FEM) Project

This repository contains the implementation of an Adaptive Finite Element Method (FEM) solver for solving differential equations, along with a complete writeup detailing all mathematical derivations, methodologies, and results.

## Overview

The project employs an adaptive finite element method with a posteriori error estimation for solving boundary value problems of the form:

-u''(x) + m * u'(x) + n * u(x) = f(x), for all x in (a, b)

The associated dissertation provides a full description of the mathematics behind the adaptive method, the a posteriori error estimation process, and the practical implementation details of the algorithm.

### Summary from Dissertation

The aim of this project was to code and analyse the efficacy of an adaptive hp finite element method in solving a 1-dimensional linear differential equation with Dirichlet boundary conditions. In order to do this, a general hp finite element method with a suitable piecewise polynomial approximation was implemented and demonstrated to be a valid model. An adaptive algorithm was then designed such that, using an error estimate based on the residual of the differential equation posed, the model could automatically be refined to suit its local behaviours across the domain. This complete model was then analysed, and it was determined that, although a global refinement method was more efficient for low accuracy approximations or particularly smooth functions, the adaptive finite element method was extremely effective and efficient for high accuracy and badly behaved function solutions.

## Project Structure

- `code/`: Contains the Octave scripts for running the adaptive FEM.
  - `AdaptiveFiniteElementMethod_A_Posteriori.m`: The main script for setting up and running the FEM solver with adaptive refinement.
  - `IterateFiniteElementMethod.m`: Function used for iteratively performing finite element analysis.
- `report/`: Contains the dissertation writeup.
  - `final_report.pdf`: Full writeup that includes the complete description of the mathematics and implementation details.

## Running the Code

The provided scripts are designed to be run in **GNU Octave**, as the author no longer has access to MATLAB. The code was originally written in MATLAB but has been adapted to work in Octave, which is a free alternative. The code was originally developed in MATLAB but is fully compatible with Octave.

To run the main script in Octave:

```sh
cd code
octave AdaptiveFiniteElementMethod_A_Posteriori.m
```

## Requirements

- [GNU Octave](https://www.gnu.org/software/octave/)


