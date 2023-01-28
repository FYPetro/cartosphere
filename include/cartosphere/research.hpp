
#ifndef __RESEARCH_HPP__
#define __RESEARCH_HPP__

#include "cartosphere/solver.hpp"

/* Build a linear system for the demo program */
void build_system(const Cartosphere::TriangularMesh& mesh, Matrix& A, Vector& b);

/* Demo */
int demo();

/* Demo Diffusion */
int demo_diffusion();

/* Quadrature rule demo */
int demo_quadrature();

/* Seminar code */
int seminar();

/* Verification of Convergence */
int convergence();

/* Generate weights */
int precompute_weights(const string& path);

/* Testing iterative refinement */
int refine(const string& path);

/* Test the coloring of obj */
int test_obj();

/* Validation, for research only */
int research_a();

// Sep 16, 2021. To gauge the error of FEM correctly.
// The goal is to start my preliminaries!!! HOPE FOR THE BEST!!!
int research_b();

// Sep 23, 2021. To test the steady state problem for the eigenfunctions
// Again, hope for the best!
int research_c(int l = 0, int m = 0, bool silent = false);

// Jan 03, 2022. To test the output of individual spherical triangles.
int research_d();

// Jan 04, 2022. Incorporate second-order time discretization.
int research_f();

// Jan 24, 2023. Read shapefiles
int research_g(const string &folder);

// Jun 06, 2022. To replicate the original paper
int benchmark();

#endif // !__RESEARCH_HPP__
