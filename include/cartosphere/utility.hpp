
#ifndef __UTILITY_HPP__
#define __UTILITY_HPP__

// Detect compiler
// Ref:	https://stackoverflow.com/a/5920028/1377770
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
	#define IS_WINDOWS
	#if defined(_DEBUG)
		#define BUILD_DEBUG
	#else
		#define BUILD_RELEASE
	#endif
#elif defined(__APPLE__) || defined(__linux__) || defined(__unix__) || defined(_POSIX_VERSION)
	#define APPLE_LIKE
#else
#   error "Unknown compiler"
#endif

#ifdef IS_WINDOWS
#define NOMINMAX
#endif

// Use std::forward
#include <utility>

// Lifecycle management for dynamically allocated objects
#include <memory>
using std::shared_ptr;
using std::make_shared;

// Obtain a few numeric limits
#include <limits>
static const double DoubleEpsilon = std::numeric_limits<double>::epsilon();
static const double DoubleMinimum = std::numeric_limits<double>::min();
static const double DoubleMaximum = std::numeric_limits<double>::max();

// OpenMP Number of threads
#ifdef _OPENMP
#include <omp.h>
static const int ThreadsMaximum = omp_get_max_threads();
#endif

// Commonly used templates
#include <vector>
using std::vector;

#include <string>
using std::string;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <sstream>
using std::stringstream;
using std::istringstream;

// Terminal output
#include <iostream>
#include <iomanip>

// Filesystem
#include <filesystem>
using std::filesystem::path;

// Timing
#include <chrono>
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

#pragma warning(push)
#pragma warning(disable: 4819)
#include <Eigen/Dense>
// Import Eigen column arrays
using ColArray = Eigen::ArrayXd;
// Import Eigen row arrays
using RowArray = Eigen::Array<double, 1, Eigen::Dynamic>;
// Import Eigen column vector
using ColVector = Eigen::VectorXd;
// Import Eigen row vector
using RowVector = Eigen::RowVectorXd;
// Import Eigen dense matrix (column-major)
using Matrix = Eigen::MatrixXd;
// Import Eigen dense matrix (row-major)
using MatrixRowMajor = Eigen::Matrix<
	double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor
>;

#include <Eigen/SparseCore>
// Import Eigen sparse matrix (row-major)
using SparseMatrixRowMajor = Eigen::SparseMatrix<double, Eigen::RowMajor>;
// Provide convenient typedef for coordinates
using SparseMatrixEntry = Eigen::Triplet<double>;

#include <Eigen/IterativeLinearSolvers>
using SolverBiCGSTAB = Eigen::BiCGSTAB<
	SparseMatrixRowMajor, Eigen::IncompleteLUT<SparseMatrixRowMajor::Scalar>
>;
#pragma warning(pop)

// Numerics, algorithms, and functionals
#include <numeric>
#include <algorithm>
#include <functional>

// Use glog to debug things
// Windows output: %APPDATA%/../Local/Temp
#define GLOG_NO_ABBREVIATED_SEVERITIES
#include <glog/logging.h>

#endif // !__UTILITY_HPP__
