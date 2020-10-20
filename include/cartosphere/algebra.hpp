
#ifndef __ALGEBRA_HPP__
#define __ALGEBRA_HPP__

// Reference: https://stackoverflow.com/a/24853652/1377770
// Remedy for disabling the reserved keyword _Complex for MSVC
#ifdef _WIN32
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#pragma warning(disable: 4190)
#include "lapacke.h"
#pragma warning(default: 4190)
#else
#include <complex>
#include "lapacke.h"
#endif // !ifdef _WIN32

#include <vector>

#include "cartosphere/predicates.hpp"

namespace Cartosphere
{
	template <typename Impl>
	class MatrixBase;

	template <typename Impl>
	bool solve(MatrixBase<Impl> const &A,
		std::vector<FLP> &x, std::vector<FLP> const &b);

	// Matrix interface
	template <typename Impl>
	class MatrixBase
	{	
	public:
		// Identifier for different sparse matrix patterns
		enum Pattern
		{
			Dense,
			Tridiagonal,
		};

		// Count the number of slots for non-zero elements in a matrix pattern
		static size_t nnzSlot(size_t row, size_t col, Pattern pattern)
		{
			size_t nnz;
			if (pattern == Tridiagonal)
			{
				size_t min = std::min(row, col);
				size_t max = std::max(row, col);
				nnz = 3 * min - 2 + ((max - min) ? 1 : 0);
			}
			else
			{
				nnz = row * col;
			}
			return nnz;
		}

	public:
		// the greatest possible value for an element of type size_t
		static size_t const npos = -1;

	public:
		// Dense constructor
		MatrixBase(size_t row, size_t col) :
			numberofRows(row), numberofColumns(col), listofEntries(row * col) {}
		// Sparse constructor
		MatrixBase(size_t row, size_t col, size_t nnz) :
			numberofRows(row), numberofColumns(col), listofEntries(nnz) {}
		// Sparse constructor using a pattern
		MatrixBase(size_t row, size_t col, Pattern pattern) :
			numberofRows(row), numberofColumns(col),
			listofEntries(nnzSlots(row, col, pattern)) {}

	public:
		// Access reference to element in row and column
		double &operator()(size_t row, size_t col)
		{
			size_t index = static_cast<Impl *>(this)->index(row, col);
			return listofEntries[index];
		}
		// Access element
		double operator()(size_t row, size_t col) const
		{
			size_t index = static_cast<Impl const *>(this)->index(row, col);
			return listofEntries[index];
		}

	public:
		// Solve single RHS
		bool solve(std::vector<FLP> &x, std::vector<FLP> const &b) const
		{
			return static_cast<Impl const*>(this)->solve(x, b);
		}
		// Return if matrix is sparse
		Pattern pattern() const
		{
			return static_cast<Impl const *>(this)->pattern();
		}
		// Return the number of rows
		size_t rows() const { return numberofRows; }
		// Return the number of rows
		size_t columns() const { return numberofColumns; }
		// Return number of elements
		size_t elements() const { return rows() * columns(); }
		// Friend function for solving equations
		friend bool solve(MatrixBase<Impl> const &A,
			std::vector<FLP> &x, std::vector<FLP> const &b)
		{
			return A.solve(x, b);
		}

	protected:
		// Main list of entries
		std::vector<FLP> listofEntries;
		// dimensions
		size_t numberofRows, numberofColumns;
	};
	
	// Row major matrix interface
	template <typename Impl>
	class RowMajor : public MatrixBase<RowMajor<Impl>>
	{
		typedef MatrixBase<RowMajor<Impl>> Base;

	public:
		// Dense constructor
		RowMajor(size_t rows, size_t cols) : Base(rows, cols) { }
		// Sparse constructor
		RowMajor(size_t rows, size_t cols, size_t nnz) : Base(rows, cols, nnz) { }

	public:
		// Solve single RHS
		bool solve(std::vector<FLP> &x, std::vector<FLP> const &b) const
		{
			return static_cast<Impl const *>(this)->solve(x, b);
		}
		// Return internal indexing
		size_t index(size_t row, size_t col) const
		{
			return static_cast<Impl const *>(this)->index(row, col);
		}
	};

	// Column major matrix interface
	template <typename Impl>
	class ColMajor : public MatrixBase<ColMajor<Impl>>
	{
		typedef MatrixBase<ColMajor<Impl>> Base;

	public:
		// Dense constructor
		ColMajor(size_t rows, size_t cols) : Base(rows, cols) { }
		// Sparse constructor
		ColMajor(size_t rows, size_t cols, size_t nnz) : Base(rows, cols, nnz) { }

	public:
		// Solve single RHS
		bool solve(std::vector<FLP> &x, std::vector<FLP> const &b) const
		{
			return static_cast<Impl const *>(this)->solve(x, b);
		}
		// Return internal indexing
		size_t index(size_t row, size_t col) const
		{
			return static_cast<Impl const *>(this)->index(row, col);
		}
	};

	// Compressed Sparse Row (CSR) Matrix
	class CompressedSparseRowMatrix : public RowMajor<CompressedSparseRowMatrix>
	{
		typedef RowMajor<CompressedSparseRowMatrix> Base;

	public:
		// Sparse Constructor
		CompressedSparseRowMatrix(size_t row, size_t col, size_t nnz) :
			Base(row, col, nnz)
		{
			listofRowCumulative.resize(row + 1);
			listofColumnIndices.resize(nnz);
		}
		// Sparse Constructor using matrix pattern
		CompressedSparseRowMatrix(size_t row, size_t col, Pattern pattern) :
			Base(row, col, pattern)
		{
			listofRowCumulative[0] = 0;
			for (size_t i = col + 2; i < col; ++i)
			if (pattern == Tridiagonal)
			{
				for (size_t i = 1; i <= row; ++i)
				{
					if (i <= col)
					{
						listofRowCumulative[i] = 3 + listofRowCumulative[i - 1];
					}
					if (i == 1 || i == col)
					{
						--listofRowCumulative[i];
					}
					if (i == col + 1)
					{
						--listofRowCumulative[i];
					}
				}
			}
		}

	protected:
		// Auxiliary indexing
		std::vector<size_t> listofRowCumulative;
		std::vector<size_t> listofColumnIndices;
	};

	class DenseColumnMajorMatrix : public ColMajor<DenseColumnMajorMatrix>
	{
		using ColMajor<DenseColumnMajorMatrix>::ColMajor;

	public:
		// Dense constructor
		DenseColumnMajorMatrix(size_t row, size_t col) : ColMajor(row, col) { }

	public:
		// Fill memory
		DenseColumnMajorMatrix &operator=(std::initializer_list<FLP> const &list)
		{
			listofEntries.assign(list);
			listofEntries.resize(elements());
			return *this;
		}
		// Solve single RHS
		bool solve(std::vector<FLP> &x, std::vector<FLP> const &b) const
		{
			// Preparation for the LAPACK dgels interface
			char T = 'N';
			lapack_int M = rows();
			lapack_int N = columns();
			lapack_int nRHS = 1;
			FLP *A = const_cast<FLP *>(listofEntries.data());
			lapack_int ldA = M;
			// Copy b to x because B will be overwritten if INFO = 0
			x = b;
			FLP *B = const_cast<FLP *>(x.data());
			lapack_int ldB = M;
			lapack_int INFO = LAPACKE_dgels(LAPACK_COL_MAJOR,
				T, M, N, nRHS, A, ldA, B, ldB);
			x.resize(N);
			// INFO == 0 means successful exit
			return INFO == 0;
		}
		// Return internal indexing
		size_t index(size_t row, size_t col) const
		{
			return row + col * columns();
		}
	};

	typedef CompressedSparseRowMatrix CSR;
	typedef DenseColumnMajorMatrix Matrix;
}

#endif //!__ALGEBRA_HPP
