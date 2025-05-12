//
// Created by Adrián on 16/4/24.
//
#include <cstdlib>
#include <iostream>
#include <string>
extern "C"
{
#include "../../../GraphBLAS/Include/GraphBLAS.h" // here put GB matirx
#include "../../../LAGraph/include/LAGraph.h"
}

#ifndef RPQ_MATRIX_BASELINEGB_HPP
#define RPQ_MATRIX_BASELINEGB_HPP
#define fullSide (~(uint)0) // (*** got this from original matrix.h include. Maybe remove later idk ***)

#define mmax(x, y) (((x) > (y)) ? (x) : (y))
namespace bm_baselinegb
{
    extern "C"
    {
#include "baseline/utilstime.h" // here all good
    }

    class wrapper
    {
        // TODO: add error handling. and think about memory leaks
    public:
        static const uint64_t full_side = fullSide;
        typedef GrB_Matrix matrix_type;   // (*** in original this is a pointer to matrix ***)
        typedef GrB_Matrix s_matrix_type; // (*** and this is matrix itself ***)

        static void time_begin()
        {
            time_beg();
        }

        // destroys matrix
        static inline void destroy(GrB_Matrix M)
        {
            GrB_Matrix_free(&M);
        };

        // creates an empty matrix
        static inline GrB_Matrix empty(uint64_t height, uint64_t width)
        {
            GrB_Matrix A;
            GrB_Matrix_new(&A, GrB_BOOL, height, width);
            return A;
        };

        // creates an identity matrix
        static inline GrB_Matrix id(uint64_t side)
        {
            GrB_Matrix E;
            GrB_Matrix_new(&E, GrB_BOOL, side, side);
            GrB_Index indices[side];
            bool values[side];
            for (GrB_Index i = 0; i < side; i++)
            {
                indices[i] = i;
                values[i] = true;
            }
            GrB_Matrix_build_BOOL(E, indices, indices, values, side, GrB_FIRST_BOOL);
            return E;
        };

        // version of identity matrix with one row or column
        static inline matrix_type id1(uint64_t side, uint64_t rc)
        {
            GrB_Matrix E;
            GrB_Matrix_new(&E, GrB_BOOL, side, side);
            GrB_Matrix_setElement_BOOL(E, true, rc - 1, rc - 1);
            return E;
        };

        static inline matrix_type mat_one(uint64_t height, uint64_t width, uint64_t r, uint64_t c)
        {
            GrB_Matrix E;
            GrB_Matrix_new(&E, GrB_BOOL, height, width);
            GrB_Matrix_setElement_BOOL(E, true, r - 1, c - 1);
            return E;
        };

        // transpose a matrix, creating a non-allocated copy that shares the
        // data. You need not (and should not) matDestroy this copy
        // use *M = matTranspose(M) to actually transpose M
        static inline GrB_Matrix transpose(GrB_Matrix M)
        {
            GrB_Matrix MT;
            GrB_Index nrows, ncols;
            GrB_Matrix_nrows(&nrows, M);
            GrB_Matrix_ncols(&ncols, M);
            GrB_Matrix_new(&MT, GrB_BOOL, nrows, ncols);
            GrB_transpose(MT, GrB_NULL, GrB_NULL, M, GrB_NULL);
            return MT;
        };

        // writes M to file, which must be opened for writing
        static inline void save(GrB_Matrix M, FILE *file)
        {
            GxB_Matrix_fprint(M, "M", GxB_COMPLETE, file);
        };

        // loads matrix from file, which must be opened for reading
        static inline GrB_Matrix load(FILE *file)
        {
            GrB_Matrix A = NULL;
            GrB_Info info = LAGraph_MMRead(&A, file, NULL);
            GrB_Matrix_wait(A, GrB_COMPLETE);
            if (info != GrB_SUCCESS)
            {
                std::cerr << "Error occured while reading the matrix: " << info << std::endl;
                return NULL;
            }

            GrB_Index nrows, ncols;
            GrB_Matrix_nrows(&nrows, A);
            GrB_Matrix_ncols(&ncols, A);

            GrB_Matrix B = NULL;
            GrB_Matrix_new(&B, GrB_BOOL, nrows + 1, ncols + 1);

            for (GrB_Index i = 0; i < nrows; ++i)
            {
                for (GrB_Index j = 0; j < ncols; ++j)
                {
                    bool val;
                    if (GrB_Matrix_extractElement_BOOL(&val, A, i, j) == GrB_SUCCESS)
                    {
                        GrB_Matrix_setElement_BOOL(B, val, i, j);
                    }
                }
            }

            GrB_Matrix_free(&A);
            return B;
        }

        // space of the matrix, in w-bit words
        static inline uint64_t space(GrB_Matrix M)
        {
            uint64_t memory_usage;
            return GxB_Matrix_memoryUsage(&memory_usage, M);
            return memory_usage;
        };

        // (boolean) sum of two matrices, assumed to be of the same side
        static inline GrB_Matrix sum(GrB_Matrix A, GrB_Matrix B)
        {
            GrB_Index nrowsA, ncolsA, nrowsB, ncolsB, nvalsA, nvalsB;
            GrB_Matrix_nrows(&nrowsA, A);
            GrB_Matrix_ncols(&ncolsA, A);
            GrB_Matrix_nrows(&nrowsB, B);
            GrB_Matrix_ncols(&ncolsB, B);

            GrB_Index nrows = std::max(nrowsA, nrowsB);
            GrB_Index ncols = std::max(ncolsA, ncolsB);

            GrB_Matrix Atmp, Btmp;
            GrB_Matrix_new(&Atmp, GrB_BOOL, nrows, ncols);
            GrB_Matrix_new(&Btmp, GrB_BOOL, nrows, ncols);

            GrB_Matrix_assign(Atmp, NULL, NULL, A, GrB_ALL, nrowsA, GrB_ALL, ncolsA, NULL);
            GrB_Matrix_assign(Btmp, NULL, NULL, B, GrB_ALL, nrowsB, GrB_ALL, ncolsB, NULL);

            GrB_Matrix C;
            GrB_Matrix_new(&C, GrB_BOOL, nrows, ncols);
            GrB_Matrix_eWiseAdd_BinaryOp(C, NULL, NULL, GrB_LOR, Atmp, Btmp, NULL);

            GrB_Matrix_free(&Atmp);
            GrB_Matrix_free(&Btmp);

            return C;
        }

        // version with one row or one column, or both
        static inline GrB_Matrix sum1(uint64_t row, GrB_Matrix A, GrB_Matrix B, uint64_t col)
        {
            GrB_Index nrows, ncols;
            GrB_Matrix C;
            GrB_Matrix_nrows(&nrows, A);
            GrB_Matrix_ncols(&ncols, A);
            GrB_Matrix_new(&C, GrB_BOOL, nrows, ncols);
            if (row == full_side && col == full_side)
            {
                return sum(A, B);
            }
            else if (row != full_side && col != full_side)
            {
                bool elA = false, elB = false;
                GrB_Matrix_extractElement_BOOL(&elA, A, row - 1, col - 1);
                GrB_Matrix_extractElement_BOOL(&elB, B, row - 1, col - 1);
                GrB_Matrix_setElement_BOOL(C, elA || elB, row - 1, col - 1);
                return C;
            }
            else if (row == full_side && col != full_side)
            {
                GrB_Vector colA, colB, colC;
                GrB_Vector_new(&colA, GrB_BOOL, nrows);
                GrB_Vector_new(&colB, GrB_BOOL, nrows);
                GrB_Vector_new(&colC, GrB_BOOL, nrows);
                GrB_Col_extract(colA, NULL, NULL, A, GrB_ALL, nrows, col - 1, NULL);
                GrB_Col_extract(colB, NULL, NULL, B, GrB_ALL, nrows, col - 1, NULL);

                GrB_Vector_eWiseAdd_BinaryOp(colC, NULL, NULL, GrB_LOR, colA, colB, NULL);
                GrB_Col_assign(C, NULL, NULL, colC, GrB_ALL, nrows, col - 1, NULL);
                GrB_Vector_free(&colA);
                GrB_Vector_free(&colB);
                GrB_Vector_free(&colC);
                return C;
            }
            else

            {
                GrB_Matrix AT, BT, CT;
                AT = transpose(A);

                BT = transpose(B);
                C = sum1(col, AT, BT, row);
                CT = transpose(C);
                destroy(C);
                destroy(AT);
                destroy(BT);
                return CT;
            }
        };

        // (boolean) product of two matrices, assumed to be of the same side
        // only rowA of A and colB of B are considered if not fullSide
        static inline GrB_Matrix mult(GrB_Matrix A, GrB_Matrix B)
        {
            GrB_Index nvalsA, nvalsB;
            GrB_Matrix_nvals(&nvalsA, A);
            GrB_Matrix_nvals(&nvalsB, B);
            GrB_Index nrows, ncols;
            GrB_Matrix C;
            GrB_Matrix_nrows(&nrows, A);
            GrB_Matrix_ncols(&ncols, B);
            GrB_Matrix_new(&C, GrB_BOOL, nrows, ncols);
            GrB_mxm(C, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, A, B, GrB_NULL);
            return C;
        };

        // version with one row or one column, or both
        static inline GrB_Matrix mult1(uint64_t row, GrB_Matrix A, GrB_Matrix B, uint64_t col)
        {
            GrB_Index nrowsA, ncolsA, nrowsB, ncolsB;
            GrB_Matrix C;
            GrB_Matrix_nrows(&nrowsA, A);
            GrB_Matrix_ncols(&ncolsA, A);
            GrB_Matrix_nrows(&nrowsB, B);
            GrB_Matrix_ncols(&ncolsB, B);
            GrB_Matrix_new(&C, GrB_BOOL, nrowsA, ncolsB);

            if (row == full_side && col == full_side)
            {
                return mult(A, B);
            }
            else if (row != full_side && col != full_side)
            {
                // <r>(A x B)<c> = (<r>A) x (B<c>)
                GrB_Matrix tmpA, tmpB;
                GrB_Matrix_new(&tmpA, GrB_BOOL, nrowsA, ncolsA);
                GrB_Matrix_new(&tmpB, GrB_BOOL, nrowsB, ncolsB);

                // rxtract column c of B
                GrB_Vector colB;
                GrB_Vector_new(&colB, GrB_BOOL, nrowsB);
                GrB_Col_extract(colB, GrB_NULL, GrB_NULL, B, GrB_ALL, nrowsB, col - 1, GrB_NULL);
                GrB_Col_assign(tmpB, GrB_NULL, GrB_NULL, colB, GrB_ALL, nrowsB, col - 1, GrB_NULL);

                // extract row r of A as vector via transpose
                GrB_Vector rowA;
                GrB_Vector_new(&rowA, GrB_BOOL, ncolsA);
                GrB_Matrix AT;
                GrB_Matrix_new(&AT, GrB_BOOL, ncolsA, nrowsA);
                GrB_transpose(AT, GrB_NULL, GrB_NULL, A, GrB_NULL);
                GrB_Col_extract(rowA, GrB_NULL, GrB_NULL, AT, GrB_ALL, ncolsA, row - 1, GrB_NULL);
                GrB_Row_assign(tmpA, GrB_NULL, GrB_NULL, rowA, row - 1, GrB_ALL, ncolsA, GrB_NULL);

                C = mult(tmpA, tmpB);
                destroy(tmpA);
                destroy(tmpB);
                GrB_Vector_free(&colB);
                GrB_Vector_free(&rowA);
                destroy(AT);
                return C;
            }
            else if (row == full_side && col != full_side)
            {
                // (A x B)<c> = A x (B<c>)
                GrB_Matrix tmpB;
                GrB_Matrix_new(&tmpB, GrB_BOOL, nrowsB, ncolsB);

                GrB_Vector colB;
                GrB_Vector_new(&colB, GrB_BOOL, nrowsB);
                GrB_Col_extract(colB, GrB_NULL, GrB_NULL, B, GrB_ALL, nrowsB, col - 1, GrB_NULL);
                GrB_Col_assign(tmpB, GrB_NULL, GrB_NULL, colB, GrB_ALL, nrowsB, col - 1, GrB_NULL);

                C = mult(A, tmpB);
                destroy(tmpB);
                GrB_Vector_free(&colB);
                return C;
            }
            else
            {
                // <r>(A x B) = (<r>A) x B
                GrB_Matrix tmpA;
                GrB_Matrix_new(&tmpA, GrB_BOOL, nrowsA, ncolsA);

                GrB_Vector rowA;
                GrB_Vector_new(&rowA, GrB_BOOL, ncolsA);
                GrB_Matrix AT;
                GrB_Matrix_new(&AT, GrB_BOOL, ncolsA, nrowsA);
                GrB_transpose(AT, GrB_NULL, GrB_NULL, A, GrB_NULL);
                GrB_Col_extract(rowA, GrB_NULL, GrB_NULL, AT, GrB_ALL, ncolsA, row - 1, GrB_NULL);
                GrB_Row_assign(tmpA, GrB_NULL, GrB_NULL, rowA, row - 1, GrB_ALL, ncolsA, GrB_NULL);

                C = mult(tmpA, B);
                destroy(tmpA);
                GrB_Vector_free(&rowA);
                destroy(AT);
                return C;
            }
        }

        static inline GrB_Matrix pow(GrB_Matrix A, uint64_t row, uint64_t col)
        {
            GrB_Index nrows, ncols;
            GrB_Matrix_nrows(&nrows, A);
            GrB_Matrix_ncols(&ncols, A);

            GrB_Matrix result;
            GrB_Matrix_new(&result, GrB_BOOL, nrows, ncols);
            GrB_Matrix_dup(&result, A);

            GrB_Matrix prev;
            GrB_Matrix_new(&prev, GrB_BOOL, nrows, ncols);

            GrB_Index prev_nvals = 0, result_nvals = 0;
            bool changed;
            do
            {

                GrB_Matrix_dup(&prev, result);

                GrB_Matrix temp;
                GrB_Matrix_new(&temp, GrB_BOOL, nrows, ncols);
                temp = mult1(row, result, A, col);

                result = sum1(row, result, temp, col);

                destroy(temp);

                GrB_Matrix_nvals(&prev_nvals, prev);
                GrB_Matrix_nvals(&result_nvals, result);

                changed = result_nvals != prev_nvals;

                destroy(prev);
                GrB_Matrix_dup(&prev, result);

            } while (changed);

            GrB_Matrix_free(&prev);
            return result;
        }

        // transitive closure of a matrix, pos says if it's + rather than *
        static inline GrB_Matrix clos(GrB_Matrix A, uint pos)
        {

            GrB_Index nrows, ncols;
            GrB_Matrix C = pow(A, full_side, full_side);

            if (!pos)
            {
                GrB_Matrix E = id(mmax(nrows, ncols));
                GrB_Matrix summ;
                GrB_Matrix_new(&summ, GrB_BOOL, nrows, ncols);
                summ = sum(C, E);
                free(C);
                free(E);
                C = summ;
            }
            return C;
        }

        static inline GrB_Matrix clos_row(uint row, GrB_Matrix ID, GrB_Matrix A, uint pos, uint *coltest)
        {
            GrB_Matrix M, P, S, E;

            GrB_Index elems;

            GrB_Index nrowsA, ncolsA;
            GrB_Matrix_nrows(&nrowsA, A);
            GrB_Matrix_ncols(&ncolsA, A);
            uint dim = mmax(nrowsA, ncolsA);
            if (ID == NULL)
            {
                if (pos)
                    E = empty(dim, dim);
                else
                {
                    E = mat_one(dim, dim, row, row);
                }
                S = sum1(row, A, E, full_side);
                destroy(E);
            }
            else
            {
                if (pos)
                {

                    S = mult1(row, ID, A, full_side);
                    GrB_Matrix_nvals(&elems, S);
                }

                else
                {
                    E = empty(dim, dim);
                    S = sum1(row, ID, E, full_side);
                    destroy(E);
                }
            }
            GrB_Matrix_nvals(&elems, S);
            bool accessS;
            if (coltest)
                GrB_Matrix_extractElement_BOOL(&accessS, S, row, *coltest);
            if (coltest && accessS)
            {
                destroy(S);
                *coltest = 1;
                return NULL;
            }
            P = mult(S, A);
            M = S;
            S = sum(P, S);

            GrB_Index elems_new;
            GrB_Matrix_nvals(&elems_new, S);

            while (elems_new != elems)
            {
                if (coltest)
                    GrB_Matrix_extractElement_BOOL(&accessS, S, row, *coltest);
                if (coltest && accessS)
                {
                    destroy(S);
                    destroy(P);
                    *coltest = 1;
                    return NULL;
                }
                elems = elems_new;
                M = P;
                P = mult(P, A);
                destroy(M);
                M = S;
                S = sum(S, P);
                destroy(M);
                GrB_Matrix_nvals(&elems_new, S);
            }
            destroy(P);
            if (coltest)
            {
                destroy(S);
                *coltest = 0;
                return NULL;
            }
            return S;
        }

        static inline GrB_Matrix clos1(uint64_t row, GrB_Matrix A, uint pos, uint64_t col)
        {
            GrB_Index nrow, ncol;

            GrB_Index nrowsA, ncolsA;
            GrB_Matrix_nrows(&nrowsA, A);
            GrB_Matrix_ncols(&ncolsA, A);
            uint side = mmax(nrowsA, ncolsA);
            uint test;
            GrB_Matrix M;
            GrB_Matrix At;

            if (row == full_side)
            {
                if (col == full_side)
                    return clos(A, pos);
                else
                {
                    GrB_Index n;
                    At = transpose(A);
                    M = clos_row(col, NULL, At, pos, NULL);
                    M = transpose(M);
                    return M;
                }
            }
            else
            {
                if (col == full_side)
                    return clos_row(row, NULL, A, pos, NULL);
            }

            GrB_Vector row_vec;
            GrB_Vector_new(&row_vec, GrB_BOOL, side);
            GrB_Col_extract(row_vec, NULL, NULL, A, GrB_ALL, side, row - 1, NULL);

            GrB_Vector_nvals(&ncol, row_vec);

            GrB_Vector col_vec;
            GrB_Vector_new(&col_vec, GrB_BOOL, side);
            GrB_Matrix AT;
            AT = transpose(A);
            GrB_Col_extract(col_vec, NULL, NULL, AT, GrB_ALL, side, col - 1, NULL);

            if (ncol < nrow)
            {
                test = row;
                At = transpose(A);
                clos_row(col, NULL, At, pos, &test);
            }
            else
            {
                test = col;
                clos_row(row, NULL, A, pos, &test);
            }
            if (test)
            {
                return mat_one(side, side, row, col);
            }
            else
                return empty(side, side);
        }

        // computes [row] A B* [col] (pos=0) or [row] A B+ [col] (pos=1)
        static inline GrB_Matrix mult_clos1(uint64_t row, GrB_Matrix A, GrB_Matrix B, uint pos, uint64_t col)
        {
            GrB_Index nrowsA, ncolsA, nrowsB, ncolsB;
            GrB_Matrix_nrows(&nrowsA, A);
            GrB_Matrix_ncols(&ncolsA, A);
            GrB_Matrix_nrows(&nrowsB, B);
            GrB_Matrix_ncols(&ncolsB, B);
            uint64_t side = mmax(nrowsA, ncolsA);

            if (row == full_side)
            {
                if (col == full_side)
                {
                    GrB_Matrix M1 = clos(B, pos);
                    GrB_Matrix M = mult(A, M1);
                    GrB_Matrix_free(&M1);
                    return M;
                }
                else
                {

                    GrB_Matrix M1 = clos1(full_side, B, pos, col);
                    GrB_Matrix M = mult(A, M1);
                    GrB_Matrix_free(&M1);
                    return M;
                }
            }
            else
            {
                if (col == full_side)
                {
                    return clos_row(row, A, B, pos, NULL);
                }
            }
            uint32_t test = col;
            clos_row(row, A, B, pos, &test);
            if (test)
                return mat_one(side, side, row, col);
            return empty(side, side);
        }

        // computes [row] A* B [col] (pos=0) or [row] A+ B [col] (pos=1)
        static inline GrB_Matrix clos_mult1(uint64_t row, GrB_Matrix A, uint pos, GrB_Matrix B, uint64_t col)
        {
            GrB_Matrix At, Bt;
            At = transpose(A);
            Bt = transpose(B);

            GrB_Matrix M = mult_clos1(col, Bt, At, pos, row);
            GrB_Matrix M_transposed = transpose(M);

            GrB_Matrix_free(&At);
            GrB_Matrix_free(&Bt);
            GrB_Matrix_free(&M);

            return M_transposed;
        }
    };
}

#endif // RPQ_MATRIX_BASELINEGB_HPP
