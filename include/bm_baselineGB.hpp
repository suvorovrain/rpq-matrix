//
// Created by Adrián on 16/4/24.
//
#include <cstdlib>
#include <iostream>
#include <string>
extern "C"
{
#include "../../../GraphBLAS/Include/GraphBLAS.h"
#include "../../../LAGraph/include/LAGraph.h"
}

#ifndef RPQ_MATRIX_BASELINEGB_HPP
#define RPQ_MATRIX_BASELINEGB_HPP
#define fullSide (~(uint)0) 

#define mmax(x, y) (((x) > (y)) ? (x) : (y))
namespace bm_baselinegb
{
    extern "C"
    {
#include "baseline/utilstime.h"
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
            printf("unreachable\n");
            exit(1);
            GrB_Matrix E;
            GrB_Matrix_new(&E, GrB_BOOL, side, side);
          
            for (GrB_Index i = 0; i < side; i++)
            {
                GrB_Matrix_setElement_BOOL(E, true, i, i);
            }

            // GrB_Index nvals;
            // GrB_Matrix_nvals(&nvals, E);
            // printf("%d E\n", nvals);
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

            if (info != GrB_SUCCESS)
            {
                std::cerr << "Error occured while reading the matrix: " << info << std::endl;
                return NULL;
            }
            GrB_Matrix_wait(A, GrB_COMPLETE);
            GrB_Index nrows, ncols;
            GrB_Matrix_nrows(&nrows, A);
            GrB_Matrix_ncols(&ncols, A);

            GrB_Matrix_resize(A, nrows + 1, ncols + 1);
            return A;
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
            GrB_Matrix C;
            GrB_Index nrows, ncols;

            GrB_Matrix_nrows(&nrows, A);
            GrB_Matrix_ncols(&ncols, A);

            GrB_Matrix_new(&C, GrB_BOOL, nrows, ncols);

            GrB_Matrix_eWiseAdd_BinaryOp(C, NULL, NULL, GrB_LOR, A, B, NULL);

            return C;
        }

        // version with one row or one column, or both
        static inline GrB_Matrix sum1(uint64_t row, GrB_Matrix A, GrB_Matrix B, uint64_t col)
        {
            // printf("DEBUG: sum1\n");
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
                // printf("DEBUG: sum1 !col\n");
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
                // GrB_Index nvalsA, nvalsB;
                // GrB_Matrix_nvals(&nvalsA,A);
                // GrB_Matrix_nvals(&nvalsB,B);

                // printf("DEBUG: sum1 !row\n");
                //         printf("DEBUG: A elems = %d and B elems = %d\n", nvalsA,nvalsB);

                GrB_Vector rowA, rowB, rowC;
                GrB_Vector_new(&rowA, GrB_BOOL, ncols);
                GrB_Vector_new(&rowB, GrB_BOOL, ncols);
                GrB_Vector_new(&rowC, GrB_BOOL, ncols);
                GrB_Col_extract(rowA, NULL, NULL, A, GrB_ALL, ncols, row-1, GrB_DESC_T0);
                GrB_Col_extract(rowB, NULL, NULL, B, GrB_ALL, ncols, row-1, GrB_DESC_T0);

                GrB_Vector_eWiseAdd_BinaryOp(rowC, NULL, NULL, GrB_LOR, rowA, rowB, NULL);
                GrB_Row_assign(C, NULL, NULL, rowC, row-1, GrB_ALL, ncols, NULL);
                GrB_Vector_free(&rowA);
                GrB_Vector_free(&rowB);
                GrB_Vector_free(&rowC);
                return C;
            }
        };

        // (boolean) product of two matrices, assumed to be of the same side
        // only rowA of A and colB of B are considered if not fullSide
        static inline GrB_Matrix mult(GrB_Matrix A, GrB_Matrix B)
        {
            // GrB_Index nvalsA, nvalsB;
            // GrB_Matrix_nvals(&nvalsA, A);
            // GrB_Matrix_nvals(&nvalsB, B);
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
            // printf("DEBUG: mult1\n");
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
                GrB_Col_extract(rowA, GrB_NULL, GrB_NULL, A, GrB_ALL, ncolsA, row - 1, GrB_DESC_T0);
                GrB_Row_assign(tmpA, GrB_NULL, GrB_NULL, rowA, row - 1, GrB_ALL, ncolsA, GrB_NULL);

                C = mult(tmpA, tmpB);
                destroy(tmpA);
                destroy(tmpB);
                GrB_Vector_free(&colB);
                GrB_Vector_free(&rowA);
                // destroy(AT);
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
                GrB_Col_extract(rowA, GrB_NULL, GrB_NULL, A, GrB_ALL, ncolsA, row - 1, GrB_DESC_T0);
                GrB_Row_assign(tmpA, GrB_NULL, GrB_NULL, rowA, row - 1, GrB_ALL, ncolsA, GrB_NULL);

                C = mult(tmpA, B);
                destroy(tmpA);
                GrB_Vector_free(&rowA);
                // destroy(AT);
                return C;
            }
        }

        static inline GrB_Matrix pow(GrB_Matrix A)
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

                GrB_Matrix temp;
                GrB_Matrix_new(&temp, GrB_BOOL, nrows, ncols);

                temp = mult(result, A);

                result = sum(result, temp);

                destroy(temp);

                GrB_Matrix_nvals(&prev_nvals, prev);
                GrB_Matrix_nvals(&result_nvals, result);
                changed = result_nvals != prev_nvals;

                GrB_Matrix_dup(&prev, result);
            } while (changed);
            destroy(prev);
            GrB_Matrix_free(&prev);
            return result;
        }

        // transitive closure of a matrix, pos says if it's + rather than *
        static inline GrB_Matrix clos(GrB_Matrix A, uint pos)
        {

            // std::cout << "unreachable" << std::endl;
            // exit(1);
            GrB_Index nrows, ncols;
            GrB_Matrix_nrows(&nrows, A);
            GrB_Matrix_ncols(&ncols, A);
            GrB_Matrix C = pow(A);

            if (!pos)
            {
                GrB_Index dim = (nrows < ncols) ? ncols : nrows;
                for (GrB_Index i = 0; i < dim; i++)
                {
                    if (i < nrows && i < ncols)
                    {
                        GrB_Matrix_setElement_BOOL(C, 1, i, i);
                    }
                }
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
                // printf("DEBUG: clos_row ID null\n");
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
                    // printf("DEBUG: clow_row !ID pos\n");
                    S = mult1(row, ID, A, full_side);
                    // GrB_Matrix_nvals(&elems, S);
                }

                else
                {
                    // printf("DEBUG: clow_row !ID !pos\n");
                    // GrB_Index nvalsID;
                    // GrB_Matrix_nvals(&nvalsID,ID);
                    // printf("DEBUG:  ID elems in closrow: %d\n", nvalsID);
                    E = empty(dim, dim);
                    S = sum1(row, ID, E, full_side);
                    destroy(E);
                }
            }
            // GrB_Matrix_nvals(&elems, S);
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
            // printf("DEBUG: clos_row before return\n");
            return S;
        }

        static inline GrB_Matrix clos1(uint64_t row, GrB_Matrix A, uint pos, uint64_t col)
        {
            GrB_Index nrow, ncol;
            // printf("DEBUG: clos1\n");
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
            GrB_Col_extract(col_vec, NULL, NULL, A, GrB_ALL, side, col - 1, GrB_DESC_T0);

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
                        // printf("DEBUG: mult_clos1\n");
            // GrB_Index nvalsA;
            // GrB_Matrix_nvals(&nvalsA,A);
                //  printf("DEBUG:  A elems in mult_clos1: %d\n", nvalsA);

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
                    // printf("DEBUG: mult_clos !row\n");
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
                        // printf("DEBUG: clos_mult1\n");

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
