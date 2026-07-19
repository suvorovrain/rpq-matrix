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
    
    public:
        struct matrix_pair
        {
            GrB_Matrix Acsr;
            GrB_Matrix Acsc;
        };
        static const uint64_t full_side = fullSide;
        typedef GrB_Matrix matrix_type;
        typedef GrB_Matrix s_matrix_type;

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

        static inline GrB_Matrix id(uint64_t side)
        {
            GrB_Matrix E;
            GrB_Matrix_new(&E, GrB_BOOL, side, side);

            GrB_Vector v;
            GrB_Vector_new(&v, GrB_BOOL, side);
            GrB_Vector_assign_BOOL(v, NULL, NULL, true, GrB_ALL, side, NULL);

            GrB_Matrix_diag(&E, v, 0);

            GrB_Vector_free(&v);

            return E;
        }

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
        static inline matrix_pair load(FILE *file)
        {
            GrB_Matrix A = NULL;
            LAGraph_MMRead(&A, file, NULL);
            GrB_Matrix_wait(A, GrB_MATERIALIZE);

            GrB_Index nrows, ncols;
            GrB_Matrix_nrows(&nrows, A);
            GrB_Matrix_ncols(&ncols, A);
            GrB_Matrix B;

            GrB_Matrix_resize(A, nrows + 1, ncols + 1);
            GrB_Matrix_dup(&B, A);
            GrB_Matrix_wait(A, GrB_MATERIALIZE);

            struct matrix_pair Acsccsr;
            Acsccsr.Acsc = B;
            Acsccsr.Acsr = A;
            GxB_Matrix_Option_set(B, GxB_FORMAT, GxB_BY_COL);
            GrB_Matrix_wait(A, GrB_MATERIALIZE);
            return Acsccsr;
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

            GrB_Index nvalsA, nvalsB;
            GrB_Matrix_nvals(&nvalsA, A);
            GrB_Matrix_nvals(&nvalsB, B);

            if (nvalsA == 0 || nvalsB == 0)
            {
                GrB_Matrix tmp = (nvalsA == 0) ? B : A;
                GrB_Matrix_dup(&C, tmp);
                return C;
            }

            GrB_Matrix_eWiseAdd_BinaryOp(C, NULL, NULL, GxB_ANY_BOOL, A, B, NULL);
            return C;
        }

        static inline GrB_Matrix sum1(uint64_t row, GrB_Matrix A, GrB_Matrix B, uint64_t col)
        {
            GrB_Index nrows, ncols;
            GrB_Matrix_nrows(&nrows, A);
            GrB_Matrix_ncols(&ncols, A);

            if (row == full_side && col == full_side)
            {
                return sum(A, B);
            }

            GrB_Matrix C;
            GrB_Matrix_new(&C, GrB_BOOL, nrows, ncols);

            if (row != full_side && col != full_side)
            {
                // <r>(A x B)<c> = (<r>A<c>) x (<r>B<c>)
                bool elA = false, elB = false;
                GrB_Matrix_extractElement_BOOL(&elA, A, row - 1, col - 1);
                GrB_Matrix_extractElement_BOOL(&elB, B, row - 1, col - 1);
                GrB_Matrix_setElement_BOOL(C, elA || elB, row - 1, col - 1);
                return C;
            }

            if (row == full_side)
            {
                // (A x B)<c> = (A<c>) x (B<c>)
                GrB_Vector colA, colB, colC;
                GrB_Vector_new(&colA, GrB_BOOL, nrows);
                GrB_Vector_new(&colB, GrB_BOOL, nrows);
                GrB_Vector_new(&colC, GrB_BOOL, nrows);

                GrB_Col_extract(colA, NULL, NULL, A, GrB_ALL, nrows, col - 1, NULL);
                GrB_Col_extract(colB, NULL, NULL, B, GrB_ALL, nrows, col - 1, NULL);

                GrB_Index nvalsVA, nvalsVB;
                GrB_Vector_nvals(&nvalsVA, colA);
                GrB_Vector_nvals(&nvalsVB, colB);

                if (nvalsVA == 0 || nvalsVB == 0)
                {
                    GrB_Vector tmp = (nvalsVA == 0) ? colB : colA;
                    GrB_Col_assign(C, NULL, NULL, tmp, GrB_ALL, nrows, col - 1, NULL);
                    GrB_Vector_free(&colA);
                    GrB_Vector_free(&colB);
                    GrB_Vector_free(&colC);
                    return C;
                }

                GrB_Vector_eWiseAdd_BinaryOp(colC, NULL, NULL, GxB_ANY_BOOL, colA, colB, NULL);
                GrB_Col_assign(C, NULL, NULL, colC, GrB_ALL, nrows, col - 1, NULL);
                GrB_Vector_free(&colA);
                GrB_Vector_free(&colB);
                GrB_Vector_free(&colC);
                return C;
            }
            else
            {
                // <r>(A x B) = (<r>A) x (<r>B)
                GrB_Vector rowA, rowB, rowC;
                GrB_Vector_new(&rowA, GrB_BOOL, ncols);
                GrB_Vector_new(&rowB, GrB_BOOL, ncols);
                GrB_Vector_new(&rowC, GrB_BOOL, ncols);

                GrB_Col_extract(rowA, NULL, NULL, A, GrB_ALL, ncols, row - 1, GrB_DESC_T0);
                GrB_Col_extract(rowB, NULL, NULL, B, GrB_ALL, ncols, row - 1, GrB_DESC_T0);

                GrB_Index nvalsVA, nvalsVB;
                GrB_Vector_nvals(&nvalsVA, rowA);
                GrB_Vector_nvals(&nvalsVB, rowB);

                if (nvalsVA == 0 || nvalsVB == 0)
                {
                    GrB_Vector tmp = (nvalsVA == 0) ? rowB : rowA;
                    GrB_Row_assign(C, NULL, NULL, tmp, row - 1, GrB_ALL, ncols, NULL);
                    GrB_Vector_free(&rowA);
                    GrB_Vector_free(&rowB);
                    GrB_Vector_free(&rowC);
                    return C;
                }

                GrB_Vector_eWiseAdd_BinaryOp(rowC, NULL, NULL, GxB_ANY_BOOL, rowA, rowB, NULL);
                GrB_Row_assign(C, NULL, NULL, rowC, row - 1, GrB_ALL, ncols, NULL);
                GrB_Vector_free(&rowA);
                GrB_Vector_free(&rowB);
                GrB_Vector_free(&rowC);
                return C;
            }
        }

        // (boolean) product of two matrices, assumed to be of the same side
        static inline GrB_Matrix mult(GrB_Matrix A, GrB_Matrix B)
        {
            GrB_Index nrows, ncols;
            GrB_Matrix_nrows(&nrows, A);
            GrB_Matrix_ncols(&ncols, B);

            GrB_Index nvalsA, nvalsB;
            GrB_Matrix_nvals(&nvalsA, A);
            GrB_Matrix_nvals(&nvalsB, B);

            if (nvalsA == 0 || nvalsB == 0)
                return empty(nrows, ncols);

            GrB_Matrix C;
            GrB_Matrix_new(&C, GrB_BOOL, nrows, ncols);
            GrB_mxm(C, GrB_NULL, GrB_NULL, GxB_ANY_PAIR_BOOL, A, B, GrB_NULL);
            return C;
        }

        static inline GrB_Matrix mult1(uint64_t row, GrB_Matrix A, GrB_Matrix B, uint64_t col)
        {
            GrB_Index nrowsA, ncolsA, nrowsB, ncolsB;
            GrB_Matrix_nrows(&nrowsA, A);
            GrB_Matrix_ncols(&ncolsA, A);
            GrB_Matrix_nrows(&nrowsB, B);
            GrB_Matrix_ncols(&ncolsB, B);

            if (row == full_side && col == full_side)
                return mult(A, B);

            GrB_Matrix C;
            GrB_Matrix_new(&C, GrB_BOOL, nrowsA, ncolsB);
            
            if (row != full_side && col != full_side)
            {
                // <r>(A x B)<c> = (<r>A) x (B<c>)
                GrB_Matrix tmpA, tmpB;
                GrB_Vector colB, rowA;

                GrB_Matrix_new(&tmpA, GrB_BOOL, nrowsA, ncolsA);
                GrB_Matrix_new(&tmpB, GrB_BOOL, nrowsB, ncolsB);
                GrB_Vector_new(&colB, GrB_BOOL, nrowsB);
                GrB_Vector_new(&rowA, GrB_BOOL, ncolsA);

                GrB_Col_extract(colB, GrB_NULL, GrB_NULL, B, GrB_ALL, nrowsB, col - 1, GrB_NULL);
                GrB_Col_assign(tmpB, GrB_NULL, GrB_NULL, colB, GrB_ALL, nrowsB, col - 1, GrB_NULL);
                GrB_Col_extract(rowA, GrB_NULL, GrB_NULL, A, GrB_ALL, ncolsA, row - 1, GrB_DESC_T0);
                GrB_Row_assign(tmpA, GrB_NULL, GrB_NULL, rowA, row - 1, GrB_ALL, ncolsA, GrB_NULL);

                GrB_Matrix_free(&C);
                C = mult(tmpA, tmpB);
                GrB_Matrix_free(&tmpA);
                GrB_Matrix_free(&tmpB);
                GrB_Vector_free(&colB);
                GrB_Vector_free(&rowA);
            }
            else if (row == full_side)
            {
                // (A x B)<c> = A x (B<c>)
                GrB_Matrix tmpB;
                GrB_Vector colB;

                GrB_Matrix_new(&tmpB, GrB_BOOL, nrowsB, ncolsB);
                GrB_Vector_new(&colB, GrB_BOOL, nrowsB);

                GrB_Col_extract(colB, GrB_NULL, GrB_NULL, B, GrB_ALL, nrowsB, col - 1, GrB_NULL);
                GrB_Col_assign(tmpB, GrB_NULL, GrB_NULL, colB, GrB_ALL, nrowsB, col - 1, GrB_NULL);

                GrB_Matrix_free(&C);
                C = mult(A, tmpB);
                GrB_Matrix_free(&tmpB);
                GrB_Vector_free(&colB);
            }
            else
            {
                // (A x B)<c> = (<r>A) x B
                GrB_Matrix tmpA;
                GrB_Vector rowA;

                GrB_Matrix_new(&tmpA, GrB_BOOL, nrowsA, ncolsA);
                GrB_Vector_new(&rowA, GrB_BOOL, ncolsA);

                GrB_Col_extract(rowA, GrB_NULL, GrB_NULL, A, GrB_ALL, ncolsA, row - 1, GrB_DESC_T0);
                GrB_Row_assign(tmpA, GrB_NULL, GrB_NULL, rowA, row - 1, GrB_ALL, ncolsA, GrB_NULL);

                GrB_Matrix_free(&C);
                C = mult(tmpA, B);
                GrB_Matrix_free(&tmpA);
                GrB_Vector_free(&rowA);
            }

            return C;
        }

        static inline GrB_Matrix pow(GrB_Matrix A)
        {
            GrB_Index nrows, ncols;
            GrB_Matrix_nrows(&nrows, A);
            GrB_Matrix_ncols(&ncols, A);

            GrB_Matrix result;
            GrB_Matrix_dup(&result, A);

            GrB_Matrix prev;
            GrB_Matrix_new(&prev, GrB_BOOL, nrows, ncols);

            GrB_Index prev_nvals = 0, result_nvals = 0;
            bool changed;

            do
            {
                GrB_Matrix temp = mult(result, A);
                GrB_Matrix sum_result = sum(result, temp);
                destroy(result);
                destroy(temp);
                result = sum_result;

                GrB_Matrix_nvals(&prev_nvals, prev);
                GrB_Matrix_nvals(&result_nvals, result);
                changed = (result_nvals != prev_nvals);

                destroy(prev);
                GrB_Matrix_dup(&prev, result);
            } while (changed);

            destroy(prev);
            return result;
        }

        // transitive closure of a matrix, pos says if it's + rather than *
        static inline GrB_Matrix clos(GrB_Matrix A, uint pos)
        {
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
            GrB_Matrix S = NULL, E = NULL, P = NULL, T = NULL;

            GrB_Index nrowsA, ncolsA;
            GrB_Matrix_nrows(&nrowsA, A);
            GrB_Matrix_ncols(&ncolsA, A);
            uint dim = mmax(nrowsA, ncolsA);

            if (ID == NULL)
            {
                if (pos)
                    E = empty(dim, dim);
                else
                    E = mat_one(dim, dim, row, row);

                S = sum1(row, A, E, full_side);
                destroy(E);
            }
            else
            {
                if (pos)
                {
                    S = mult1(row, ID, A, full_side);
                }
                else
                {
                    E = empty(dim, dim);
                    S = sum1(row, ID, E, full_side);
                    destroy(E);
                }
            }

            bool accessS = false;
            if (coltest)
                GrB_Matrix_extractElement_BOOL(&accessS, S, row, *coltest);
            if (coltest && accessS)
            {
                destroy(S);
                *coltest = 1;
                return NULL;
            }

            P = mult(S, A);
            T = S;
            S = sum(P, T);
            destroy(T);
            destroy(P);

            GrB_Index elems = 0, elems_new = 0;
            GrB_Matrix_nvals(&elems_new, S);

            while (elems_new != elems)
            {
                elems = elems_new;

                if (coltest)
                    GrB_Matrix_extractElement_BOOL(&accessS, S, row, *coltest);
                if (coltest && accessS)
                {
                    destroy(S);
                    *coltest = 1;
                    return NULL;
                }

                P = mult(S, A);
                T = S;
                S = sum(S, P);
                destroy(T);
                destroy(P);

                GrB_Matrix_nvals(&elems_new, S);
            }

            if (coltest)
            {
                destroy(S);
                *coltest = 0;
                return NULL;
            }

            return S;
        }

        static inline GrB_Matrix clos_col(uint64_t col, GrB_Matrix ID, GrB_Matrix A, uint pos, uint *rowtest)
        {
            GrB_Matrix S = NULL, E = NULL, P = NULL, T = NULL;

            GrB_Index nrowsA, ncolsA;
            GrB_Matrix_nrows(&nrowsA, A);
            GrB_Matrix_ncols(&ncolsA, A);
            uint dim = mmax(nrowsA, ncolsA);

            if (ID == NULL)
            {
                if (pos)
                    E = empty(dim, dim);
                else
                    E = mat_one(dim, dim, col, col);

                S = sum1(full_side, A, E, col);
                destroy(E);
            }
            else
            {
                if (pos)
                {
                    S = mult1(full_side, A,ID, col);
                }
                else
                {
                    E = empty(dim, dim);
                    S = sum1(full_side, ID, E, col);
                    destroy(E);
                }
            }

            bool accessS = false;
            if (rowtest)
                GrB_Matrix_extractElement_BOOL(&accessS, S, *rowtest, col);
            if (rowtest && accessS)
            {
                destroy(S);
                *rowtest = 1;
                return NULL;
            }

            P = mult(A, S);
            T = S;
            S = sum(P, T);
            destroy(T);
            destroy(P);

            GrB_Index elems = 0, elems_new = 0;
            GrB_Matrix_nvals(&elems_new, S);

            while (elems_new != elems)
            {
                elems = elems_new;

                if (rowtest)
                    GrB_Matrix_extractElement_BOOL(&accessS, S, *rowtest, col);
                if (rowtest && accessS)
                {
                    destroy(S);
                    *rowtest = 1;
                    return NULL;
                }

                P = mult(A, S);
                T = S;
                S = sum(S, P);
                destroy(T);
                destroy(P);

                GrB_Matrix_nvals(&elems_new, S);
            }

            if (rowtest)
            {
                destroy(S);
                *rowtest = 0;
                return NULL;
            }

            return S;
        }

        // versions to choose one row or one column, or both
        static inline GrB_Matrix clos1(uint64_t row, GrB_Matrix A, uint pos, uint64_t col)
        {
            GrB_Index nrowsA, ncolsA;
            GrB_Matrix_nrows(&nrowsA, A);
            GrB_Matrix_ncols(&ncolsA, A);
            uint side = mmax(nrowsA, ncolsA);
            uint test = 0;

            if (row == full_side)
            {
                if (col == full_side)
                    return clos(A, pos);
                else
                    return clos_col(col, NULL, A, pos, NULL);
            }
            else if (col == full_side)
            {
                return clos_row(row, NULL, A, pos, NULL);
            }

            GrB_Vector row_vec = NULL, col_vec = NULL;
            GrB_Index nrow = 0, ncol = 0;

            GrB_Vector_new(&row_vec, GrB_BOOL, side);
            GrB_Col_extract(row_vec, NULL, NULL, A, GrB_ALL, side, row - 1, NULL);
            GrB_Vector_nvals(&ncol, row_vec);
            GrB_Vector_free(&row_vec);

            GrB_Vector_new(&col_vec, GrB_BOOL, side);
            GrB_Col_extract(col_vec, NULL, NULL, A, GrB_ALL, side, col - 1, GrB_DESC_T0);
            GrB_Vector_nvals(&nrow, col_vec);
            GrB_Vector_free(&col_vec);

            if (ncol < nrow)
            {
                test = row;
                clos_col(col, NULL, A, pos, &test);
            }
            else
            {
                test = col;
                clos_row(row, NULL, A, pos, &test);
            }

            return test ? mat_one(side, side, row, col) : empty(side, side);
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
                GrB_Matrix M1 = NULL, M = NULL;
                if (col == full_side)
                {
                    M1 = clos(B, pos);
                    M = mult(A, M1);
                }
                else
                {
                    M1 = clos1(full_side, B, pos, col);
                    M = mult(A, M1);
                }
                GrB_Matrix_free(&M1);
                return M;
            }
            else if (col == full_side)
            {
                return clos_row(row, A, B, pos, NULL);
            }

            uint32_t test = col;
            clos_row(row, A, B, pos, &test);
            return test ? mat_one(side, side, row, col) : empty(side, side);
        }

        // computes [row] A* B [col] (pos=0) or [row] A+ B [col] (pos=1)
        static inline GrB_Matrix clos_mult1(uint64_t row, GrB_Matrix A, uint pos, GrB_Matrix B, uint64_t col)
        {
            GrB_Index nrowsA, ncolsA, nrowsB, ncolsB;
            GrB_Matrix_nrows(&nrowsA, A);
            GrB_Matrix_ncols(&ncolsA, A);
            GrB_Matrix_nrows(&nrowsB, B);
            GrB_Matrix_ncols(&ncolsB, B);
            uint64_t side = mmax(nrowsA, ncolsA);

            if (col == full_side)
            {
                GrB_Matrix M1 = NULL, M = NULL;
                if (row == full_side)
                {
                    M1 = clos(B, pos);
                    M = mult(M1, A);
                }
                else
                {
                    M1 = clos1(full_side, B, pos, row);
                    M = mult(M1, A);
                }
                GrB_Matrix_free(&M1);
                return M;
            }
            else if (row == full_side)
            {
                return clos_col(col, B, A, pos, NULL);
            }

            uint32_t test = row;
            clos_col(col, B, A, pos, &test);
            return test ? mat_one(side, side, row, col) : empty(side, side);
        }
    };
}

#endif // RPQ_MATRIX_BASELINEGB_HPP
