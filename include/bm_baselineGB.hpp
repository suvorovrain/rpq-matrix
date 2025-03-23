//
// Created by Adrián on 16/4/24.
//

#ifndef RPQ_MATRIX_BASELINEGB_HPP
#define RPQ_MATRIX_BASELINEGB_HPP
#define fullSide (~(uint)0) // (*** got this from original matrix.h include. Maybe remove later idk ***)
namespace bm_baselinegb
{

    extern "C"
    {
#include "../../../GraphBLAS/Include/GraphBLAS.h" // here put GB matirx
#include "../../../LAGraph/include/LAGraph.h"
#include "baseline/utilstime.h" // here all good
    }

    class wrapper
    {
        // TODO: add error handling. and think about memory leaks
    public:
        static const uint64_t full_side = fullSide;
        typedef GrB_Matrix matrix_type;    // (*** in original this is a pointer to matrix ***)
        typedef GrB_Matrix *s_matrix_type; // (*** and this is matrix itself ***)

        /*    (*** this used only in baseline_build.cpp. I'll just changed it to GB matrix creation in that file.***)

         // creates matrix of width x height with n cells (2n ints row,col)
         // reorders cells array
         static inline matrix_type create(uint64_t height, uint64_t width, uint64_t n, uint *cells){
             return matCreate(height, width, n, cells);
         };
        */

        static void time_begin()
        {
            time_beg();
        }

        // destroys matrix (*** it's just a free for matrix. changed it to GB free***)
        static inline void destroy(GrB_Matrix M)
        {
            GrB_Matrix_free(&M);
        };

        /* (*** this function is using only in lib matrix logic and utils. GB by default creates empty matrix (i hope so).***) */
        // creates an empty matrix
        static inline GrB_Matrix empty(uint64_t height, uint64_t width)
        {
            GrB_Matrix A;
            GrB_Matrix_new(&A, GrB_BOOL, height, width);
            return A; // not sure about matrix type, but i think it's bool
        };

        // creates an identity matrix (*** i moved here logic of creation of id matrix via GB tools ***)
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

        // (*** this is using in rpq_solver. i just made matrix with one non zero element***)
        // version of identity matrix with one row or column
        static inline matrix_type id1(uint64_t side, uint64_t rc)
        {
            GrB_Matrix E;
            GrB_Matrix_new(&E, GrB_BOOL, side, side);
            GrB_Matrix_setElement_BOOL(E, true, rc, rc);
            return E;
        };

        /* (*** seems like this is not used in rpq solver logic***)
        // creates a new copy of A, with its own data
        static inline matrix copy(matrix A){
            return matCopy(A);
        };
        */

        // (*** ATTENTION: danger place here. in original this function return struct, not a pointer.***)
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

        //  (*** ATTENTION: this is simple funciton for result saving. I need to use GB method for this. ***)
        // writes M to file, which must be opened for writing
        static inline void save(GrB_Matrix M, FILE *file)
        {
            GxB_Matrix_fprint(M, "M", GxB_COMPLETE, file);
        };

        // (*** ATTENTION: used in rpq_solver. this is maybe bad. cuz LAGraph_MMRead read mtx file ***)
        // loads matrix from file, which must be opened for reading
        static inline GrB_Matrix load(FILE *file)
        {
            GrB_Matrix A;
            LAGraph_MMRead(&A, file, "message");
            return A;
        };

        /* (*** ATTENTION: this method returns the number of bytes that matrix occupies.
        maybe GB method works other way ***) */
        // space of the matrix, in w-bit words
        static inline uint64_t space(GrB_Matrix M)
        {
            uint64_t memory_usage;
            return GxB_Matrix_memoryUsage(&memory_usage, M);
            return memory_usage;
        };

        /* (*** i guess i don't need it***)
        // dimensions of M, returns elems
        static inline uint64_t dims(matrix M, uint *logside, uint *width, uint *height){
            return matDims(M, width, height);
        };
        */
        /* (*** i guess i don't need it***)
            // accesses a cell
            static inline uint access(matrix M, uint64_t row, uint64_t col){
                return matAccess(M, row, col);
            }
        */

        // (*** simple A + B from authors paper. Did it with GB ***)
        // (boolean) sum of two matrices, assumed to be of the same side
        static inline GrB_Matrix sum(GrB_Matrix A, GrB_Matrix B)
        {
            GrB_Index nrows, ncols;
            GrB_Matrix C;
            GrB_Matrix_nrows(&nrows, A);
            GrB_Matrix_ncols(&ncols, A);
            GrB_Matrix_new(&C, GrB_BOOL, nrows, ncols);
            GrB_Matrix_eWiseAdd_BinaryOp(C, GrB_NULL, GrB_NULL, GrB_LOR, A, B, GrB_NULL);
            return C;
            // return matSum(A, B);
        };

        //(*** this method is used in cases when we apply restriction optimizations. I hope i made it right ***)
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
                // if we don't have both restrictions, just call common matrix addition
                return sum(A, B);
            }
            else if (row != full_side && col != full_side)
            {
                // if we got both restrictions, just do logical OR on two elements from A and B
                bool elA, elB;
                GrB_Matrix_extractElement_BOOL(&elA, A, row, col);
                GrB_Matrix_extractElement_BOOL(&elB, B, row, col);
                GrB_Matrix_setElement_BOOL(C, elA || elB, row, col);
                return C;
            }
            else if (row == full_side && col != full_side)
            {
                // get corresponding columns
                GrB_Vector colA, colB, colC;
                GrB_Col_extract(colA, GrB_NULL, GrB_NULL, A, GrB_ALL, row, col, GrB_NULL); // is "row" on the right place?
                GrB_Col_extract(colB, GrB_NULL, GrB_NULL, B, GrB_ALL, row, col, GrB_NULL); // is "row" on the right place?

                // A + B columns
                GrB_Vector_eWiseAdd_BinaryOp(colC, GrB_NULL, GrB_NULL, GrB_LOR, colA, colB, GrB_NULL);

                // save result column in C matrix
                GrB_Col_assign(C, GrB_NULL, GrB_NULL, colC, GrB_ALL, row, col, GrB_NULL); // is "row" on the right place?
                return C;
            }
            else
            {
                GrB_Matrix AT, BT, CT;
                AT = transpose(A);
                BT = transpose(B);
                C = sum1(col, A, B, row);
                CT = transpose(C);
                return CT;
            }
        };

        // (boolean) product of two matrices, assumed to be of the same side
        // only rowA of A and colB of B are considered if not fullSide
        static inline GrB_Matrix mult(GrB_Matrix A, GrB_Matrix B)
        {
            GrB_Index nrows, ncols;
            GrB_Matrix C;
            GrB_Matrix_nrows(&nrows, A);
            GrB_Matrix_ncols(&ncols, B);
            GrB_Matrix_new(&C, GrB_BOOL, nrows, ncols);
            GrB_mxm(C, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, A, B, GrB_NULL); // i hope i chose right semiring :)
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
                // if we don't have both restrictions, just call common matrix multiplication
                return mult(A, B);
            }
            else if (row != full_side && col != full_side)
            {
                // <r>(A x B)<c> = (<r>A) x (B<c>)
                GrB_Matrix tmpA, tmpB;
                GrB_Matrix_new(&tmpA, GrB_BOOL, nrowsA, ncolsA);
                GrB_Matrix_new(&tmpB, GrB_BOOL, nrowsB, ncolsB);

                // (B<c>)
                GrB_Vector colB;
                GrB_Col_extract(colB, GrB_NULL, GrB_NULL, B, GrB_ALL, row, col, GrB_NULL);
                GrB_Col_assign(tmpB, GrB_NULL, GrB_NULL, colB, GrB_ALL, row, col, GrB_NULL);

                //(<r>A) need to double check correct usage of operations below
                GrB_Vector rowA;
                GrB_Matrix AT;
                AT = transpose(A);
                GrB_Col_extract(rowA, GrB_NULL, GrB_NULL, AT, GrB_ALL, col, row, GrB_NULL);
                GrB_Row_assign(tmpA, GrB_NULL, GrB_NULL, rowA, row, GrB_ALL, col, GrB_NULL);

                C = mult(tmpA, tmpB);
                return C;
            }
            else if (row == full_side && col != full_side)
            {
                // (A x B)<c> = A x (B<c>)
                GrB_Matrix tmpB;
                GrB_Matrix_new(&tmpB, GrB_BOOL, nrowsB, ncolsB);

                // (B<c>)
                GrB_Vector colB;
                GrB_Col_extract(colB, GrB_NULL, GrB_NULL, B, GrB_ALL, row, col, GrB_NULL);
                GrB_Col_assign(tmpB, GrB_NULL, GrB_NULL, colB, GrB_ALL, row, col, GrB_NULL);
                C = mult(A, tmpB);
                return C;
            }
            else
            {
                // <r>(A x B) = (<r>A) x B
                GrB_Matrix tmpA;
                GrB_Matrix_new(&tmpA, GrB_BOOL, nrowsA, ncolsA);

                //(<r>A) need to double check correct usage of operations below
                GrB_Vector rowA;
                GrB_Matrix AT;
                AT = transpose(A);
                GrB_Col_extract(rowA, GrB_NULL, GrB_NULL, AT, GrB_ALL, col, row, GrB_NULL);
                GrB_Row_assign(tmpA, GrB_NULL, GrB_NULL, rowA, row, GrB_ALL, col, GrB_NULL);

                C = mult(tmpA, B);
                return C;
            }
        };

        static inline GrB_Matrix pow(GrB_Matrix A, uint64_t row, uint64_t col)
        {
            GrB_Matrix plusA, powA; // plusA stores sum of matrix powers
            GrB_Index nrows, ncols;
            GrB_Matrix_nrows(&nrows, A);
            GrB_Matrix_ncols(&ncols, A);
            GrB_Matrix_new(&plusA, GrB_BOOL, nrows, ncols);
            GrB_Matrix_new(&powA, GrB_BOOL, nrows, ncols);
            GrB_Matrix_dup(&powA, A);
            GrB_Matrix_dup(&plusA, A);
            bool is_equal = false;
            while (!is_equal)
            {
                // A^{k} = A^{k-1} * A
                powA = mult1(row, powA, A, col);
                // GrB_mxm(powA, GrB_NULL, GrB_NULL, GrB_LOR_LAND_SEMIRING_BOOL, powA, A, GrB_NULL);

                // compute (A + ... + A^{k-1}) + A^k
                GrB_Matrix summ;
                summ = sum1(row, plusA, powA, col);

                // compare (A + ... + A^{k-1}) and (A + ... + A^{k-1} + A^k)
                GrB_Matrix compare;
                GrB_Matrix_new(&compare, GrB_BOOL, nrows, ncols);
                GrB_Matrix_eWiseAdd_BinaryOp(compare, GrB_NULL, GrB_NULL, GrB_EQ_BOOL, summ, plusA, GrB_NULL);
                GrB_Matrix_reduce_BOOL(&is_equal, GrB_NULL, GrB_LAND_MONOID_BOOL, compare, GrB_NULL);

                plusA = summ;
            }
            return plusA;
        }

        // transitive closure of a matrix, pos says if it's + rather than *
        static inline GrB_Matrix clos(GrB_Matrix A, uint pos)
        {
            GrB_Matrix C;
            C = pow(A, full_side, full_side);
            if (!pos)
            { // (*** Kleene Star ***)
                GrB_Matrix E;
                GrB_Index nrows, ncols;
                GrB_Matrix_nrows(&nrows, A);
                GrB_Matrix_ncols(&ncols, A);
                E = id(nrows);
                C = sum(E, C); // can i pass C as argument hmmmmmm...
            }
            return C;
        };

        // versions to choose one row or one column, or both
        static inline GrB_Matrix clos1(uint64_t row, GrB_Matrix A, uint pos, uint64_t col)
        {
            GrB_Matrix C;
            C = pow(A, row, col);
            if (!pos)
            { // (*** Kleene Star ***)
                GrB_Matrix E;
                GrB_Index nrows, ncols;
                GrB_Matrix_nrows(&nrows, A);
                GrB_Matrix_ncols(&ncols, A);
                E = id(nrows);
                C = sum1(row,E, C,col); // can i pass C as argument hmmmmmm...
            }
            return C;
        };

        // computes [row] A B* [col] (pos=0) or [row] A B+ [col] (pos=1)
        static inline GrB_Matrix mult_clos1 (uint64_t row, GrB_Matrix A, GrB_Matrix B, uint pos, uint64_t col){
            
            if (row == full_side && col == full_side){
                GrB_Matrix M1,M;
                M1 = clos(B,pos);
                M = mult(A,M1);
                GrB_Matrix_free(&M1);
                return M;
            } else if (row != full_side && col != full_side){
                
            }
            // return matMultClos1(row, A, B, pos, col);
        };

        //     // computes [row] A* B [col] (pos=0) or [row] A+ B [col] (pos=1)
        //     static inline matrix clos_mult1 (uint64_t row, matrix A, uint pos, matrix B, uint64_t col){
        //         return matClosMult1(row, A, pos, B, col);
        //     };
    };
}

#endif // RPQ_MATRIX_BASELINEGB_HPP
