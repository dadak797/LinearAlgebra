#pragma once
#pragma warning (disable:4819)

#include "mkl.h"
#include <omp.h>


class BLAS_Interface
{
private:
    static int s_NumThread;

public:
    static void SetNumThreads(int nt);
    static int GetNumThreads();

    // Dot Product with BLAS
    //template <typename DataType>
    //static DataType DOT(int n, const DataType* x, int incx, const DataType* y, int incy);
};


inline void BLAS_Interface::SetNumThreads(int nt)
{
    s_NumThread = nt;
    MKL_Set_Num_Threads(nt);
    omp_set_num_threads(nt);
}

inline int BLAS_Interface::GetNumThreads()
{
    return s_NumThread;
}