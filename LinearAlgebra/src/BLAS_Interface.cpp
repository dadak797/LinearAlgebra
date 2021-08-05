#pragma warning (disable:4819)

#include "pch.h"
#include "BLAS_Interface.h"

#include "mkl.h"
#include <omp.h>


inline void BLAS_Interface::SetNumThreads(int nt)
{
    m_NumThread = nt;
    MKL_Set_Num_Threads(nt);
    omp_set_num_threads(nt);
}