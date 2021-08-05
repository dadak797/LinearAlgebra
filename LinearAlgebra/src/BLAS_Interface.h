#pragma once


class BLAS_Interface
{
private:
    static int m_NumThread;

public:
    static void SetNumThreads(int nt);
};