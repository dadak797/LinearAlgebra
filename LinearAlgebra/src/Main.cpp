#include "pch.h"

#include "Vector.h"
#include <iostream>
#include "BLASInterface.h"


int main()
{
    Vector<DynamicAllocator<double>> matA(9);

    matA = 1.0, -1.0, 0.0
        , -1.0, 5.0, -4.0;

    int* ipiv = new int[10];

    int res = BLASInterface::GETRF(LAPACK_ROW_MAJOR, 2, 3, matA.Data(), 3, ipiv);

    std::cout << "return value: " << res << std::endl;
    for (int i = 0; i < 10; i++)
        std::cout << "ipiv: " << ipiv[i] << std::endl;
    
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            std::cout << matA[i * 3 + j] << " ";
        }
        std::cout << std::endl;
    }

    //Vector<DynamicAllocator<double>> matB(3);
    //matB = -2.0, -6.0, 38.0;

    //res = BLASInterface::SYTRS(LAPACK_ROW_MAJOR, 'L', 3, 1, matA.Data(), 3, ipiv, matB.Data(), 1);

    //for (const auto& value : matB)
    //    std::cout << value << std::endl;

    delete ipiv;

    return 0;
}