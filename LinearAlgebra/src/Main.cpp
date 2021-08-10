#include "pch.h"

#include "Vector.h"
#include <iostream>
#include "BLASInterface.h"


int main()
{
    Vector<DynamicAllocator<double>> matA(9);

    matA = 3.0, 1.0, 1.0
        , 1.0, 2.0, 2.0
        , 1.0, 2.0, 2.0;

    Vector<DynamicAllocator<double>> vecW(3);

    vecW = 0.0, 0.0, 0.0;

    int res = BLASInterface::SYEV(LAPACK_ROW_MAJOR, 'V', 'L', 3, matA.Data(), 3, vecW.Data());

    std::cout << "Return: " << res << std::endl;

   
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            std::cout << matA[i * 3 + j] << " ";
        }
        std::cout << std::endl;
    }

    //for (const auto& value : matA)
      //  std::cout << value << std::endl;

    for (const auto& value : vecW)
        std::cout << value << std::endl;

    return 0;
}