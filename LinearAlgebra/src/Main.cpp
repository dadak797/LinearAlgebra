#include "pch.h"

#include "Vector.h"
#include <iostream>
#include "Matrix.h"


int main()
{
    Matrix matrix(2, 4);

    matrix(0, 0) = 1;
    matrix(0, 1) = 2;
    matrix(0, 2) = 3;
    matrix(0, 3) = 4;

    matrix(1, 0) = 5;
    matrix(1, 1) = 6;
    matrix(1, 2) = 7;
    matrix(1, 3) = 8;

    //matrix(2, 0) = 9;
    //matrix(2, 1) = 10;
    //matrix(2, 2) = 11;
    //matrix(2, 3) = 12;

    for (int i = 0; i < matrix.RowSize(); i++)
    {
        for (int j = 0; j < matrix.ColSize(); j++)
        {
            std::cout << matrix(i, j) << " ";
        }
        std::cout << std::endl;
    }

    Matrix matB(matrix);

    for (int i = 0; i < matB.RowSize(); i++)
    {
        for (int j = 0; j < matB.ColSize(); j++)
        {
            std::cout << matB(i, j) << " ";
        }
        std::cout << std::endl;
    }

    for (int i = 0; i < matrix.RowSize(); i++)
    {
        for (int j = 0; j < matrix.ColSize(); j++)
        {
            std::cout << matrix(i, j) << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}