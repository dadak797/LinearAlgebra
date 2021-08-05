#pragma once

#include "ArrayAllocator.h"
#include "VectorCommaAssignment.h"
#include "VectorEvaluation.h"
#include "Range.h"


#ifdef _CHECK_RANGE
  #define VectorAssert(cond) assert(cond)
#else
  #define VectorAssert(cond)
#endif


template <typename AllocatorType = DynamicAllocator<double>>
class Vector : public AllocatorType
{
public:
    using DataType = typename AllocatorType::MemberType;
    using TNumType = typename AllocatorType::MemberType;

private:
    int m_Size;

    template <typename AllocatorType1>
    void CopyVector(const Vector<AllocatorType1>& src)
    {
        const DataType* v1 = src.Data();
        DataType* v0 = AllocatorType::Data();  // "AllocatorType::" should not be omitted and it can be changed to "this->".

        for (int i = 0; i < m_Size; i++)
        {
            v0[i] = v1[i];
        }
    }

public:
    using Iterator = DataType*;
    using ConstIterator = const DataType*;
    //    typedef DataType* Iterator;
    //    typedef const DataType* ConstIterator;

    Iterator Begin() { return AllocatorType::Data(); }
    Iterator End() { return AllocatorType::Data() + Size(); }
    ConstIterator Begin() const { return AllocatorType::Data(); }
    ConstIterator End() const { return AllocatorType::Data() + Size(); }

    // Constructors
    Vector()
        : m_Size(AllocatorType::DEFAULT_SIZE), AllocatorType() {}

    Vector(int size)
        : m_Size(size), AllocatorType(size) {}

    Vector(const DataType* pDataType, int size)
        : m_Size(size), AllocatorType(const_cast<DataType*>(pDataType)) {}

    Vector(const Vector<AllocatorType>& other)
        : m_Size(other.m_Size), AllocatorType(other.m_Size)
    {
        CopyVector(other);
    }

    ~Vector() { AllocatorType::DeAllocate(); }

    void Resize(int size)
    {
        if (m_Size == size) return;
        AllocatorType::DeAllocate();
        m_Size = size;
        AllocatorType::Allocate(m_Size);
    }

    void Clear() { Resize(0); }
    int Size() const { return m_Size; }
    int Stride() const { return 1; }

    void Swap(Vector& other)
    { 
        std::swap(m_Size, other.m_Size);
        AllocatorType::Swap(other);
    }

    void Reserve(int size) 
    { 
        if (m_Size < size)
        {
            Resize(size);
        }
    }

    const DataType& operator[](int i) const
    { 
        VectorAssert(i < m_Size);
        return *AllocatorType::Data(i);
    }

    DataType& operator[](int i)
    {
        VectorAssert(i < m_Size);
        return *AllocatorType::Data(i);
    }

    const DataType& operator()(int i) const
    {
        VectorAssert(i < m_Size);
        return *AllocatorType::Data(i);
    }

    DataType& operator()(int i)
    {
        VectorAssert(i < m_Size);
        return *AllocatorType::Data(i);
    }
};


//
//    //Vector<SharedAllocator<DataType>> operator()(const Range& range) const
//    //{
//    //    return Vector<SharedAllocator<DataType>>(const_cast<DataType*>(Data(range.First())), range.Count(m_Size - 1));
//    //}
//
//    //Vector<SharedAllocator<DataType>> operator()(const Range& range)
//    //{
//    //    return Vector<SharedAllocator<DataType>>(const_cast<DataType*>(Data(range.First())), range.Count(m_Size - 1));
//    //}
//
//    //VectorCommaAssignment<AllocatorType> operator=(const DataType& value)
//    //{
//    //    if (m_Size > 0)
//    //    {
//    //        *Data() = value;
//    //    }
//
//    //    return VectorCommaAssignment<AllocatorType>(this, 1, value);
//    //}
//
//    //Vector& operator=(const Vector& other)
//    //{
//    //    Resize(other.Size());
//    //    CopyVector(other);
//    //    return *this;
//    //}
//
//    //template <typename AllocaterType1>
//    //Vector& operator=(const Vector<AllocaterType1>& other);
//    ////{
//    ////    Resize(other.Size());
//    ////    CopyVector(other);
//    ////    return *this;
//    ////}
//
//    //template <typename VectorType1, typename VectorType2, typename OprType>
//    //Vector& operator=(const Vec2Expr<VectorType1, VectorType2, OprType>& vov);
//    ////{
//    ////    DataType* v3 = Data();
//    ////    for (int i = 0; i < m_Size; i++)
//    ////    {
//    ////        v3[i] = vov[i];
//    ////    }
//
//    ////    return *this;
//    ////}
//
//    //Vector& operator+=(const DataType& v1);
//    ////{
//    ////    DataType* v2 = Data();
//    ////    for (int i = 0; i < m_Size; i++)
//    ////    {
//    ////        v2[i] += v1;
//    ////    }
//    ////    return *this;
//    ////}
//    //
//    //template <typename AllocatorType> 
//    //Vector& operator+=(const Vector<AllocatorType>& other)
//    //{
//    //    assert(m_Size == other.m_Size);
//    //    const DataType* v1 = other.Data();
//    //    DataType* v2 = Data();
//    //    for (int i = 0; i < m_Size; i++)
//    //    {
//    //        v2[i] += v1[i];
//    //    }
//
//    //    return *this;
//    //}
//
//    //template <typename VectorType1, typename VectorType2, typename OprType>
//    //Vector& operator+=(const Vec2Expr<VectorType1, VectorType2, OprType>& vov)
//    //{
//    //    DataType* v = Data();
//    //    for (int i = 0; i < m_Size; i++)
//    //    {
//    //        v[i] += vov[i];
//    //    }
//
//    //    return *this;
//    //}
//
//    //Vector& operator-=(const DataType& v1)
//    //{
//    //    DataType* v2 = Data();
//    //    for (int i = 0; i < m_Size; i++)
//    //    {
//    //        v2[i] -= v1;
//    //    }
//    //    return *this;
//    //}
//
//    //template <typename AllocatorType>
//    //Vector& operator-=(const Vector<AllocatorType>& other)
//    //{
//    //    assert(m_Size == other.m_Size);
//    //    const DataType* v1 = other.Data();
//    //    DataType* v2 = Data();
//    //    for (int i = 0; i < m_Size; i++)
//    //    {
//    //        v2[i] -= v1[i];
//    //    }
//
//    //    return *this;
//    //}
//
//    //template <typename VectorType1, typename VectorType2, typename OprType>
//    //Vector& operator-=(const Vec2Expr<VectorType1, VectorType2, OprType>& vov)
//    //{
//    //    DataType* v = Data();
//    //    for (int i = 0; i < m_Size; i++)
//    //    {
//    //        v[i] -= vov[i];
//    //    }
//
//    //    return *this;
//    //}
//
//    //Vector& operator*=(const DataType& v1)
//    //{
//    //    DataType* v2 = Data();
//    //    for (int i = 0; i < m_Size; i++)
//    //    {
//    //        v2[i] *= v1;
//    //    }
//    //    return *this;
//    //}
//
//    //template <typename AllocatorType>
//    //Vector& operator*=(const Vector<AllocatorType>& other)
//    //{
//    //    assert(m_Size == other.m_Size);
//    //    const DataType* v1 = other.Data();
//    //    DataType* v2 = Data();
//    //    for (int i = 0; i < m_Size; i++)
//    //    {
//    //        v2[i] *= v1[i];
//    //    }
//
//    //    return *this;
//    //}
//
//    //template <typename VectorType1, typename VectorType2, typename OprType>
//    //Vector& operator*=(const Vec2Expr<VectorType1, VectorType2, OprType>& vov)
//    //{
//    //    DataType* v = Data();
//    //    for (int i = 0; i < m_Size; i++)
//    //    {
//    //        v[i] *= vov[i];
//    //    }
//
//    //    return *this;
//    //}
//
//    //Vector& operator/=(const DataType& v1)
//    //{
//    //    DataType* v2 = Data();
//    //    for (int i = 0; i < m_Size; i++)
//    //    {
//    //        v2[i] /= v1;
//    //    }
//    //    return *this;
//    //}
//
//    //template <typename AllocatorType>
//    //Vector& operator/=(const Vector<AllocatorType>& other)
//    //{
//    //    assert(m_Size == other.m_Size);
//    //    const DataType* v1 = other.Data();
//    //    DataType* v2 = Data();
//    //    for (int i = 0; i < m_Size; i++)
//    //    {
//    //        v2[i] /= v1[i];
//    //    }
//
//    //    return *this;
//    //}
//
//    //template <typename VectorType1, typename VectorType2, typename OprType>
//    //Vector& operator/=(const Vec2Expr<VectorType1, VectorType2, OprType>& vov)
//    //{
//    //    DataType* v = Data();
//    //    for (int i = 0; i < m_Size; i++)
//    //    {
//    //        v[i] /= vov[i];
//    //    }
//
//    //    return *this;
//    //}
//};


//template <typename AllocatorType>
//template <typename AllocatorType1>
//void Vector<AllocatorType>::CopyVector(const Vector<AllocatorType1>& src)
//{
//    const DataType* v1 = src.Data();
//    DataType* v0 = Data();
//
//    for (int i = 0; i < m_Size; i++)
//    {
//        v0[i] = v1[i];
//    }
//}


//template <typename AllocatorType>
//template <typename AllocaterType1>
//Vector<AllocatorType>& Vector<AllocatorType>::operator=(const Vector<AllocaterType1>& other)
//{
//    Resize(other.Size());
//    CopyVector(other);
//    return *this;
//}
//
//
//template <typename AllocatorType>
//template <typename VectorType1, typename VectorType2, typename OprType>
//Vector<AllocatorType>& Vector<AllocatorType>::operator=(const Vec2Expr<VectorType1, VectorType2, OprType>& vov)
//{
//    DataType* v3 = Data();
//    for (int i = 0; i < m_Size; i++)
//    {
//        v3[i] = vov[i];
//    }
//
//    return *this;
//}
//
//
//template <typename AllocatorType>
//Vector<AllocatorType>& Vector<AllocatorType>::operator+=(const DataType& v1)
//{
//    DataType* v2 = Data();
//    for (int i = 0; i < m_Size; i++)
//    {
//        v2[i] += v1;
//    }
//    return *this;
//}