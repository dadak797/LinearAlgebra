#pragma once


template<typename AllocatorType>
class Vector;

template<typename MemberType>
class Scal2Vec 
{
public:
    //typedef MemberType DataType;
    using DataType = MemberType;

    const DataType& m_Scalar;
    Scal2Vec(const DataType& scalar)
        : m_Scalar(scalar) {}
    DataType operator[](int i) const
    {
        return m_Scalar;
    }
};

template<typename VectorType1, typename VectorType2, typename OprType>
class Vec2Expr
{
private:
    const VectorType1& m_Vec1;
    const VectorType2& m_Vec2;

public:
    using DataType = typename VectorType1::DataType;
    //typedef typename VectorType1::DataType DataType;
    
    Vec2Expr(const VectorType1& vec1, const VectorType2& vec2)
        : m_Vec1(vec1), m_Vec2(vec2) {}
    DataType operator[](int i) const
    {
        return OprType::Compute(m_Vec1[i], m_Vec2[i]);
    }
    int Size() const { return m_Vec1.Size(); }
};

// OprType VecAdd
// opr +

struct VecAdd 
{
    template<typename DataType> 
    static DataType Compute(const DataType& v1, const DataType& v2) 
    { 
        return v1 + v2; 
    }
};

template<typename VectorTypeR1, typename VectorTypeL1, typename OprType1, 
    typename VectorTypeR2, typename VectorTypeL2, typename OprType2>
inline Vec2Expr<Vec2Expr<VectorTypeR1, VectorTypeL1, OprType1>, Vec2Expr<VectorTypeR2, VectorTypeL2, OprType2>, VecAdd> operator+(
    const Vec2Expr<VectorTypeR1, VectorTypeL1, OprType1>& v1, const Vec2Expr<VectorTypeR2, VectorTypeL2, OprType2>& v2)
{
    return Vec2Expr<Vec2Expr<VectorTypeR1, VectorTypeL1, OprType1>, Vec2Expr<VectorTypeR2, VectorTypeL2, OprType2>, VecAdd>(v1, v2);
}

template<typename VectorTypeR1, typename VectorTypeL1, typename OprType1, typename AllocatorType>
inline Vec2Expr<Vec2Expr<VectorTypeR1, VectorTypeL1, OprType1>, Vector<AllocatorType>, VecAdd> operator+(
    const Vec2Expr<VectorTypeR1, VectorTypeL1, OprType1>& v1, const Vector<AllocatorType>& v2)
{
    return Vec2Expr<Vec2Expr<VectorTypeR1, VectorTypeL1, OprType1>, Vector<AllocatorType>, VecAdd>(v1, v2);
}

template<typename VectorTypeR1, typename VectorTypeL1, typename OprType1, typename AllocatorType>
inline Vec2Expr<Vector<AllocatorType>, Vec2Expr<VectorTypeR1, VectorTypeL1, OprType1>, VecAdd> operator+(
    const Vector<AllocatorType>& v1, const Vec2Expr<VectorTypeR1, VectorTypeL1, OprType1>& v2)
{
    return Vec2Expr<Vector<AllocatorType>, Vec2Expr<VectorTypeR1, VectorTypeL1, OprType1>, VecAdd>(v1, v2); 
}

template<typename AllocatorType1, typename AllocatorType2>
inline Vec2Expr<Vector<AllocatorType1>, Vector<AllocatorType2>, VecAdd> operator+(
    const Vector<AllocatorType1>& v1, const Vector<AllocatorType2>& v2)
{
    return Vec2Expr<Vector<AllocatorType1>, Vector<AllocatorType2>, VecAdd>(v1, v2); 
}

template<typename AllocatorType>
inline Vec2Expr<Scal2Vec<typename AllocatorType::MemberType>, Vector<AllocatorType>, VecAdd> operator+(
    const Scal2Vec<typename AllocatorType::MemberType>& v1, const Vector<AllocatorType>& v2)
{
    return Vec2Expr<Scal2Vec<typename AllocatorType::member_type>, Vector<AllocatorType>, VecAdd>(v1, v2); 
}

template<typename VectorTypeR1, typename VectorTypeL1, typename OprType1>
inline Vec2Expr<Scal2Vec<typename VectorTypeR1::DataType>, Vec2Expr<VectorTypeR1, VectorTypeL1, OprType1>, VecAdd> operator+(
    const Scal2Vec<typename VectorTypeR1::DataType>& v1, const Vec2Expr<VectorTypeR1, VectorTypeL1, OprType1>& v2)
{
    return Vec2Expr<Scal2Vec<typename VectorTypeR1::data_type>, Vec2Expr<VectorTypeR1, VectorTypeL1, OprType1>, VecAdd>(v1, v2); 
}