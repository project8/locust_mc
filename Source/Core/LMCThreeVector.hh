#ifndef LMCTHREEVECTOR_H_
#define LMCTHREEVECTOR_H_

#include <cmath>
#include <iostream>

namespace locust
{

    class LMCThreeVector
    {
        public:
            static const LMCThreeVector sInvalid;
            static const LMCThreeVector sZero;

            static const LMCThreeVector sXUnit;
            static const LMCThreeVector sYUnit;
            static const LMCThreeVector sZUnit;

        public:
            LMCThreeVector();
            //~LMCThreeVector();

            //assignment

    LMCThreeVector( const LMCThreeVector& aVector );
    LMCThreeVector& operator=( const LMCThreeVector& aVector );

    LMCThreeVector( const double anArray[ 3 ] );
    LMCThreeVector& operator=( const double anArray[ 3 ] );

    LMCThreeVector( const double& aX, const double& aY, const double& aZ );
    void SetComponents( const double& aX, const double& aY, const double& aZ );
    void SetComponents( const double* aData );
    void SetMagnitude( const double& aMagnitude );
    void SetX( const double& aX );
    void SetY( const double& aY );
    void SetZ( const double& aZ );
    void SetPolarAngle( const double& anAngle );
    void SetAzimuthalAngle( const double& anAngle );

    //cast

    operator double*();
    operator const double*() const;

    //access

    double& operator[]( int anIndex );
    const double& operator[]( int anIndex ) const;

    double& X();
    const double& X() const;
    const double& GetX() const;
    double& Y();
    const double& Y() const;
    const double& GetY() const;
    double& Z();
    const double& Z() const;
    const double& GetZ() const;

    const double* Components() const;

    //comparison

    bool operator==( const LMCThreeVector& aVector ) const;
    bool operator!=( const LMCThreeVector& aVector ) const;
    bool operator<( const LMCThreeVector& aVector ) const;

    //properties

    double Dot( const LMCThreeVector& aVector ) const;
    double Magnitude() const;
    double MagnitudeSquared() const;
    double Perp() const;
    double PerpSquared() const;
    double PolarAngle() const;
    double AzimuthalAngle() const;
    LMCThreeVector Unit() const;
    LMCThreeVector Orthogonal() const;
    LMCThreeVector Cross( const LMCThreeVector& aVector ) const;

private:
    double fData[ 3 ];
};

inline LMCThreeVector::LMCThreeVector( const LMCThreeVector& aVector )
{
    fData[ 0 ] = aVector.fData[ 0 ];
    fData[ 1 ] = aVector.fData[ 1 ];
    fData[ 2 ] = aVector.fData[ 2 ];
}
inline LMCThreeVector& LMCThreeVector::operator=( const LMCThreeVector& aVector )
{
    fData[ 0 ] = aVector.fData[ 0 ];
    fData[ 1 ] = aVector.fData[ 1 ];
    fData[ 2 ] = aVector.fData[ 2 ];
    return *this;
}

inline LMCThreeVector::LMCThreeVector( const double anArray[ 3 ] )
{
    fData[ 0 ] = anArray[ 0 ];
    fData[ 1 ] = anArray[ 1 ];
    fData[ 2 ] = anArray[ 2 ];
}
inline LMCThreeVector& LMCThreeVector::operator=( const double anArray[ 3 ] )
{
    fData[ 0 ] = anArray[ 0 ];
    fData[ 1 ] = anArray[ 1 ];
    fData[ 2 ] = anArray[ 2 ];
    return *this;
}

inline LMCThreeVector::LMCThreeVector( const double& aX, const double& aY, const double& aZ )
{
    fData[ 0 ] = aX;
    fData[ 1 ] = aY;
    fData[ 2 ] = aZ;
}
inline void LMCThreeVector::SetComponents( const double& aX, const double& aY, const double& aZ )
{
    fData[ 0 ] = aX;
    fData[ 1 ] = aY;
    fData[ 2 ] = aZ;
}
inline void LMCThreeVector::SetComponents( const double* aData )
{
    fData[ 0 ] = aData[ 0 ];
    fData[ 1 ] = aData[ 1 ];
    fData[ 2 ] = aData[ 2 ];
}
inline void LMCThreeVector::SetMagnitude( const double& aMagnitude )
{
    const double tMagnitude = Magnitude();
    const double tRatio = aMagnitude / tMagnitude;
    fData[ 0 ] *= tRatio;
    fData[ 1 ] *= tRatio;
    fData[ 2 ] *= tRatio;
    return;
}
inline void LMCThreeVector::SetX( const double& aX )
{
    fData[ 0 ] = aX;
}
inline void LMCThreeVector::SetY( const double& aY )
{
    fData[ 1 ] = aY;
}
inline void LMCThreeVector::SetZ( const double& aZ )
{
    fData[ 2 ] = aZ;
}
inline void LMCThreeVector::SetAzimuthalAngle( const double &anAngle )
{
    const double tRadius = Perp();
    SetComponents( tRadius * cos( anAngle ), tRadius * sin( anAngle ), Z() );
}
inline void LMCThreeVector::SetPolarAngle( const double &anAngle )
{
    const double tMagnitude = Magnitude();
    const double tRadius = Perp();
    SetComponents( tMagnitude * X() / tRadius * sin( anAngle ), tMagnitude * Y() / tRadius * sin( anAngle ), tMagnitude * cos( anAngle ) );
}
inline LMCThreeVector::operator double *()
{
    return fData;
}
inline LMCThreeVector::operator const double *() const
{
    return fData;
}

inline double& LMCThreeVector::operator[]( int anIndex )
{
    return fData[ anIndex ];
}
inline const double& LMCThreeVector::operator[]( int anIndex ) const
{
    return fData[ anIndex ];
}

inline double& LMCThreeVector::X()
{
    return fData[ 0 ];
}
inline const double& LMCThreeVector::X() const
{
    return fData[ 0 ];
}
inline const double& LMCThreeVector::GetX() const
{
    return fData[ 0 ];
}
inline double& LMCThreeVector::Y()
{
    return fData[ 1 ];
}
inline const double& LMCThreeVector::Y() const
{
    return fData[ 1 ];
}
inline const double& LMCThreeVector::GetY() const
{
    return fData[ 1 ];
}
inline double& LMCThreeVector::Z()
{
    return fData[ 2 ];
}
inline const double& LMCThreeVector::Z() const
{
    return fData[ 2 ];
}
inline const double& LMCThreeVector::GetZ() const
{
    return fData[ 2 ];
}
inline const double* LMCThreeVector::Components() const
{
    return (const double*)fData;
}

inline double LMCThreeVector::Dot( const LMCThreeVector& aVector ) const
{
    return (fData[ 0 ] * aVector.fData[ 0 ] + fData[ 1 ] * aVector.fData[ 1 ] + fData[ 2 ] * aVector.fData[ 2 ]);
}

inline double LMCThreeVector::Magnitude() const
{
    return sqrt( fData[ 0 ] * fData[ 0 ] + fData[ 1 ] * fData[ 1 ] + fData[ 2 ] * fData[ 2 ] );
}
inline double LMCThreeVector::MagnitudeSquared() const
{
    return fData[ 0 ] * fData[ 0 ] + fData[ 1 ] * fData[ 1 ] + fData[ 2 ] * fData[ 2 ];
}

inline double LMCThreeVector::Perp() const
{
    return sqrt( fData[ 0 ] * fData[ 0 ] + fData[ 1 ] * fData[ 1 ] );
}
inline double LMCThreeVector::PerpSquared() const
{
    return fData[ 0 ] * fData[ 0 ] + fData[ 1 ] * fData[ 1 ];
}

inline double LMCThreeVector::PolarAngle() const
{
    return atan2( sqrt( fData[ 0 ] * fData[ 0 ] + fData[ 1 ] * fData[ 1 ] ), fData[ 2 ] );
}
inline double LMCThreeVector::AzimuthalAngle() const
{
    return atan2( fData[ 1 ], fData[ 0 ] );
}

inline LMCThreeVector LMCThreeVector::Unit() const
{
    const double tMagnitude = Magnitude();
    return tMagnitude > 0.0 ? LMCThreeVector( fData[ 0 ] / tMagnitude, fData[ 1 ] / tMagnitude, fData[ 2 ] / tMagnitude ) : LMCThreeVector( fData[ 0 ], fData[ 1 ], fData[ 2 ] );
}
inline LMCThreeVector LMCThreeVector::Orthogonal() const
{
    const double tX = fData[ 0 ] < 0.0 ? -fData[ 0 ] : fData[ 0 ];
    const double tY = fData[ 1 ] < 0.0 ? -fData[ 1 ] : fData[ 1 ];
    const double tZ = fData[ 2 ] < 0.0 ? -fData[ 2 ] : fData[ 2 ];
    if( tX < tY )
    {
        return tX < tZ ? LMCThreeVector( 0., fData[ 2 ], -fData[ 1 ] ) : LMCThreeVector( fData[ 1 ], -fData[ 0 ], 0. );
    }
    else
    {
        return tY < tZ ? LMCThreeVector( -fData[ 2 ], 0., fData[ 0 ] ) : LMCThreeVector( fData[ 1 ], -fData[ 0 ], 0. );
    }
}
inline LMCThreeVector LMCThreeVector::Cross( const LMCThreeVector& aVector ) const
{
    return LMCThreeVector( fData[ 1 ] * aVector.fData[ 2 ] - fData[ 2 ] * aVector.fData[ 1 ], fData[ 2 ] * aVector.fData[ 0 ] - fData[ 0 ] * aVector.fData[ 2 ], fData[ 0 ] * aVector.fData[ 1 ] - fData[ 1 ] * aVector.fData[ 0 ] );
}
inline bool LMCThreeVector::operator ==( const LMCThreeVector& aVector ) const
{
    return (aVector.fData[ 0 ] == fData[ 0 ] && aVector.fData[ 1 ] == fData[ 1 ] && aVector.fData[ 2 ] == fData[ 2 ]) ? true : false;
}
inline bool LMCThreeVector::operator !=( const LMCThreeVector& aVector ) const
{
    return (aVector.fData[ 0 ] != fData[ 0 ] || aVector.fData[ 1 ] != fData[ 1 ] || aVector.fData[ 2 ] != fData[ 2 ]) ? true : false;
}
inline bool LMCThreeVector::operator <( const LMCThreeVector& aVector ) const
{
    return Magnitude() < aVector.Magnitude();
}

inline LMCThreeVector operator+( const LMCThreeVector& aLeft, const LMCThreeVector& aRight )
{
    LMCThreeVector aResult( aLeft );
    aResult[ 0 ] += aRight[ 0 ];
    aResult[ 1 ] += aRight[ 1 ];
    aResult[ 2 ] += aRight[ 2 ];
    return aResult;
}
inline LMCThreeVector& operator+=( LMCThreeVector& aLeft, const LMCThreeVector& aRight )
{
    aLeft[ 0 ] += aRight[ 0 ];
    aLeft[ 1 ] += aRight[ 1 ];
    aLeft[ 2 ] += aRight[ 2 ];
    return aLeft;
}

inline LMCThreeVector operator-( const LMCThreeVector& aLeft, const LMCThreeVector& aRight )
{
    LMCThreeVector aResult( aLeft );
    aResult[ 0 ] -= aRight[ 0 ];
    aResult[ 1 ] -= aRight[ 1 ];
    aResult[ 2 ] -= aRight[ 2 ];
    return aResult;
}
inline LMCThreeVector& operator-=( LMCThreeVector& aLeft, const LMCThreeVector& aRight )
{
    aLeft[ 0 ] -= aRight[ 0 ];
    aLeft[ 1 ] -= aRight[ 1 ];
    aLeft[ 2 ] -= aRight[ 2 ];
    return aLeft;
}

inline double operator*( const LMCThreeVector& aLeft, const LMCThreeVector& aRight )
{
    return aLeft[ 0 ] * aRight[ 0 ] + aLeft[ 1 ] * aRight[ 1 ] + aLeft[ 2 ] * aRight[ 2 ];
}

inline LMCThreeVector operator*( const double aScalar, const LMCThreeVector& aVector )
{
    LMCThreeVector aResult( aVector );
    aResult[ 0 ] *= aScalar;
    aResult[ 1 ] *= aScalar;
    aResult[ 2 ] *= aScalar;
    return aResult;
}
inline LMCThreeVector operator*( const LMCThreeVector& aVector, const double aScalar )
{
    LMCThreeVector aResult( aVector );
    aResult[ 0 ] *= aScalar;
    aResult[ 1 ] *= aScalar;
    aResult[ 2 ] *= aScalar;
    return aResult;
}
inline LMCThreeVector& operator*=( LMCThreeVector& aVector, const double aScalar )
{
    aVector[ 0 ] *= aScalar;
    aVector[ 1 ] *= aScalar;
    aVector[ 2 ] *= aScalar;
    return aVector;
}
inline LMCThreeVector operator/( const LMCThreeVector& aVector, const double aScalar )
{
    LMCThreeVector aResult( aVector );
    aResult[ 0 ] /= aScalar;
    aResult[ 1 ] /= aScalar;
    aResult[ 2 ] /= aScalar;
    return aResult;
}
inline LMCThreeVector& operator/=( LMCThreeVector& aVector, const double aScalar )
{
    aVector[ 0 ] /= aScalar;
    aVector[ 1 ] /= aScalar;
    aVector[ 2 ] /= aScalar;
    return aVector;
}

}

#endif
