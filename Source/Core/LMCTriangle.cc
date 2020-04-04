#include "LMCTriangle.hh"

namespace locust
{
    const LMCTriangle gZeroTriangle(LMCThreeVector::sInvalid,LMCThreeVector::sInvalid,LMCThreeVector::sInvalid);
    const LMCTriangle gInvalidTriangle(LMCThreeVector::sZero,LMCThreeVector::sZero,LMCThreeVector::sZero);

    LMCTriangle::LMCTriangle()
    {
        for (int i=0; i<3; ++i)
        {
            fVertices[i]=LMCThreeVector();
        }
    }

    LMCTriangle::LMCTriangle(const LMCTriangle& aTriangle)
    {
        for (int i=0; i<3; ++i)
        {
            fVertices[i]=aTriangle.GetVertex(i);
        }
    }

    LMCTriangle::LMCTriangle(const LMCThreeVector& vertex0, const LMCThreeVector& vertex1, const LMCThreeVector& vertex2)
    {
        fVertices[0]=vertex0;
        fVertices[1]=vertex1;
        fVertices[2]=vertex2;
    }

    LMCTriangle::LMCTriangle(const double vertex0[3], const double vertex1[3], const double vertex2[3])
    {
        fVertices[0]=vertex0;
        fVertices[1]=vertex1;
        fVertices[2]=vertex2;
    }

    LMCTriangle::~LMCTriangle(){}

    LMCTriangle& LMCTriangle::operator=(const LMCTriangle& aTriangle)
    {
        for (int i=0; i<3; ++i)
        {
            fVertices[i]=aTriangle.GetVertex(i);
        }
        return *this;
    }
            
    LMCThreeVector& LMCTriangle::operator[](int vertexIndex)
    {
        return fVertices[vertexIndex];
    }

    const LMCThreeVector& LMCTriangle::operator[](int vertexIndex) const
    {
        return fVertices[vertexIndex];
    }

    bool LMCTriangle::operator==(const LMCTriangle& aTriangle) const
    {
        return ((aTriangle.GetVertex(0)==fVertices.at(0)) && (aTriangle.GetVertex(1)==fVertices.at(1)) && (aTriangle.GetVertex(2)==fVertices.at(2)));
    }

    bool LMCTriangle::operator!=(const LMCTriangle& aTriangle) const
    {
        return ((aTriangle.GetVertex(0)!=fVertices.at(0)) || (aTriangle.GetVertex(1)!=fVertices.at(1)) || (aTriangle.GetVertex(2)!=fVertices.at(2)));
    }

    LMCTriangle& LMCTriangle::operator+=(const LMCTriangle& aTriangle) 
    {
        for (int i=0; i<3; ++i)
        {
            fVertices.at(i)+=aTriangle.GetVertex(i);
        }
        return *this;
    }

    LMCTriangle& LMCTriangle::operator+=(const LMCThreeVector& aVector) 
    {
        for (int i=0; i<3; ++i)
        {
            fVertices.at(i)+=aVector;
        }
        return *this;
    }

    LMCTriangle& LMCTriangle::SetMagnitude(double magnitude)
    {
        for (int i=0; i<3; ++i)
        {
            std::cout<<"Before operator* "<<fVertices.at(i).Magnitude()<< " "<<fVertices.at(i).X()<< " "<<fVertices.at(i).Y()<<" "<<fVertices.at(i).Z() <<std::endl;
            fVertices.at(i).SetMagnitude(magnitude);
            std::cout<<"After operator* "<<fVertices.at(i).Magnitude()<< " "<<fVertices.at(i).X()<< " "<<fVertices.at(i).Y()<<" "<<fVertices.at(i).Z() <<std::endl;
        }
        return *this;
    }

    LMCTriangle& LMCTriangle::SetVertex(int vertexIndex,const LMCThreeVector& vertex)
    {
        fVertices[vertexIndex]=vertex;
        return *this;
    }

    const LMCThreeVector& LMCTriangle::GetVertex(int vertexIndex) const
    {
        return fVertices.at(vertexIndex);
    }
            
    const std::array<LMCThreeVector,3>& LMCTriangle::GetVertices(int vertexIndex) const
    {
        return fVertices;
    }

    LMCThreeVector LMCTriangle::GetSide01() const 
    {
        return fVertices.at(0)-fVertices.at(1);
    }

    LMCThreeVector LMCTriangle::GetSide02() const 
    {
        return fVertices.at(0)-fVertices.at(2);
    }

    LMCThreeVector LMCTriangle::GetSide12() const 
    {
        return fVertices.at(1)-fVertices.at(2);
    }

    void LMCTriangle::GetSides(LMCThreeVector& side01,LMCThreeVector& side02,LMCThreeVector& side12) const 
    {
        side01=GetSide01();
        side02=GetSide02();
        side12=GetSide12();
    }

    LMCThreeVector LMCTriangle::GetCenter() const
    {
        LMCThreeVector center((fVertices.at(0)+fVertices.at(1)+fVertices.at(2))/3.0);
        return center;
    }

    LMCThreeVector LMCTriangle::GetSideMidPoint01() const
    {
        LMCThreeVector point((fVertices.at(0)+fVertices.at(1))/2.0);
        return point;
    }

    LMCThreeVector LMCTriangle::GetSideMidPoint02() const
    {
        LMCThreeVector point((fVertices.at(0)+fVertices.at(2))/2.0);
        return point;
    }

    LMCThreeVector LMCTriangle::GetSideMidPoint12() const
    {
        LMCThreeVector point((fVertices.at(1)+fVertices.at(2))/2.0);
        return point;
    }

    void LMCTriangle::GetSideMidPoints(LMCThreeVector& midPoint01,LMCThreeVector& midPoint02,LMCThreeVector& midPoint12) const
    {
        midPoint01=GetSideMidPoint01();
        midPoint02=GetSideMidPoint02();
        midPoint12=GetSideMidPoint12();
        return;
    }

    double LMCTriangle::GetArea() const
    {
        LMCThreeVector side01=GetSide01();
        LMCThreeVector side02=GetSide02();
        return (side01.Cross(side02)).Magnitude()/2.0;
    }

    const LMCThreeVector LMCTriangle::GetNormal() const
    {
        LMCThreeVector side01=GetSide01();
        LMCThreeVector side02=GetSide02();
        // need to conform the direction, there are two possibilities
        return (side01.Cross(side02)).Unit();
    }

    bool LMCTriangle::Quadrasect(std::vector<LMCTriangle>& dividedTriangles) const
    {
        dividedTriangles.push_back(LMCTriangle(fVertices.at(0),GetSideMidPoint01(),GetSideMidPoint02()));
        dividedTriangles.push_back(LMCTriangle(GetSideMidPoint01(),fVertices.at(1),GetSideMidPoint12()));
        dividedTriangles.push_back(LMCTriangle(GetSideMidPoint02(),GetSideMidPoint12(),fVertices.at(2)));
        dividedTriangles.push_back(LMCTriangle(GetSideMidPoint01(),GetSideMidPoint01(),GetSideMidPoint02()));
        return true;
    }

    bool LMCTriangle::IsEquilateralTriangle() const
    {
        return (GetSide01().Magnitude()==GetSide02().Magnitude() && GetSide01().Magnitude()==GetSide12().Magnitude()); 
    }
}
