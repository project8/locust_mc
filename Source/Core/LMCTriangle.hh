#ifndef LMCTRIANGLE_H_
#define LMCTRIANGLE_H_

#include "LMCThreeVector.hh"

#include <array>
#include <vector>

namespace locust
{
    class LMCTriangle 
    {
        public:
            static const LMCTriangle gZeroTriangle;
            static const LMCTriangle gInvalidTriangle;

            LMCTriangle();
            //Constructors
            LMCTriangle(const LMCTriangle& aTriangle);
            LMCTriangle(const LMCThreeVector& vertex0, const LMCThreeVector& vertex1, const LMCThreeVector& vertex2);
            LMCTriangle(const double vertex0[0],const double vertex1[3],const double vertex2[3]);
            
            //Destructors
            ~LMCTriangle();

            // operators 
            LMCTriangle& operator=(const LMCTriangle& aTriangle);
            LMCThreeVector& operator[](int vertexIndex);
            const LMCThreeVector& operator[](int vertexIndex) const;
            bool operator==(const LMCTriangle& aTriangle) const;
            bool operator!=(const LMCTriangle& aTriangle) const;
            LMCTriangle& operator+=(const LMCTriangle& aTriangle);
            LMCTriangle& operator+=(const LMCThreeVector& aVector);

            // Setters for vertices
            LMCTriangle& SetVertex(int vertexIndex,const LMCThreeVector& vertex);
            LMCTriangle& SetMagnitude(double magnitude);
            LMCTriangle& Move(const LMCThreeVector& displacementVector);
            // Getter for vertices
            const LMCThreeVector& GetVertex(int vertexIndex) const;
            const std::array<LMCThreeVector,3>& GetVertices(int vertexIndex) const;

            // Other getters 
            LMCThreeVector GetCenter() const;
            LMCThreeVector GetSide01() const;
            LMCThreeVector GetSide02() const;
            LMCThreeVector GetSide12() const;
            void GetSides(LMCThreeVector& side01,LMCThreeVector& side02,LMCThreeVector& side12) const;
            LMCThreeVector GetSideMidPoint01() const;
            LMCThreeVector GetSideMidPoint02() const;
            LMCThreeVector GetSideMidPoint12() const;
            void GetSideMidPoints(LMCThreeVector& midPoint01,LMCThreeVector& midPoint02,LMCThreeVector& midPoint12) const;
            double GetLength01() const;
            double GetLength02() const;
            double GetLength12() const;

            double GetArea() const; 
            const LMCThreeVector GetNormal() const; 
            
            bool Quadrasect(std::vector<LMCTriangle>& dividedTriangles) const;
            bool IsEquilateralTriangle() const;

        private:
            std::array<LMCThreeVector,3> fVertices;
    };
}

#endif
