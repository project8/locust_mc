#ifndef LMCTRIANGLE_H_
#define LMCTRIANGLE_H_

#include "LMCThreeVector.hh"

#include <array>

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

            // Getter for vertices
            const LMCThreeVector& GetVertex(int vertexIndex) const;
            const std::array<LMCThreeVector,3>& GetVertices(int vertexIndex) const;

            // Other getters 
            const LMCThreeVector GetMidPoint01() const;
            const LMCThreeVector GetMidPoint02() const;
            const LMCThreeVector GetMidPoint12() const;
            const void GetMidPoints(LMCThreeVector& midPoint01,LMCThreeVector& midPoint02,LMCThreeVector& midPoint12) const;

            double GetArea() const; 
            const LMCThreeVector GetNormal() const; 

        private:
            std::array<LMCThreeVector,3> fVertices;
    };
}

#endif
