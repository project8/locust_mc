#ifndef LMCSPHERE_H_
#define LMCSPHERE_H_

#include "LMCTriangle.hh"
#include "LMCConst.hh"

#include <vector>

namespace locust
{
    class LMCSphere 
    {
        public:
            LMCSphere();
            LMCSphere(const LMCSphere& aSphere);
            ~LMCSphere(){}

            virtual LMCSphere& operator=(const LMCSphere& aSphere);
            virtual bool operator==(const LMCSphere& aSphere) const;
            virtual bool operator!=(const LMCSphere& aSphere) const;

            virtual const LMCThreeVector GetCenter() const; 
            virtual double GetRadius() const; 
            double GetSurfaceArea() const; 
            double GetVolume() const; 
            
            bool IsPointOnSphere(const LMCThreeVector& aPoint) const;
            bool IsPointInside(const LMCThreeVector& aPoint) const;

        protected:
            double fRadius;
            LMCThreeVector fCenter;
            virtual bool ConstructSphere()=0;
            virtual void SetRadius(double radius);
            virtual void SetCenter(const LMCThreeVector& aPoint);
    };

    class LMCIcoSphere:public LMCSphere 
    {
        public:
            /// Constructs and destructors, also generates the full sphere when a constructor is called 
            LMCIcoSphere();
            LMCIcoSphere(double radius);
            LMCIcoSphere(double radius, LMCThreeVector center);
            LMCIcoSphere(double radius, LMCThreeVector center, unsigned minFaces);
            LMCIcoSphere(const LMCIcoSphere& anIcosphere);
            ~LMCIcoSphere(){}
            
            /// Setters
            /// Set the radius of the sphere, overrides the parent function because this function also need to call the ConstructSphere function 
            void SetRadius(double radius) override;
            /// Set the center of the sphere, overrides the parent function because this function also need to call the ConstructSphere function 
            /// To modify the center after ConstructSphere is already called, call DisplaceCenter instead of this function
            void SetCenter(const LMCThreeVector& aPoint) override;
            /// Set the minimum number of vertices to be created when ConstructSphere is called, overrides the parent function because this function also need to call the ConstructSphere function 
            void SetMinVertices(unsigned minVertices);
            /// Displaces the center of the sphere once the ConstructSphere method is already called 
            bool DisplaceCenter();
            
            /// Get the number of faces 
            unsigned GetNFaces() const;
            /// Get a vector of LMCTriangles corresponding to the faces 
            const std::vector<LMCTriangle>& GetFaces() const;
            /// Get a vector of LMCThreeVectors corresponding to the face centers
            bool GetFaceCenters(std::vector<LMCThreeVector>& faceCenters) const;
            /// Get a vector of LMCThreeVectors corresponding to the vertices 
            bool GetVertices(std::vector<LMCThreeVector>& vertices) const;

        protected:
            /// Member variables
            /// The minimum number of vertices the sphere can be represented by
            //The actual number of vertices will always be higher than this number unless this number is higher than fMaximumVertices
            unsigned fMinimumVertices;
            // The actual meat of the class, the location of faces of the triangle
            std::vector<LMCTriangle> fTriangleFaces;
            // The actual meat of the class, the location of vertices
            std::vector<LMCThreeVector> fSphereVertices;

            /// Member functions
            /// Once a tringle is subdivided into 4 triangles, the radius has to be modified to match with the fRadius otherwise, it will just be an icosahedron and not sphere
            bool AmendToRadius(std::vector<LMCTriangle>& triangles);
            /// The index of the vertex, primarily used to make sure there is no repition of the vertices
            int GetVertexIndex(LMCThreeVector& point, double threshold) ;
            /// Divide each triangle into 4 sub triangles making sure all the vertices are at a unit distance 
            bool DivideTriangleFaces();
            void GenerateIcosahedron();
            bool WriteVertices() const; 
            bool WriteFaceCenters() const; 
            bool ConstructSphere() override;
            const unsigned fMaximumVertices=1000000;
    };
}
#endif
