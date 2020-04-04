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
            LMCIcoSphere();
            LMCIcoSphere(double radius);
            LMCIcoSphere(double radius, LMCThreeVector center);
            LMCIcoSphere(double radius, LMCThreeVector center, unsigned minFaces);
            LMCIcoSphere(const LMCIcoSphere& anIcosphere);
            ~LMCIcoSphere(){}
            
            unsigned GetNFaces() const;
            void SetRadius(double radius) override;
            void SetCenter(const LMCThreeVector& aPoint) override;
            void SetMinFaces(unsigned minFaces);
            
            const std::vector<LMCTriangle>& GetFaces() const;
            bool GetFaceCenters(std::vector<LMCThreeVector>& faceCenters) const;

        protected:
            unsigned fMinimumFaces;
            std::vector<LMCTriangle> fTriangleFaces;
            bool AmendToRadius(std::vector<LMCTriangle>&triangles);
            bool DivideTriangleFaces();
            void GenerateIcosahedron();
            bool ConstructSphere() override;
            const unsigned fMaxFaces=1000;
    };
}

#endif
