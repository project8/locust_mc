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

            const LMCThreeVector GetCenter() const; 
            double GetRadius() const; 
            double GetSurfaceArea() const; 
            double GetVolume() const; 
            
            bool IsPointOnSphere(const LMCThreeVector& aPoint) const;
            bool IsPointInside(const LMCThreeVector& aPoint) const;

        protected:
            double fRadius;
            LMCThreeVector fCenter;
            virtual bool ConstructSphere()=0;
            void SetRadius(double radius);
            void SetCenter(const LMCThreeVector& aPoint);
    };

    class LMCIcosphere:public LMCSphere 
    {
        public:
            LMCIcosphere();
            LMCIcosphere(const LMCIcosphere& anIcosphere);
            ~LMCIcosphere(){}

        protected:
            std::vector<LMCTriangle> fTriangles;
            int nIterations=0;
            bool ConstructSphere() override;
    };
}

#endif
