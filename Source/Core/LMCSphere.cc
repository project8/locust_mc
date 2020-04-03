#include "LMCSphere.hh"

namespace locust
{
    LMCSphere::LMCSphere()
    {
        fRadius=0.0;
        fCenter=LMCThreeVector();
    }

    LMCSphere::LMCSphere(const LMCSphere& aSphere)
    {
        fRadius=aSphere.GetRadius();
        fCenter=aSphere.GetCenter();
    }

    LMCSphere& LMCSphere::operator=(const LMCSphere& aSphere)
    {
        fRadius=aSphere.GetRadius();
        fCenter=aSphere.GetCenter();
        return *this;
    }
            
    bool LMCSphere::operator==(const LMCSphere& aSphere) const
    {
        return (fRadius==aSphere.GetRadius() && fCenter==aSphere.GetCenter());
    }

    bool LMCSphere::operator!=(const LMCSphere& aSphere) const
    {
        return (fRadius!=aSphere.GetRadius() || fCenter!=aSphere.GetCenter());
    }

    void LMCSphere::SetRadius(double radius)
    {
        fRadius=radius;
    }

    void LMCSphere::SetCenter(const LMCThreeVector& aPoint)
    {
        fCenter=aPoint;
    }

    double LMCSphere::GetRadius() const 
    {
        return fRadius;
    }

    const LMCThreeVector LMCSphere::GetCenter() const
    {
        return fCenter;
    }

    double LMCSphere::GetSurfaceArea() const
    {
        return 4.0*LMCConst::Pi()*pow(fRadius,2.0);
    }

    double LMCSphere::GetVolume() const
    {
        return 4.0*LMCConst::Pi()*pow(fRadius,3.0)/3.0;
    }

    bool LMCSphere::IsPointOnSphere(const LMCThreeVector& aPoint) const
    {
        return ((aPoint-fCenter).Magnitude()==fRadius);
    }

    bool LMCSphere::IsPointInside(const LMCThreeVector& aPoint) const
    {
        return ((aPoint-fCenter).Magnitude()<fRadius);
    }

    LMCIcosphere::LMCIcosphere()
    {
        //Implement me
    }

    LMCIcosphere::LMCIcosphere(const LMCIcosphere& anIcosphere)
    {
        //Implement me
    }

    bool LMCIcosphere::ConstructSphere() 
    {
        return true;
        //Implement me
    }
}
