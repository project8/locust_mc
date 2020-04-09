#include "LMCSphere.hh"
#include "logger.hh"
#include <iostream>
#include <fstream>

namespace locust
{
    LOGGER( lmclog, "LMCSphere" );
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
        fRadius=fabs(radius);
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

    LMCIcoSphere::LMCIcoSphere()
    {
        fRadius=0.0;
        fCenter=LMCThreeVector();
        fMinimumVertices=0.0;
        ConstructSphere();
    }

    LMCIcoSphere::LMCIcoSphere(double radius)
    {
        fRadius=fabs(radius);
        fCenter=LMCThreeVector();
        fMinimumVertices=0.0;
    }

    LMCIcoSphere::LMCIcoSphere(double radius, LMCThreeVector center)
    {
        fRadius=fabs(radius);
        fCenter=center;
        fMinimumVertices=0.0;
        ConstructSphere();
    }

    LMCIcoSphere::LMCIcoSphere(double radius, LMCThreeVector center, unsigned minVertices)
    {
        fRadius=fabs(radius);
        fCenter=center;
        fMinimumVertices=minVertices;
        ConstructSphere();
    }

    LMCIcoSphere::LMCIcoSphere(const LMCIcoSphere& anIcosphere)
    {
        fRadius=anIcosphere.GetRadius();
        fCenter=anIcosphere.GetCenter();
        fTriangleFaces=anIcosphere.GetFaces();
    }

    void LMCIcoSphere::SetRadius(double radius)
    {
        fRadius=fabs(radius);
        ConstructSphere();
    }

    void LMCIcoSphere::SetCenter(const LMCThreeVector& aPoint)
    {
        fCenter=aPoint;
        ConstructSphere();
    }

    void LMCIcoSphere::SetMinVertices(unsigned minVertices)
    {
        fMinimumVertices=minVertices;
        ConstructSphere();
    }

    const std::vector<LMCTriangle>& LMCIcoSphere::GetFaces() const 
    {
        return fTriangleFaces;
    }

    bool LMCIcoSphere::GetFaceCenters(std::vector<LMCThreeVector>& faceCenters) const 
    {
        for(auto&& triangle:fTriangleFaces)
        {
            faceCenters.push_back(triangle.GetCenter());
        }
        if(faceCenters.size()!=fTriangleFaces.size()) return false;
        return true;
    }

    bool LMCIcoSphere::GetVertices(std::vector<LMCThreeVector>& vertices) const 
    {
        for(auto&& vertex:fSphereVertices)
        {
            vertices.push_back(vertex);
        }
        return true;
    }

    unsigned LMCIcoSphere::GetNFaces() const 
    {
        return fTriangleFaces.size();
    }

    void LMCIcoSphere::GenerateIcosahedron()
    {
        double s = fRadius*sqrt((5.0-sqrt(5.0))/10.0);  
        double t = fRadius*sqrt((5.0+sqrt(5.0))/10.0);  

        //First build a unit icosahedron
        fSphereVertices.push_back(LMCThreeVector(-s,t,0));
        fSphereVertices.push_back(LMCThreeVector(s,t,0));
        fSphereVertices.push_back(LMCThreeVector(-s,-t,0));
        fSphereVertices.push_back(LMCThreeVector(s,-t,0));

        fSphereVertices.push_back(LMCThreeVector(0,-s,t));
        fSphereVertices.push_back(LMCThreeVector(0,s,t));
        fSphereVertices.push_back(LMCThreeVector(0,-s,-t));
        fSphereVertices.push_back(LMCThreeVector(0,s,-t));

        fSphereVertices.push_back(LMCThreeVector(t,0,-s));
        fSphereVertices.push_back(LMCThreeVector(t,0,s));
        fSphereVertices.push_back(LMCThreeVector(-t,0,-s));
        fSphereVertices.push_back(LMCThreeVector(-t,0,s));

        //First build the triangle faces joining the vertices always going in the anticlockwise direction looking from outside into the icosahedron
        fTriangleFaces.clear();
        // Start with all the triangles having vvertex 0 in common
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(0),fSphereVertices.at(11),fSphereVertices.at(5)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(0),fSphereVertices.at(5),fSphereVertices.at(1)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(0),fSphereVertices.at(1),fSphereVertices.at(7)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(0),fSphereVertices.at(7),fSphereVertices.at(10)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(0),fSphereVertices.at(10),fSphereVertices.at(11)));

        // Add triangles that share a side or vertex with the previous ones 
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(11),fSphereVertices.at(10),fSphereVertices.at(2)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(2),fSphereVertices.at(4),fSphereVertices.at(11)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(5),fSphereVertices.at(11),fSphereVertices.at(4)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(4),fSphereVertices.at(9),fSphereVertices.at(5)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(1),fSphereVertices.at(5),fSphereVertices.at(9)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(9),fSphereVertices.at(8),fSphereVertices.at(1)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(7),fSphereVertices.at(1),fSphereVertices.at(8)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(6),fSphereVertices.at(2),fSphereVertices.at(10)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(10),fSphereVertices.at(7),fSphereVertices.at(6)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(8),fSphereVertices.at(6),fSphereVertices.at(7)));

        //All triangles having vertix 3 common 
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(3),fSphereVertices.at(2),fSphereVertices.at(6)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(3),fSphereVertices.at(4),fSphereVertices.at(2)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(3),fSphereVertices.at(9),fSphereVertices.at(4)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(3),fSphereVertices.at(8),fSphereVertices.at(9)));
        fTriangleFaces.push_back(LMCTriangle(fSphereVertices.at(3),fSphereVertices.at(6),fSphereVertices.at(8)));

    }

    int LMCIcoSphere::GetVertexIndex(LMCThreeVector& point,double threshold) 
    {
        for(int i=0;i<fSphereVertices.size();++i)
        {
            if((fSphereVertices.at(i)-point).Magnitude()<threshold)
            {
                return i;
            }
        }
        point.SetMagnitude(fRadius);
        fSphereVertices.push_back(point);
        return fSphereVertices.size()-1;
    }

    /*Logic in DivideTriangleFaces:
      For each triangle:
      For each set of two points
      Find the midpoint
      If the midpoint exists in SphereVertices
      Get the index and add it to the triangle 
      If doesn't exist, add it to SphereVertices and use the index to add to the triangle
      */
    bool LMCIcoSphere::DivideTriangleFaces()
    {
        unsigned nInitialTriangles=fTriangleFaces.size();
        if(fSphereVertices.size()>fMinimumVertices) return true;
        std::vector<LMCTriangle> tempTriangleFaces;
        //adhoc value to make sure there are no 
        double lengthThreshold=fTriangleFaces.at(0).GetLength01()/100.0;
        for(unsigned i=0;i<nInitialTriangles;++i)
        {
            std::vector<LMCTriangle> dividedTriangles;
            /*
               if(! fTriangleFaces.at(i).Quadrasect(dividedTriangles)) return false;
               AmendToRadius(dividedTriangles);
               tempTriangleFaces.insert(tempTriangleFaces.end(),dividedTriangles.begin(),dividedTriangles.end());
               if(tempTriangleFaces.size()>=fMaxFaces) return false;
               */
            LMCThreeVector point01=fTriangleFaces.at(i).GetSideMidPoint01();
            LMCThreeVector point02=fTriangleFaces.at(i).GetSideMidPoint02();
            LMCThreeVector point12=fTriangleFaces.at(i).GetSideMidPoint12();
            int index01=GetVertexIndex(point01,lengthThreshold);
            int index02=GetVertexIndex(point02,lengthThreshold);
            int index12=GetVertexIndex(point12,lengthThreshold);
            dividedTriangles.push_back(LMCTriangle(fTriangleFaces.at(i).GetVertex(0),fSphereVertices.at(index01),fSphereVertices.at(index02)));
            dividedTriangles.push_back(LMCTriangle(fSphereVertices.at(index01),fTriangleFaces.at(i).GetVertex(1),fSphereVertices.at(index12)));
            dividedTriangles.push_back(LMCTriangle(fSphereVertices.at(index02),fSphereVertices.at(index12),fTriangleFaces.at(i).GetVertex(2)));
            dividedTriangles.push_back(LMCTriangle(fSphereVertices.at(index01),fSphereVertices.at(index12),fSphereVertices.at(index02)));
    //        AmendToRadius(dividedTriangles);
            tempTriangleFaces.insert(tempTriangleFaces.end(),dividedTriangles.begin(),dividedTriangles.end());
            if(fSphereVertices.size()>=fMaximumVertices) return false;
        }
        fTriangleFaces.erase(fTriangleFaces.begin(),fTriangleFaces.end());
        fTriangleFaces=tempTriangleFaces;
        DivideTriangleFaces();
        return true;
    }

    bool LMCIcoSphere::WriteVertices() const 
    {
        std::ofstream myfile;
        myfile.open("SphereVertices.txt");
        myfile<<"# Vertices of the spheres \n";
        myfile<<"# Each row has three sets of three numbers for the vertices of the faces of the triangles \n";
        for(auto&& vertex:fSphereVertices)
        {
            myfile<<vertex.X()<< ","<<vertex.Y()<<","<<vertex.Z()<<"\n";
        }
        myfile.close();
        return true;
    }

    bool LMCIcoSphere::WriteFaceCenters() const 
    {
        std::ofstream myfile;
        myfile.open("SphereFaceCenters.txt");
        myfile<<"# Face centers of triangles forming the spheres\n";
        myfile<<"# Each row has three numbers for the centers of rhe faces of the triangles \n";
        for(unsigned i=0;i<fTriangleFaces.size();++i)
        {
            LMCThreeVector faceCenter=fTriangleFaces.at(i).GetCenter();
            myfile<<faceCenter.X()<< ","<<faceCenter.Y()<<","<<faceCenter.Z()<<"\n";
        }
        myfile.close();
        return true;
    }
    
    bool LMCIcoSphere::DisplaceCenter()
    {
        for(int i=0;i<fSphereVertices.size();i++)
        {
            fSphereVertices.at(i)+=fCenter;
        }
        for(int i=0;i<fTriangleFaces.size();i++)
        {
            for(int i=0;i<3;i++)
            {
                fTriangleFaces.at(i).Move(fCenter);
            }
        }
        return true;
    }

    //Logic used from: (similar construction alos used in Eigen )
    //http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
    bool LMCIcoSphere::ConstructSphere() 
    {
        if(fRadius==0) 
        {
            LERROR(lmclog,"Radius of sphere can't be 0");
            exit(-1);
        }
        GenerateIcosahedron();
        DivideTriangleFaces();
        DisplaceCenter();
        WriteVertices();
        WriteFaceCenters();
        return true;
    }
}
