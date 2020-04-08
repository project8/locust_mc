#include "LMCSphere.hh"
#include <iostream>
#include <fstream>

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
        fRadius=abs(radius);
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
        fMinimumFaces=0.0;
        ConstructSphere();
    }

    LMCIcoSphere::LMCIcoSphere(double radius)
    {
        fRadius=abs(radius);
        fCenter=LMCThreeVector();
        fMinimumFaces=0.0;
    }

    LMCIcoSphere::LMCIcoSphere(double radius, LMCThreeVector center)
    {
        fRadius=abs(radius);
        fCenter=center;
        fMinimumFaces=0.0;
        ConstructSphere();
    }

    LMCIcoSphere::LMCIcoSphere(double radius, LMCThreeVector center, unsigned minFaces)
    {
        fRadius=abs(radius);
        fCenter=center;
        fMinimumFaces=minFaces;
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
        fRadius=abs(radius);
        ConstructSphere();
    }

    void LMCIcoSphere::SetCenter(const LMCThreeVector& aPoint)
    {
        fCenter=aPoint;
        ConstructSphere();
    }

    void LMCIcoSphere::SetMinFaces(unsigned minFaces)
    {
        fMinimumFaces=minFaces;
        ConstructSphere();
    }

    const std::vector<LMCTriangle>& LMCIcoSphere::GetFaces() const 
    {
        return fTriangleFaces;
    }

    bool LMCIcoSphere::GetFaceCenters(std::vector<LMCThreeVector>& faceCenters) const 
    {
        for(auto&& vertex:fTriangleFaces)
        {
            faceCenters.push_back(vertex.GetCenter());
        }
        if(faceCenters.size()!=fTriangleFaces.size()) false;
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
        std::vector<LMCThreeVector> sphereVertices;
        sphereVertices.push_back(LMCThreeVector(-s,t,0));
        sphereVertices.push_back(LMCThreeVector(s,t,0));
        sphereVertices.push_back(LMCThreeVector(-s,-t,0));
        sphereVertices.push_back(LMCThreeVector(s,-t,0));

        sphereVertices.push_back(LMCThreeVector(0,-s,t));
        sphereVertices.push_back(LMCThreeVector(0,s,t));
        sphereVertices.push_back(LMCThreeVector(0,-s,-t));
        sphereVertices.push_back(LMCThreeVector(0,s,-t));

        sphereVertices.push_back(LMCThreeVector(t,0,-s));
        sphereVertices.push_back(LMCThreeVector(t,0,s));
        sphereVertices.push_back(LMCThreeVector(-t,0,-s));
        sphereVertices.push_back(LMCThreeVector(-t,0,s));

        for(auto&& vertex:sphereVertices)
        {
            std::cout<<"vertex.Magnitude() "<<vertex.Magnitude()<<std::endl;
        }
        //First build the triangle faces joining the vertices always going in the anticlockwise direction looking from outside into the icosahedron
        fTriangleFaces.clear();
        // Start with all the triangles having vvertex 0 in common
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(0),sphereVertices.at(11),sphereVertices.at(5)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(0),sphereVertices.at(5),sphereVertices.at(1)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(0),sphereVertices.at(1),sphereVertices.at(7)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(0),sphereVertices.at(7),sphereVertices.at(10)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(0),sphereVertices.at(10),sphereVertices.at(11)));

        // Add triangles that share a side or vertex with the previous ones 
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(11),sphereVertices.at(10),sphereVertices.at(2)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(2),sphereVertices.at(4),sphereVertices.at(11)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(5),sphereVertices.at(11),sphereVertices.at(4)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(4),sphereVertices.at(9),sphereVertices.at(5)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(1),sphereVertices.at(5),sphereVertices.at(9)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(9),sphereVertices.at(8),sphereVertices.at(1)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(7),sphereVertices.at(1),sphereVertices.at(8)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(6),sphereVertices.at(2),sphereVertices.at(10)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(10),sphereVertices.at(7),sphereVertices.at(6)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(8),sphereVertices.at(6),sphereVertices.at(7)));

        //All triangles having vertix 3 common 
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(3),sphereVertices.at(2),sphereVertices.at(6)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(3),sphereVertices.at(4),sphereVertices.at(2)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(3),sphereVertices.at(9),sphereVertices.at(4)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(3),sphereVertices.at(8),sphereVertices.at(9)));
        fTriangleFaces.push_back(LMCTriangle(sphereVertices.at(3),sphereVertices.at(6),sphereVertices.at(8)));

    }

    bool LMCIcoSphere::AmendToRadius(std::vector<LMCTriangle>&triangles)
    {
        for(auto&& triangle:triangles)
        {
            triangle.SetMagnitude(fRadius);
        }
        return true;
    }

     int GetVertexIndex(LMCThreeVector& point) const
     {

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
        if(nInitialTriangles>fMinimumFaces) return true;
        std::vector<LMCTriangle> tempTriangleFaces;
        for(unsigned i=0;i<nInitialTriangles;++i)
        {
            std::vector<LMCTriangle> dividedTriangles;
            /*
            if(! fTriangleFaces.at(i).Quadrasect(dividedTriangles)) return false;
            AmendToRadius(dividedTriangles);
            tempTriangleFaces.insert(tempTriangleFaces.end(),dividedTriangles.begin(),dividedTriangles.end());
            if(tempTriangleFaces.size()>=fMaxFaces) return false;
            */
            point01=fTriangleFaces.at(i).GetSideMidPoint01();
            point02=fTriangleFaces.at(i).GetSideMidPoint02();
            point03=fTriangleFaces.at(i).GetSideMidPoint12();
            int index01=GetVertexIndex(point01);
            int index02=GetVertexIndex(point02);
            int index02=GetVertexIndex(point12);
            dividedTriangles.push_back(fTriangleFaces.at(i).GetVertex(0),sphereVertices.at(index01),sphereVertices.at(index02));
            dividedTriangles.push_back(sphereVertices.at(index01),fTriangleFaces.at(i).GetVertex(1),sphereVertices.at(index12));
            dividedTriangles.push_back(sphereVertices.at(index02),sphereVertices.at(index12),fTriangleFaces.at(i).GetVertex(2));
            dividedTriangles.push_back(sphereVertices.at(index01),sphereVertices.at(index12),sphereVertices.at(index02));
            AmendToRadius(dividedTriangles);
            tempTriangleFaces.insert(tempTriangleFaces.end(),dividedTriangles.begin(),dividedTriangles.end());
            if(tempTriangleFaces.size()>=fMaxFaces) return false;
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
        for(unsigned i=0;i<fTriangleFaces.size();++i)
        {
            for(unsigned k=0;k<3;++k)
            {
                LMCThreeVector vertex=fTriangleFaces.at(i).GetVertex(k);
                myfile<<vertex.X()<< ","<<vertex.Y()<<","<<vertex.Z()<<"\t";
            }
            myfile<<"\n";
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

    //Logic used from: (similar construction alos used in Eigen )
    //http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
    bool LMCIcoSphere::ConstructSphere() 
    {
        GenerateIcosahedron();
        DivideTriangleFaces();
        WriteVertices();
        WriteFaceCenters();
        //MoveCenter();
        std::cout<< " --------------------------------CENTER NOT SET YET --------------------------------"<<std::endl; 
        return true;
    }
}
