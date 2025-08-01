#include <gmsh.h>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <fstream>

// Write temperatures to a Gmsh .pos file for visualization
static void write_pos(const std::string &name,
                      const std::vector<double> &coords,
                      const Eigen::VectorXd &temp)
{
    std::ofstream file(name);
    file << "View \"Temperature\" {\n";
    for (std::size_t i = 0; i < temp.size(); ++i)
    {
        file << "SP(" << coords[3*i] << "," << coords[3*i+1] << "," << coords[3*i+2]
             << "){" << temp[i] << "};\n";
    }
    file << "};\n";
}

static double tri_area(double x1,double y1,double x2,double y2,double x3,double y3)
{
    return 0.5*std::abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
}

static Eigen::Matrix3d local_stiffness(const double c[6], double alpha)
{
    double x1=c[0],y1=c[1],x2=c[2],y2=c[3],x3=c[4],y3=c[5];
    double area = tri_area(x1,y1,x2,y2,x3,y3);
    double b1=y2-y3,c1=x3-x2;
    double b2=y3-y1,c2=x1-x3;
    double b3=y1-y2,c3=x2-x1;
    Eigen::Matrix3d k; k.setZero();
    double f = alpha/(4.0*area);
    k(0,0)=f*(b1*b1+c1*c1); k(0,1)=f*(b1*b2+c1*c2); k(0,2)=f*(b1*b3+c1*c3);
    k(1,0)=k(0,1); k(1,1)=f*(b2*b2+c2*c2); k(1,2)=f*(b2*b3+c2*c3);
    k(2,0)=k(0,2); k(2,1)=k(1,2); k(2,2)=f*(b3*b3+c3*c3);
    return k;
}

int main()
{
    // Parameters
    double alpha = 1e-4;           // thermal diffusivity
    double dt = 0.1;               // time step
    double total_time = 10.0;      // total simulation time
    double source_temp = 200.0;    // temperature of heaters
    double initial_temp = 20.0;    // room temperature
    double heater_radius = 0.25;   // radius of heated region

    gmsh::initialize();
    gmsh::open("stove.msh");

    // Nodes
    std::vector<std::size_t> nodeTags; std::vector<double> nodeCoords, param;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, param);
    std::size_t nNodes = nodeTags.size();
    std::map<std::size_t,std::size_t> tagToIndex;
    for(std::size_t i=0;i<nNodes;++i) tagToIndex[nodeTags[i]] = i;

    // Triangular elements
    std::vector<int> types; std::vector<std::vector<std::size_t>> elemTags, elemNodeTags;
    gmsh::model::mesh::getElements(types, elemTags, elemNodeTags);
    std::vector<std::size_t> triNodes;
    for(std::size_t i=0;i<types.size();++i)
        if(types[i]==2){ triNodes = elemNodeTags[i]; break; }
    std::size_t nTri = triNodes.size()/3;

    Eigen::SparseMatrix<double> K(nNodes,nNodes);
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::VectorXd mass = Eigen::VectorXd::Zero(nNodes);

    for(std::size_t e=0;e<nTri;++e){
        std::size_t n1=triNodes[3*e];
        std::size_t n2=triNodes[3*e+1];
        std::size_t n3=triNodes[3*e+2];
        std::size_t i1=tagToIndex[n1];
        std::size_t i2=tagToIndex[n2];
        std::size_t i3=tagToIndex[n3];
        double c[6]={nodeCoords[3*i1],nodeCoords[3*i1+1],
                     nodeCoords[3*i2],nodeCoords[3*i2+1],
                     nodeCoords[3*i3],nodeCoords[3*i3+1]};
        Eigen::Matrix3d k = local_stiffness(c,alpha);
        double area = tri_area(c[0],c[1],c[2],c[3],c[4],c[5]);
        mass(i1)+=area/3.0; mass(i2)+=area/3.0; mass(i3)+=area/3.0;
        std::array<std::size_t,3> idx={i1,i2,i3};
        for(int a=0;a<3;++a)
            for(int b=0;b<3;++b)
                triplets.emplace_back(idx[a],idx[b],k(a,b));
    }
    K.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd T = Eigen::VectorXd::Constant(nNodes, initial_temp);
    Eigen::VectorXd Tnew(nNodes);
    std::vector<std::pair<double,double>> centers={{3,1},{7,1},{3,5},{7,5}};

    write_pos("temperature_0.pos", nodeCoords, T);
    int steps = static_cast<int>(total_time/dt);
    for(int s=1;s<=steps;++s){
        Eigen::VectorXd lap = K*T;
        Tnew = T + dt * (-lap.array()/mass.array()).matrix();
        for(std::size_t i=0;i<nNodes;++i){
            double x=nodeCoords[3*i], y=nodeCoords[3*i+1];
            for(const auto &c:centers){
                double dx=x-c.first, dy=y-c.second;
                if(std::sqrt(dx*dx+dy*dy)<heater_radius){ Tnew(i)=source_temp; break; }
            }
        }
        T.swap(Tnew);
        if(s%10==0) write_pos("temperature_"+std::to_string(s)+".pos", nodeCoords, T);
    }

    gmsh::finalize();
    return 0;
}
