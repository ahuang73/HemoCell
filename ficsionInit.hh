#ifndef FICSIONINIT_HH
#define FICSIONINIT_HH

#include <limits>

#include "ficsionInit.h"
#include "palabos3D.hh"

//using namespace plb;
//using namespace std;

/* ************* Functions poiseuillePressure and poiseuilleVelocity ******************* */
static T poiseuillePressure(IncomprFlowParam<T> const &parameters, plint maxN)
{
    const T a = parameters.getNy()-1;
    const T b = parameters.getNz()-1;
    const T nu = parameters.getLatticeNu();
    const T uMax = parameters.getLatticeU();
    T sum = T();
    for (plint iN = 0; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum += ((T)1 / (pow(twoNplusOne,3)*cosh(twoNplusOne*pi*b/((T)2*a))));
    }
    for (plint iN = 1; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum -= ((T)1 / (pow(twoNplusOne,3)*cosh(twoNplusOne*pi*b/((T)2*a))));
    }
    T alpha = -(T)8 * uMax * pi * pi * pi / (a*a*(pi*pi*pi-(T)32*sum)); // alpha = -dp/dz / mu
    T deltaP = - (alpha * nu);
    return deltaP;
}

T poiseuilleVelocity(plint iY, plint iZ, IncomprFlowParam<T> const& parameters, plint maxN)
{
    const T a = parameters.getNy()-1;
    const T b = parameters.getNz()-1;
    const T y = (T)iY - a / (T)2;
    const T z = (T)iZ - b / (T)2;
    const T alpha = - poiseuillePressure(parameters,maxN) / parameters.getLatticeNu();
    T sum = T();
    for (plint iN = 0; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum += (cos(twoNplusOne*pi*y/a)*cosh(twoNplusOne*pi*z/a)
            / ( pow(twoNplusOne,3)*cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }
    for (plint iN = 1; iN < maxN; iN += 2)
    {
        T twoNplusOne = (T)2*(T)iN+(T)1;
        sum -= (cos(twoNplusOne*pi*y/a)*cosh(twoNplusOne*pi*z/a)
            / ( pow(twoNplusOne,3)*cosh(twoNplusOne*pi*b/((T)2*a)) ));
    }
    sum *= ((T)4 * alpha * a *a /pow(pi,3));
    sum += (alpha / (T)2 * (y * y - a*a / (T)4));
    return sum;
}

template <typename T>
SquarePoiseuilleDensityAndVelocity<T>::SquarePoiseuilleDensityAndVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
: parameters(parameters_), maxN(maxN_)
{ }

template <typename T>
void SquarePoiseuilleDensityAndVelocity<T>::operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const {
        rho = (T)1;
        u[0] = poiseuilleVelocity(iY, iZ, parameters, maxN);
        u[1] = T();
        u[2] = T();
}


template <typename T>
SquarePoiseuilleVelocity<T>::SquarePoiseuilleVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
: parameters(parameters_), maxN(maxN_)
{ }

template <typename T>
void SquarePoiseuilleVelocity<T>::operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const  {
    u[0] = poiseuilleVelocity(iY, iZ, parameters, maxN);
    u[1] = T();
    u[2] = T();
}

/* ************* iniLattice ******************* */
void iniLattice( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                 IncomprFlowParam<T> const& parameters,
                 OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();

    Box3D top    = Box3D(0,    nx-1, 0, ny-1, 0,    0);
    Box3D bottom = Box3D(0,    nx-1, 0, ny-1, nz-1, nz-1);

    Box3D left   = Box3D(0, nx-1, 0,    0,    1, nz-2);
    Box3D right  = Box3D(0, nx-1, ny-1, ny-1, 1, nz-2);

    Box3D inlet  = Box3D(0,    0,    1,    ny-2, 1, nz-2);
    Box3D outlet = Box3D(nx-1, nx-1, 1,    ny-2, 1, nz-2);

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, inlet );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, outlet );

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, top );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, bottom );

    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, left );
    boundaryCondition.setVelocityConditionOnBlockBoundaries ( lattice, right );

    setBoundaryVelocity(lattice, inlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));
    setBoundaryVelocity(lattice, outlet, SquarePoiseuilleVelocity<T>(parameters, NMAX));

    setBoundaryVelocity(lattice, top, Array<T,3>(0.0,0.0,0.0));
    setBoundaryVelocity(lattice, bottom, Array<T,3>(0.0,0.0,0.0));
    setBoundaryVelocity(lattice, left, Array<T,3>(0.0,0.0,0.0));
    setBoundaryVelocity(lattice, right, Array<T,3>(0.0,0.0,0.0));

    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), SquarePoiseuilleDensityAndVelocity<T>(parameters, NMAX));

    setExternalVector( lattice, lattice.getBoundingBox(),
                       DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T,DESCRIPTOR<T>::d>(0.0,0.0,0.0));

    lattice.initialize();
}

/* ************* writeFicsionLogFile ******************* */
template<typename T>
void writeFicsionLogFile(IncomprFlowParam<T> const& parameters,
                  std::string const& title)
{
    std::string fullName = global::directories().getLogOutDir() + "plbLog.dat";
    plb_ofstream ofile(fullName.c_str());
    ofile << title << "\n\n";
    ofile << "Velocity in lattice units: u=" << parameters.getLatticeU() << "\n";
    ofile << "Reynolds number:           Re=" << parameters.getRe() << "\n";
    ofile << "Lattice resolution:        N=" << parameters.getResolution() << "\n";
    ofile << "Relaxation frequency:      omega=" << parameters.getOmega() << "\n";
    ofile << "Relaxation time:           tau=" << parameters.getTau() << "\n";
    ofile << "Grid spacing deltaX:       dx=" << parameters.getDeltaX() << "\n";
    ofile << "Time step deltaT:          dt=" << parameters.getDeltaT() << "\n";
    ofile << "Lattice  Nu:               nu=" << parameters.getLatticeNu() << "\n";
    ofile << "Physical Nu:               nu_p=" << parameters.getLatticeNu()*parameters.getDeltaX()*(parameters.getDeltaX()/parameters.getDeltaT()) << "\n";
    ofile << "Extent of the system:      lx=" << parameters.getLx() << "\n";
    ofile << "Extent of the system:      ly=" << parameters.getLy() << "\n";
    ofile << "Extent of the system:      lz=" << parameters.getLz() << "\n";
}


/* ************* Class GetTensorFieldFromExternalVectorFunctional3D ******************* */
template<typename T, template<typename U> class Descriptor, int nDim>
class GetTensorFieldFromExternalVectorFunctional3D : public BoxProcessingFunctional3D_LT<T,Descriptor, T, nDim> {
public:
    GetTensorFieldFromExternalVectorFunctional3D (
        int vectorStartsAt_ ) : vectorStartsAt(vectorStartsAt_)
    {
        PLB_ASSERT( vectorStartsAt+nDim <=
        Descriptor<T>::ExternalField::numScalars );
    }
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice, TensorField3D<T,nDim>& tensor) {
        Dot3D offset = computeRelativeDisplacement(lattice, tensor);
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            plint oX = iX + offset.x;
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint oY = iY + offset.y;
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    plint oZ = iZ + offset.z;
                    Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                    Array<T,nDim> externalVector;

                    for (plint iD=0; iD<nDim; ++iD) {
                        externalVector[iD] = *cell.getExternal(vectorStartsAt+iD);
                    }
                    tensor.get(oX,oY,oZ) = externalVector;
                }
            }
        }
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }
    virtual GetTensorFieldFromExternalVectorFunctional3D<T,Descriptor,nDim>* clone() const {
        return new GetTensorFieldFromExternalVectorFunctional3D<T,Descriptor,nDim>(*this);
    }

private:
    int vectorStartsAt;
};


/* ************* writeVTK ******************* */
template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();

    MultiTensorField3D<T,3> force(lattice);
    applyProcessingFunctional(new GetTensorFieldFromExternalVectorFunctional3D<T,DESCRIPTOR,3>(
        DESCRIPTOR<T>::ExternalField::forceBeginsAt), lattice.getBoundingBox(), lattice, force);

    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<T>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,T>(force, "force",  (T) dx/dt/dt);
    vtkOut.writeData<T>(*computeNorm(force, force.getBoundingBox()), "forceNorm",  dx/dt/dt);

//    ImageWriter<T> imageWriter("leeloo");
//    add(force, forceScalar, force.getBoundingBox());
//    imageWriter.writeScaledPpm(scalarField, createFileName("PPM", iter, 6));

//     vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
//     vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
}


#endif  // FICSIONINIT_HH
