#ifndef IMMERSEDBOUNDARYMETHOD_3D_H
#define IMMERSEDBOUNDARYMETHOD_3D_H

#include "hemocell.h"
#include <vector>

namespace plb {

double phi2 (double x) ;

template<typename T>
T phi3 (T x) ;

template<typename T>
T phi4 (T x) ;

template<typename T>
T phi4c (T x) ;

template<typename T, template<typename U> class Descriptor>
void interpolationCoefficients (
        BlockLattice3D<T,Descriptor> const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights);

template<typename T, template<typename U> class Descriptor>
void interpolationCoefficientsPhi1 (
        BlockLattice3D<T,Descriptor> const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights);

void interpolationCoefficientsPhi2 (
        BlockLattice3D<double,DESCRIPTOR> const& block, SurfaceParticle3D * particle);

template<typename T, template<typename U> class Descriptor>
void interpolationCoefficientsPhi3 (
        BlockLattice3D<T,Descriptor> const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights);

void interpolationCoefficientsPhi4 (
        BlockLattice3D<double,DESCRIPTOR> const& block, Array<double,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<double>& weights);

template<typename T, template<typename U> class Descriptor>
void interpolationCoefficientsPhi4c (
        BlockLattice3D<T,Descriptor> const& block, Array<T,3> const& position,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights);

/*
 * In case one of the interpolating boundary nodes is a boundary,
 * the force is spread to the rest of the nodes
 * */
template<typename T, template<typename U> class Descriptor>
void curateInterpolationCoefficients (BlockLattice3D<T,Descriptor>& fluid,
        std::vector<Dot3D>& cellPos, std::vector<T>& weights) ;


}
#endif  // IMMERSEDBOUNDARYMETHOD_3D_H

