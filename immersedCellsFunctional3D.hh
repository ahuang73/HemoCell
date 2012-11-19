#ifndef IMMERSED_CELLS_FUNCTIONAL_3D_HH
#define IMMERSED_CELLS_FUNCTIONAL_3D_HH

#include "immersedCellsFunctional3D.h"

#include "particles/particleField3D.h"
#include "dataProcessors/metaStuffFunctional3D.h"
#include "core/plbDebug.h"
#include "core/blockStatistics.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"


/* ******** countCellVolume *********************************** */
template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
void countCellVolume (TriangleBoundary3D<T> Cells,
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, std::vector<plint> cellIds,
                std::vector<T>& cellVolumes) //Perhaps add TAGS
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);

    ComputeCellVolumeParticlesFunctional3D<T,Descriptor> functional(Cells, cellIds);
    applyProcessingFunctional(functional, domain, particleArg);
    functional.getCellVolumeArray(cellVolumes, cellIds);
}

/* ******** ComputeCellVolumeParticlesFunctional3D *********************************** */
template<typename T, template<typename U> class Descriptor>
ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::ComputeCellVolumeParticlesFunctional3D(
        TriangleBoundary3D<T> const& triangleBoundary_, std::vector<plint> cellIds_)
    : triangleBoundary(triangleBoundary_)
{
    numberOfCells = 0;
    for (pluint var = 0; var < cellIds_.size(); ++var)
        if (cellIds_[var] > numberOfCells)
            numberOfCells = cellIds_[var];

    for (pluint i=0; i< (pluint) numberOfCells+1; ++i)
        volumeIds.push_back(this->getStatistics().subscribeSum());
//    pcout << "Done subscribing" << std::endl;
}

template<typename T, template<typename U> class Descriptor>
void ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    std::map<plint, T> tagsVolume;
    typename std::map<plint, T>::iterator tagsVolumeIterator;

//lmount: Find particles within this domain
    ParticleField3D<T,Descriptor>& particleField
        = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    
    std::map<plint, plint> tagsNr;
    typename std::map<plint, plint>::iterator tagsNrIterator;
    
    plint meshID = 1;
    triangleBoundary.pushSelect(0,meshID);
    TriangularSurfaceMesh<T> triangleMesh = triangleBoundary.getMesh();
    for (pluint iA = 0; iA < particles.size(); ++iA) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
        ImmersedCellParticle3D<T,Descriptor>* particle =
                dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
        
        plint iVertex = particle->getTag();
        std::vector<plint> neighbors = triangleMesh.getNeighborTriangleIds(iVertex);
        for (pluint iB = 0; iB < neighbors.size(); ++iB) {
            Array<T,3> v0 = triangleMesh.getVertex(neighbors[iB],0);
            Array<T,3> v1 = triangleMesh.getVertex(neighbors[iB],1);
            Array<T,3> v2 = triangleMesh.getVertex(neighbors[iB],2);
            
//             /* ************* Other Calculation ********* */
//             Array<T,3> areaTimesNormal = triangleMesh.computeTriangleNormal(neighbors[iB], true);
//             T triangleVolumeT6 = 2.0 * VectorTemplate<T,Descriptor>::scalarProduct(areaTimesNormal, ((v0+v1+v2)/3.0)) ;
//             /* ********************************************* */
            
            /* Calculating the volume contibution of a face based on the formula:
             * V[j] = 1.0/6.0 * (X3[j] cross X2[j])*X1[j]  */
            Array<T,3> tmp;
            crossProduct(v1, v2, tmp);
            T triangleVolumeT6 =  VectorTemplate<T,Descriptor>::scalarProduct(v0,tmp); // * (1.0/6.0)
            
            this->getStatistics().gatherSum(volumeIds[particle->get_cellId()], triangleVolumeT6/6.0/3.0); // every volume is evaluated 3 times
        }
        
    }    
    triangleBoundary.popSelect();
    
}


template<typename T, template<typename U> class Descriptor>
void ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::getCellVolumeArray(std::vector<T>& cellVolumes, std::vector<plint> cellIds) const {
    for (pluint i = 0; i < (pluint) volumeIds.size(); ++i)
        cellVolumes.push_back(this->getStatistics().getSum(volumeIds[i]));
}


template<typename T, template<typename U> class Descriptor>
ComputeCellVolumeParticlesFunctional3D<T,Descriptor>* ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::clone() const {
    return new ComputeCellVolumeParticlesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ComputeCellVolumeParticlesFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
}


/* ******** TranslateTaggedParticlesFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
TranslateTaggedParticlesFunctional3D<T,Descriptor>::
    TranslateTaggedParticlesFunctional3D(Array<T,3> const& translation_, plint tag_) : tag(tag_), translation(translation_)
{ }

template<typename T, template<typename U> class Descriptor>
void TranslateTaggedParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    std::vector<plint> cellIds;
    for (pluint iP = 0; iP < particles.size(); ++iP) {
        Particle3D<T,Descriptor>* nonTypedParticle = particles[iP];
        ImmersedCellParticle3D<T,Descriptor>* particle =
            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);

        if (particle->get_cellId() == tag) {
            particle->getPosition() = particle->getPosition() + translation;
        }
    }

//     pcout << cellIds.size() << std::endl;
//     particleField.removeParticles(domain,tag);
}

template<typename T, template<typename U> class Descriptor>
TranslateTaggedParticlesFunctional3D<T,Descriptor>* TranslateTaggedParticlesFunctional3D<T,Descriptor>::clone() const {
    return new TranslateTaggedParticlesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void TranslateTaggedParticlesFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}


//template<typename T, template<typename U> class Descriptor>
//CountCellVolumeFunctional3D<T,Descriptor>::CountCellVolumeFunctional3D(plint tag_)
//    : cellVolumeId(this->getStatistics().subscribeIntSum()), tag(tag_)
//{ }
//
//template<typename T, template<typename U> class Descriptor>
//void CountCellVolumeFunctional3D<T,Descriptor>::processGenericBlocks (
//        Box3D domain, std::vector<AtomicBlock3D*> blocks )
//{
//    PLB_PRECONDITION( blocks.size()==1 );
//    ParticleField3D<T,Descriptor>& particleField
//        = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
//    std::vector<Particle3D<T,Descriptor>*> particles;
//    particleField.findParticles(domain, particles);
//    T cellVolume = 0;
//    for (pluint iA = 0; iA < particles.size(); ++iA) {
//        Particle3D<T,Descriptor>* nonTypedParticle = particles[iA];
//        ImmersedCellParticle3D<T,Descriptor>* particle =
//            dynamic_cast<ImmersedCellParticle3D<T,Descriptor>*> (nonTypedParticle);
//
//        if (particle->get_cellId() == tag) {
////         if (particle->getTag() == tag) {
//        	cellVolume += 0.;
//        }
//    }
//    this->getStatistics().gatherIntSum(cellVolumeId, cellVolume);
//}
//
//template<typename T, template<typename U> class Descriptor>
//CountCellVolumeFunctional3D<T,Descriptor>* CountCellVolumeFunctional3D<T,Descriptor>::clone() const {
//    return new CountCellVolumeFunctional3D<T,Descriptor>(*this);
//}
//
//template<typename T, template<typename U> class Descriptor>
//void CountCellVolumeFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
//    modified[0] = modif::nothing;
//}
//
//template<typename T, template<typename U> class Descriptor>
//T CountCellVolumeFunctional3D<T,Descriptor>::getVollumeCell() const {
//    return this->getStatistics().getSum(cellVolumeId);
//}
//
//template< typename T, template<typename U> class Descriptor,
//          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
//T countCellVolume (
//                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, plint tag)
//{
//    std::vector<MultiBlock3D*> particleArg;
//    particleArg.push_back(&particles);
//
//    CountCellVolumeFunctional3D<T,Descriptor> functional(tag);
//    applyProcessingFunctional(functional, domain, particleArg);
//    return functional.getVollumeCell();
//}
//
//
//





#endif  // IMMERSED_CELLS_FUNCTIONAL_3D_H
