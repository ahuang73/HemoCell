#ifndef HEMOCELL_PARTICLE_FIELD_CPP
#define HEMOCELL_PARTICLE_FIELD_CPP

#include "hemoCellParticleField.h"
#include "hemocell.h"

/* *************** class HemoParticleField3D ********************** */

HemoCellParticleField::HemoCellParticleField(plint nx, plint ny, plint nz)
    : AtomicBlock3D(nx,ny,nz, &this->dataTransfer)
{ 
    dataTransfer.setBlock(*this);
    AddOutputMap(); 
}

HemoCellParticleField::HemoCellParticleField(HemoCellParticleField const& rhs)
    : AtomicBlock3D(rhs)
{
    dataTransfer.setBlock(*this);
    for (pluint i=0; i<rhs.particles.size(); ++i) {
        addParticle(this->getBoundingBox(),rhs.particles[i]->clone());
    }
    AddOutputMap();
}

HemoCellParticleField::~HemoCellParticleField()
{
  AtomicBlock3D::dataTransfer = new HemoCellParticleDataTransfer();
    for (pluint i=0; i<particles.size(); ++i) {
        delete particles[i];
    }
}

HemoCellParticleField& HemoCellParticleField::operator=(HemoCellParticleField const& rhs){
 HemoCellParticleField *copy = new HemoCellParticleField(rhs);
 delete this;
  return *copy;
}

HemoCellParticleField* HemoCellParticleField::clone() const
{
    return new HemoCellParticleField(*this);
}

void HemoCellParticleField::addParticle(Box3D domain, HemoCellParticle* particle) {
    //Box3D finalDomain;
    //plint x,y,z;
    HemoCellParticle * local_sparticle;
    Array<double,3> pos = particle->position;


    while (particles_per_type.size()<=particle->celltype) {
      particles_per_type.push_back(std::vector<HemoCellParticle*>());
    }

    if( this->isContainedABS(pos, this->getBoundingBox()) )
    {
      //check if we have particle already, if so, we must overwrite but not
      //forget to delete the old entry
      if ((!(particles_per_cell.find(particle->cellId) == 
           particles_per_cell.end())) && particles_per_cell[particle->cellId][particle->vertexId]) {
          local_sparticle =  particles_per_cell[particle->cellId][particle->vertexId];

          //If our particle is local, do not replace it, envelopes are less important
          if (isContainedABS(local_sparticle->position, getBoundingBox())) {
           
          } else {
            //We have the particle already, replace it
            *local_sparticle = *particle;
            local_sparticle->repoint_force_vectors();
          }
          //If the pointers are the same we shouldn't delete either
          if(!(particle==local_sparticle)) {
            delete particle;
            particle = local_sparticle;
          }
      } else {
        //new entry
        particles.push_back(particle);
          particles_per_type[particle->celltype].push_back(particle); //TODO, not accurate for already existing particles
          if(this->isContainedABS(pos, localDomain)) {
            lpc[particle->cellId] = true;
          }
          
        insert_ppc(particle);
      }

      particle->setTag(-1);
    }
    else {
        delete particle;
    }
}

void HemoCellParticleField::insert_ppc(HemoCellParticle* sparticle) {
  if (particles_per_cell.find(sparticle->cellId) == particles_per_cell.end()) {
    particles_per_cell[sparticle->cellId].resize((*cellFields)[sparticle->celltype]->numVertex);
    for (unsigned int i = 0; i < particles_per_cell[sparticle->cellId].size(); i++) {
      particles_per_cell[sparticle->cellId][i] = NULL;
    }
  }
  particles_per_cell[sparticle->cellId][sparticle->vertexId] = sparticle;

}


void HemoCellParticleField::removeParticles(plint tag) {
//Almost the same, but we save a lot of branching by making a seperate function
    vector<HemoCellParticle*> remainingParticles = particles;

    particles_per_type.clear();
    particles_per_cell.clear();
    lpc.clear();
    particles.clear();

    for (pluint i=0; i<remainingParticles.size(); ++i) {
       if (remainingParticles[i]->getTag() == tag) {
         delete remainingParticles[i];
       }
       else {
           addParticle(this->getBoundingBox(), remainingParticles[i]);
       }
    }
}

void HemoCellParticleField::removeParticles(Box3D domain, plint tag) {
//Almost the same, but we save a lot of branching by making a seperate function
    std::vector<HemoCellParticle*> remainingParticles = particles;
    Box3D finalDomain;
    Array<double,3> pos; 

    particles_per_type.clear();
    particles_per_cell.clear();
    lpc.clear();
    particles.clear();
    
    for (pluint i=0; i<remainingParticles.size(); ++i) {
       pos = remainingParticles[i]->position;
       intersect(domain, this->getBoundingBox(), finalDomain);
       if (this->isContainedABS(pos,finalDomain) && remainingParticles[i]->getTag() == tag) {
         delete remainingParticles[i];
       }
       else {
           addParticle(this->getBoundingBox(), remainingParticles[i]);
       }
    }
}

void HemoCellParticleField::removeParticles(Box3D domain) {
//Almost the same, but we save a lot of branching by making a seperate function
    std::vector<HemoCellParticle*> remainingParticles = particles;
    Box3D finalDomain;
    Array<double,3> pos; 

    particles_per_type.clear();
    particles_per_cell.clear();
    lpc.clear();
    particles.clear();

    for (pluint i=0; i<remainingParticles.size(); ++i) {
       pos = remainingParticles[i]->position;
       intersect(domain, this->getBoundingBox(), finalDomain);
       if (this->isContainedABS(pos,finalDomain)) {
         delete remainingParticles[i];
       }
       else {
           addParticle(this->getBoundingBox(), remainingParticles[i]);
       }
    }
}

void HemoCellParticleField::findParticles (
        Box3D domain, std::vector<HemoCellParticle*>& found ) 
{
    found.clear();
    PLB_ASSERT( contained(domain, this->getBoundingBox()) );
    Array<double,3> pos; 
    for (pluint i=0; i<particles.size(); ++i) {
        pos = particles[i]->position;
        if (this->isContainedABS(pos,domain)) {
            found.push_back(particles[i]);
        }
    }
}
void HemoCellParticleField::findParticles (
        Box3D domain, std::vector<HemoCellParticle*>& found, pluint type) const
{
    
    found.clear();
    PLB_ASSERT( contained(domain, this->getBoundingBox()) );
    Array<double,3> pos; 
    if (!(particles_per_type.size() > type)) {return;} else {
    for (pluint i=0; i<particles_per_type[type].size(); ++i) {
        pos = particles_per_type[type][i]->position;
        if (this->isContainedABS(pos,domain)) {
            found.push_back(particles_per_type[type][i]);
        }
    }
    }
    
}

bool HemoCellParticleField::isContainedABS(Array<double,3> pos, Box3D box) const {
		Dot3D const& location = this->getLocation();
    double x = pos[0]-location.x;
    double y = pos[1]-location.y;
    double z = pos[2]-location.z;
    if (box.z1 < -10000) {
      exit(0);
    } 
    if (location.x < -10000) {
      exit(0);
    } 

    return (x > box.x0-0.5) && (x <= box.x1+0.5) &&
           (y > box.y0-0.5) && (y <= box.y1+0.5) &&
           (z > box.z0-0.5) && (z <= box.z1+0.5);

}

inline plint HemoCellParticleField::nearestCell(double const pos) const {
  return int(pos + 0.5);
}

inline void HemoCellParticleField::computeGridPosition (
            Array<double,3> const& position,
                    plint* iX, plint* iY, plint* iZ ) const
{
      Dot3D const& location = this->getLocation();
      *iX = nearestCell(position[0]) - location.x;
      *iY = nearestCell(position[1]) - location.y;
      *iZ = nearestCell(position[2]) - location.z;
}

void HemoCellParticleField::computeGridPosition (
            Array<double,3> const& position,
                    plint& iX, plint& iY, plint& iZ ) const
{
      Dot3D const& location = this->getLocation();
      iX = nearestCell(position[0]) - location.x;
      iY = nearestCell(position[1]) - location.y;
      iZ = nearestCell(position[2]) - location.z;
}

int HemoCellParticleField::deleteIncompleteCells(pluint ctype, bool twice) {
  //Function must be called twice since addParticle can remove a particle
  //unintentionally, for now, catch it here; TODO, this can be done better
  int deleted = 0;

  //Warning, TODO, high complexity, should be rewritten 
  //For now abuse tagging and the remove function
  for ( const auto &lpc_it : particles_per_cell ) {
    int cellid = lpc_it.first;
    bool broken = false;
    for (pluint i = 0; i < particles_per_cell[cellid].size() ; i++) {
      if (!particles_per_cell[cellid][i]) {
        broken = true;
        break;
      }
    }
    if (!broken) {continue;}

    //actually add to tobedeleted list
    for (pluint i = 0; i < particles_per_cell[cellid].size() ; i++) {
      if (particles_per_cell[cellid][i] == NULL) {continue;}
      if (particles_per_cell[cellid][i]->celltype != ctype) {break;} //certainly a entry, therefore we check here if it is the right type, if not, exit

      particles_per_cell[cellid][i]->setTag(1);
      deleted++;
    }
  } 

  //We have our list, now abuse the removeall function
  removeParticles(1);

  return deleted; 
}

int HemoCellParticleField::deleteIncompleteCells(bool twice) {
  //Function must be called twice since addParticle can remove a particle
  //unintentionally, for now, catch it here; TODO, this can be done better
  int deleted = 0;

  //Warning, TODO, high complexity, should be rewritten 
  //For now abuse tagging and the remove function
  for ( const auto &lpc_it : particles_per_cell ) {
    int cellid = lpc_it.first;
    bool broken = false;
    for (pluint i = 0; i < particles_per_cell[cellid].size() ; i++) {
      if (!particles_per_cell[cellid][i]) {
        broken = true;
        break;
      }
    }
    if (!broken) {continue;}

    //actually add to tobedeleted list
    for (pluint i = 0; i < particles_per_cell[cellid].size() ; i++) {
      if (particles_per_cell[cellid][i] == NULL) {continue;}
      particles_per_cell[cellid][i]->setTag(1);
      deleted++;
    }
  } 

  //We have our list, now abuse the removeall function

  removeParticles(1);

  return deleted; 
}

void HemoCellParticleField::setlocalDomain(Box3D & localDomain_) {
  localDomain = localDomain_;
  localDomain.x0 -= this->getLocation().x;
  localDomain.x1 -= this->getLocation().x;
  localDomain.y0 -= this->getLocation().y;
  localDomain.y1 -= this->getLocation().y;
  localDomain.z0 -= this->getLocation().z;
  localDomain.z1 -= this->getLocation().z;
}


void HemoCellParticleField::advanceParticles() {
  for(HemoCellParticle* particle:particles){
    particle->advance();
  }
  removeParticles(1);
}

void HemoCellParticleField::separateForceVectors() {
  //Also save the total force, therfore recalculate in advance
  applyConstitutiveModel();

  for (HemoCellParticle* sparticle : particles) {
    //Save Total Force
    sparticle->force_total = sparticle->force + sparticle->force_repulsion;

    //Just repoint all possible outputs for now //TODO only repoint the ones we
    //want
    ////TODO this can leak if particle is deleted between seperate and unify,
    //rewrite to reference
    sparticle->force_volume = new Array<double,3>(0.0,0.0,0.0);
    sparticle->force_link = new Array<double,3>(0.0,0.0,0.0);
    sparticle->force_area = new Array<double,3>(0.0,0.0,0.0);
    sparticle->force_bending = new Array<double,3>(0.0,0.0,0.0);
  }
}

void HemoCellParticleField::unifyForceVectors() {
  for (HemoCellParticle* sparticle : particles) {
    //Just repoint all possible outputs for now //TODO only repoint the ones we
    //want
    delete sparticle->force_volume;
    delete sparticle->force_link;
    delete sparticle->force_area;
    delete sparticle->force_bending;
    sparticle->repoint_force_vectors();
  }
}

void HemoCellParticleField::applyConstitutiveModel(bool forced) {
  deleteIncompleteCells();
  
  for (pluint ctype = 0; ctype < (*cellFields).size(); ctype++) {
    if ((*cellFields).hemocell.iter % (*cellFields)[ctype]->timescale == 0 || forced) {
      vector<HemoCellParticle*> found;
      findParticles(getBoundingBox(),found,ctype);
      for (HemoCellParticle* particle : found) {
        particle->force = {0.,0.,0.};
      }
      (*cellFields)[ctype]->mechanics->ParticleMechanics(particles_per_cell,lpc,ctype);
    }
  }
  
}

void HemoCellParticleField::applyRepulsionForce(bool forced) {
  const double r_const = cellFields->repulsionConstant;
  const double r_cutoff = cellFields->repulsionCutoff;
  
  vector<HemoCellParticle *> found;
  Box3D localDomainplusOne = localDomain;
  localDomainplusOne.enlarge(1);
  findParticles(localDomainplusOne,found);
  
  if (cellFields->hemocell.iter % cellFields->repulsionTimescale == 0 || forced) {  
    for (HemoCellParticle * particle : found) {
      particle->force_repulsion = {0.,0.,0.};
    }
    
    if (found.size() > 1) {
      for (unsigned int i = 0 ; i < found.size() - 1 ; i ++) {
        for (unsigned int j = i + 1 ; j < found.size() ; j++ ) {
          if (found[i]->cellId == found[j]->cellId) { continue; } //No incell repulsion
          const Array<double,3> dv = found[i]->position - found[j]->position;
          const double distance = sqrt(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2]);
          if (distance < r_cutoff) {
            const Array<double, 3> rfm = r_const * (1/(distance/r_cutoff))  * (dv/distance);
            found[i]->force_repulsion = found[i]->force_repulsion + rfm;
            found[j]->force_repulsion = found[j]->force_repulsion - rfm;
          } 
        }
      }
    }
  }
}

void HemoCellParticleField::spreadParticleForce(Box3D domain) {
  vector<HemoCellParticle*> localParticles;
  findParticles(domain,localParticles);
  HemoCellParticle * sparticle;

  for (pluint i = 0; i < localParticles.size(); i++ ) {
    sparticle = localParticles[i];
    //Clever trick to allow for different kernels for different particle types.
    (*cellFields)[sparticle->celltype]->kernelMethod(*atomicLattice,sparticle);

    //Directly change the force on a node , Palabos developers hate this one
    //quick non-functional trick.
    for (pluint j = 0; j < sparticle->kernelLocations.size(); j++) {
      //Yay for direct access
      sparticle->kernelLocations[j]->external.data[0] += ((sparticle->force_repulsion[0] + sparticle->force[0]) * sparticle->kernelWeights[j]);
      sparticle->kernelLocations[j]->external.data[1] += ((sparticle->force_repulsion[1] + sparticle->force[1]) * sparticle->kernelWeights[j]);
      sparticle->kernelLocations[j]->external.data[2] += ((sparticle->force_repulsion[2] + sparticle->force[2]) * sparticle->kernelWeights[j]);
    }

  }
}

void HemoCellParticleField::interpolateFluidVelocity(Box3D domain) {
  vector<HemoCellParticle*> localParticles;
  findParticles(domain,localParticles);
  //TODO, remove casting
  HemoCellParticle * sparticle;
  //Prealloc is nice
  Array<double,3> velocity;
  Array<double,3> velocity_comp;


  for (pluint i = 0; i < localParticles.size(); i++ ) {
    sparticle = localParticles[i];
    //Clever trick to allow for different kernels for different particle types.
    (*cellFields)[sparticle->celltype]->kernelMethod(*atomicLattice,sparticle);

    //We have the kernels, now calculate the velocity of the particles.
    //Palabos developers, sorry for not using a functional...
    velocity = {0.0,0.0,0.0};
    for (pluint j = 0; j < sparticle->kernelLocations.size(); j++) {
      //Yay for direct access
      sparticle->kernelLocations[j]->computeVelocity(velocity_comp);
      velocity += velocity_comp * sparticle->kernelWeights[j];
    }
    sparticle->v = velocity;
  }

}

HemoCellParticleDataTransfer& HemoCellParticleField::getDataTransfer() {
    return dataTransfer;
}
HemoCellParticleDataTransfer const& HemoCellParticleField::getDataTransfer() const {
    return dataTransfer;
}

std::string HemoCellParticleField::getBlockName() {
    return std::string("HemoParticleField3D");
}

HemoCellFields* HemoCellParticleField::cellFields=0;

#endif
