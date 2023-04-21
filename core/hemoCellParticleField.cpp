/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab
in the University of Amsterdam. Any questions or remarks regarding this library
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "hemoCellParticleField.h"
#include "hemocell.h"
#include "octree.h"
#include "mollerTrumbore.h"
#include "bindingField.h"
#include "interiorViscosity.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#include <Eigen3/Eigenvalues>
#pragma GCC diagnostic pop
#include "cellInfo.h"
namespace hemo
{
  /* *************** class HemoParticleField3D ********************** */

  HemoCellParticleField::HemoCellParticleField(plint nx, plint ny, plint nz)
      : AtomicBlock3D(nx, ny, nz, 0), particleDataTransfer(*(new HemoCellParticleDataTransfer()))
  {
    boundingBox = Box3D(0, this->getNx() - 1, 0, this->getNy() - 1, 0, this->getNz() - 1);
    dataTransfer = &particleDataTransfer;
    particleDataTransfer.setBlock(*this);
    AddOutputMap();
  }

  HemoCellParticleField::HemoCellParticleField(HemoCellParticleField const &rhs)
      : AtomicBlock3D(rhs), particleDataTransfer(*(new HemoCellParticleDataTransfer()))
  {
    boundingBox = Box3D(0, this->getNx() - 1, 0, this->getNy() - 1, 0, this->getNz() - 1);
    dataTransfer = &particleDataTransfer;
    particleDataTransfer.setBlock(*this);
    for (const HemoCellParticle &particle : rhs.particles)
    {
      addParticle(particle.sv);
    }
    ppc_up_to_date = false;
    lpc_up_to_date = false;
    ppt_up_to_date = false;
    pg_up_to_date = false;
    AddOutputMap();
  }

  HemoCellParticleField::~HemoCellParticleField()
  {
    // AtomicBlock3D::dataTransfer = new HemoCellParticleDataTransfer();
    if (particle_grid)
    {
      delete[] particle_grid;
      particle_grid = 0;
    }
    if (particle_grid_size)
    {
      delete[] particle_grid_size;
      particle_grid_size = 0;
    }

    // Sanitize for MultiBlockLattice destructor (releasememory). It can't handle releasing non-background dynamics that are not singular
    if (global.enableInteriorViscosity)
    {
      for (const Dot3D &internalPoint : internalPoints)
      {
        atomicLattice->get(internalPoint.x, internalPoint.y, internalPoint.z);
        atomicLattice->get(internalPoint.x, internalPoint.y, internalPoint.z).attributeDynamics(&atomicLattice->getBackgroundDynamics());
      }
      InteriorViscosityHelper::get(*cellFields).empty(*this);
    }
  }

  HemoCellParticleField &HemoCellParticleField::operator=(HemoCellParticleField const &rhs)
  {
    HemoCellParticleField *copy = new HemoCellParticleField(rhs);
    delete this;
    return *copy;
  }

  HemoCellParticleField *HemoCellParticleField::clone() const
  {
    return new HemoCellParticleField(*this);
  }

  const vector<vector<unsigned int>> &HemoCellParticleField::get_particles_per_type()
  {
    if (!ppt_up_to_date)
    {
      update_ppt();
    }
    return _particles_per_type;
  }
  const map<int, vector<int>> &HemoCellParticleField::get_particles_per_cell()
  {
    if (!ppc_up_to_date)
    {
      update_ppc();
    }
    return _particles_per_cell;
  }

  const map<int, bool> &HemoCellParticleField::get_lpc()
  {
    if (!lpc_up_to_date)
    {
      update_lpc();
    }
    return _lpc;
  }
  void HemoCellParticleField::update_lpc()
  {
    _lpc.clear();
    for (const HemoCellParticle &particle : particles)
    {
      if (isContainedABS(particle.sv.position, localDomain))
      {
        _lpc[particle.sv.cellId] = true;
      }
    }
    lpc_up_to_date = true;
  }
  void HemoCellParticleField::update_ppt()
  {
    _particles_per_type.clear();
    _particles_per_type.resize(cellFields->size());

    for (unsigned int i = 0; i < particles.size(); i++)
    {
      _particles_per_type[particles[i].sv.celltype].push_back(i);
    }
    ppt_up_to_date = true;
  }
  void HemoCellParticleField::update_ppc()
  {
    _particles_per_cell.clear();

    for (unsigned int i = 0; i < particles.size(); i++)
    {
      insert_ppc(&particles[i], i);
    }
    ppc_up_to_date = true;
  }

  void HemoCellParticleField::update_pg()
  {
    // Check if map exists, otherwise create it
    if (!particle_grid)
    {
      particle_grid = new hemo::Array<unsigned int, 10>[this->atomicLattice->getNx() * this->atomicLattice->getNy() * this->atomicLattice->getNz()];
    }
    if (!particle_grid_size)
    {
      particle_grid_size = new unsigned int[this->atomicLattice->getNx() * this->atomicLattice->getNy() * this->atomicLattice->getNz()];
    }
    if (!this->atomicLattice)
    {
      return;
    }

    memset(particle_grid_size, 0, sizeof(unsigned int) * this->atomicLattice->getNx() * this->atomicLattice->getNy() * this->atomicLattice->getNz());
    Dot3D const &location = this->atomicLattice->getLocation();
    hemo::Array<T, 3> *pos;

    for (unsigned int i = 0; i < particles.size(); i++)
    {
      pos = &particles[i].sv.position;
      int x = pos->operator[](0) - location.x + 0.5;
      int y = pos->operator[](1) - location.y + 0.5;
      int z = pos->operator[](2) - location.z + 0.5;
      if ((x >= 0) && (x < this->atomicLattice->getNx()) &&
          (y >= 0) && (y < this->atomicLattice->getNy()) &&
          (z >= 0) && (z < this->atomicLattice->getNz()))
      {
        unsigned int index = grid_index(x, y, z);
        particle_grid[index][particle_grid_size[index]] = i;
        particle_grid_size[index]++;
      }
    }
    pg_up_to_date = true;
  }

  void HemoCellParticleField::addParticle(HemoCellParticle *particle)
  {
    addParticle(particle->sv);
  }
  void HemoCellParticleField::addParticle(const HemoCellParticle::serializeValues_t &sv)
  {
    HemoCellParticle *local_sparticle, *particle;
    const hemo::Array<T, 3> &pos = sv.position;
    const map<int, vector<int>> &particles_per_cell = get_particles_per_cell();

    if (this->isContainedABS(pos, this->getBoundingBox()))
    {
      // check if we have particle already, if so, we must overwrite but not
      // forget to delete the old entry
      if ((!(particles_per_cell.find(sv.cellId) == particles_per_cell.end())))
      {
        if (particles_per_cell.at(sv.cellId)[sv.vertexId] != -1)
        {
          local_sparticle = &particles[particles_per_cell.at(sv.cellId)[sv.vertexId]];

          // If our particle is local, do not replace it, envelopes are less important
          if (isContainedABS(local_sparticle->sv.position, localDomain))
          {
            return;
          }
          else
          {
            // We have the particle already, replace it
            local_sparticle->sv = sv;
            particle = local_sparticle;
            particle->setTag(-1);

            // Invalidate lpc hemo::Array
            lpc_up_to_date = false;
            pg_up_to_date = false;
          }
        }
        else
        {
          goto outer_else;
        }
      }
      else
      {
      outer_else:
        // new entry
        particles.emplace_back(sv);
        particle = &particles.back();

        // invalidate ppt
        ppt_up_to_date = false;
        if (this->isContainedABS(pos, localDomain))
        {
          _lpc[particle->sv.cellId] = true;
        }
        if (ppc_up_to_date)
        { // Otherwise its rebuild anyway
          insert_ppc(particle, particles.size() - 1);
        }

        if (pg_up_to_date)
        {
          Dot3D const &location = this->atomicLattice->getLocation();
          hemo::Array<T, 3> &pos = particle->sv.position;
          int x = pos[0] - location.x + 0.5;
          int y = pos[1] - location.y + 0.5;
          int z = pos[2] - location.z + 0.5;
          if ((x >= 0) && (x <= this->atomicLattice->getNx()) &&
              (y >= 0) && (y <= this->atomicLattice->getNy()) &&
              (z >= 0) && (z <= this->atomicLattice->getNz()))
          {
            unsigned int index = grid_index(x, y, z);
            particle_grid[index][particle_grid_size[index]] = particles.size() - 1;
            particle_grid_size[index]++;
          }
        }
      }
    }
  }

  void HemoCellParticleField::addParticlePreinlet(const HemoCellParticle::serializeValues_t &sv)
  {
    HemoCellParticle *local_sparticle, *particle;
    const hemo::Array<T, 3> &pos = sv.position;
    const map<int, vector<int>> &particles_per_cell = get_particles_per_cell();

    if (this->isContainedABS(pos, this->getBoundingBox()))
    {
      // check if we have particle already, if so, we must overwrite but not
      // forget to delete the old entry
      if ((!(particles_per_cell.find(sv.cellId) == particles_per_cell.end())))
      {
        if (particles_per_cell.at(sv.cellId)[sv.vertexId] != -1)
        {
          return;
        }
        else
        {
          goto outer_else;
        }
      }
      else
      {
      outer_else:
        // new entry
        particles.emplace_back(sv);
        particle = &particles.back();

        // invalidate ppt
        ppt_up_to_date = false;
        if (this->isContainedABS(pos, localDomain))
        {
          _lpc[particle->sv.cellId] = true;
        }
        if (ppc_up_to_date)
        { // Otherwise its rebuild anyway
          insert_ppc(particle, particles.size() - 1);
        }

        if (pg_up_to_date)
        {
          Dot3D const &location = this->atomicLattice->getLocation();
          hemo::Array<T, 3> &pos = particle->sv.position;
          int x = pos[0] - location.x + 0.5;
          int y = pos[1] - location.y + 0.5;
          int z = pos[2] - location.z + 0.5;
          if ((x >= 0) && (x <= this->atomicLattice->getNx()) &&
              (y >= 0) && (y <= this->atomicLattice->getNy()) &&
              (z >= 0) && (z <= this->atomicLattice->getNz()))
          {
            unsigned int index = grid_index(x, y, z);
            particle_grid[index][particle_grid_size[index]] = particles.size() - 1;
            particle_grid_size[index]++;
          }
        }
      }
    }
  }

  void inline HemoCellParticleField::insert_ppc(HemoCellParticle *sparticle, unsigned int index)
  {
    if (_particles_per_cell.find(sparticle->sv.cellId) == _particles_per_cell.end())
    {
      _particles_per_cell[sparticle->sv.cellId].resize((*cellFields)[sparticle->sv.celltype]->numVertex, -1);
    }
    _particles_per_cell.at(sparticle->sv.cellId)[sparticle->sv.vertexId] = index;
  }
  void inline HemoCellParticleField::insert_preinlet_ppc(HemoCellParticle *sparticle, unsigned int index)
  {
    if (_preinlet_particles_per_cell.find(sparticle->sv.cellId) == _preinlet_particles_per_cell.end())
    {
      _preinlet_particles_per_cell[sparticle->sv.cellId].resize((*cellFields)[sparticle->sv.celltype]->numVertex);
      for (unsigned int i = 0; i < _preinlet_particles_per_cell[sparticle->sv.cellId].size(); i++)
      {
        _preinlet_particles_per_cell[sparticle->sv.cellId][i] = -1;
      }
    }
    _preinlet_particles_per_cell.at(sparticle->sv.cellId)[sparticle->sv.vertexId] = index;
  }

  void HemoCellParticleField::removeParticles(plint tag)
  {
    // Almost the same, but we save a lot of branching by making a seperate function

    const unsigned int old_size = particles.size();
    for (unsigned int i = 0; i < particles.size(); i++)
    {
      if (particles[i].getTag() == tag)
      {
        particles[i] = particles.back();
        particles.pop_back();
        i--;
      }
    }
    if (particles.size() != old_size)
    {
      lpc_up_to_date = false;
      ppt_up_to_date = false;
      ppc_up_to_date = false;
      pg_up_to_date = false;
    }
  }

  void HemoCellParticleField::removeParticles(Box3D domain, plint tag)
  {
    // Almost the same, but we save a lot of branching by making a seperate function
    Box3D finalDomain;

    intersect(domain, this->getBoundingBox(), finalDomain);

    const unsigned int old_size = particles.size();
    for (unsigned int i = 0; i < particles.size(); i++)
    {
      if (particles[i].getTag() == tag && this->isContainedABS(particles[i].sv.position, finalDomain))
      {
        particles[i] = particles.back();
        particles.pop_back();
        i--;
      }
    }
    if (particles.size() != old_size)
    {
      lpc_up_to_date = false;
      ppt_up_to_date = false;
      ppc_up_to_date = false;
      pg_up_to_date = false;
    }
  }

  void HemoCellParticleField::removeParticles(Box3D domain)
  {
    // Almost the same, but we save a lot of branching by making a seperate function

    Box3D finalDomain;

    intersect(domain, this->getBoundingBox(), finalDomain);

    const unsigned int old_size = particles.size();
    for (unsigned int i = 0; i < particles.size(); i++)
    {
      if (this->isContainedABS(particles[i].sv.position, finalDomain))
      {
        particles[i] = particles.back();
        particles.pop_back();
        i--;
      }
    }
    if (particles.size() != old_size)
    {
      lpc_up_to_date = false;
      ppt_up_to_date = false;
      ppc_up_to_date = false;
      pg_up_to_date = false;
    }
  }

  // remove everything outside this domain
  void HemoCellParticleField::removeParticles_inverse(Box3D domain)
  {
    // Almost the same, but we save a lot of branching by making a seperate function

    Box3D finalDomain;

    intersect(domain, this->getBoundingBox(), finalDomain);

    const unsigned int old_size = particles.size();
    for (unsigned int i = 0; i < particles.size(); i++)
    {
      if (!this->isContainedABS(particles[i].sv.position, finalDomain))
      {
        particles[i] = particles.back();
        particles.pop_back();
        i--;
      }
    }
    if (particles.size() != old_size)
    {
      lpc_up_to_date = false;
      ppt_up_to_date = false;
      ppc_up_to_date = false;
      pg_up_to_date = false;
    }
  }

  void HemoCellParticleField::syncEnvelopes()
  {
    removeParticles_inverse(localDomain);
  }

  void HemoCellParticleField::findParticles(
      Box3D domain, std::vector<HemoCellParticle *> &found)
  {
    found.clear();
    PLB_ASSERT(contained(domain, this->getBoundingBox()));
    for (HemoCellParticle &particle : particles)
    {
      if (this->isContainedABS(particle.sv.position, domain))
      {
        found.push_back(&particle);
      }
    }
  }

  void HemoCellParticleField::findParticles(
      Box3D domain, std::vector<const HemoCellParticle *> &found) const
  {
    found.clear();
    PLB_ASSERT(contained(domain, this->getBoundingBox()));
    for (const HemoCellParticle &particle : particles)
    {
      if (this->isContainedABS(particle.sv.position, domain))
      {
        found.push_back(&particle);
      }
    }
  }
  void HemoCellParticleField::findParticles(
      Box3D domain, std::vector<HemoCellParticle *> &found, pluint type)
  {

    found.clear();
    PLB_ASSERT(contained(domain, this->getBoundingBox()));
    // hemo::Array<T,3> pos;
    const vector<vector<unsigned int>> &particles_per_type = get_particles_per_type();
    if (!(particles_per_type.size() > type))
    {
      return;
    }
    else
    {
      for (const unsigned int i : particles_per_type[type])
      {
        if (this->isContainedABS(particles[i].sv.position, domain))
        {
          found.push_back(&(particles[i]));
        }
      }
    }
  }

  inline plint HemoCellParticleField::nearestCell(T const pos) const
  {
    return int(pos + 0.5);
  }

  inline void HemoCellParticleField::computeGridPosition(
      hemo::Array<T, 3> const &position,
      plint *iX, plint *iY, plint *iZ) const
  {
    Dot3D const &location = this->getLocation();
    *iX = nearestCell(position[0]) - location.x;
    *iY = nearestCell(position[1]) - location.y;
    *iZ = nearestCell(position[2]) - location.z;
  }

  void HemoCellParticleField::computeGridPosition(
      hemo::Array<T, 3> const &position,
      plint &iX, plint &iY, plint &iZ) const
  {
    Dot3D const &location = this->getLocation();
    iX = nearestCell(position[0]) - location.x;
    iY = nearestCell(position[1]) - location.y;
    iZ = nearestCell(position[2]) - location.z;
  }

  void HemoCellParticleField::issueWarning(HemoCellParticle &p)
  {
    cout << "(HemoCell) (Delete Cells) WARNING! Particle deleted from local domain. This means the whole cell will be deleted!" << endl;
    cout << "\t Particle ID:" << p.sv.cellId << endl;
    cout << "\t Position: " << p.sv.position[0] << ", " << p.sv.position[1] << ", " << p.sv.position[2] << "; vel.: " << p.sv.v[0] << ", " << p.sv.v[1] << ", " << p.sv.v[2] << "; force: " << p.sv.force[0] << ", " << p.sv.force[1] << ", " << p.sv.force[2] << endl;
  }

  int HemoCellParticleField::deleteIncompleteCells(pluint ctype, bool verbose)
  {
    int deleted = 0;

    const map<int, vector<int>> &particles_per_cell = get_particles_per_cell();
    // Warning, TODO, high complexity, should be rewritten
    // For now abuse tagging and the remove function
    for (const auto &lpc_it : particles_per_cell)
    {
      int cellid = lpc_it.first;
      bool broken = false;
      for (pluint i = 0; i < particles_per_cell.at(cellid).size(); i++)
      {
        if (particles_per_cell.at(cellid)[i] == -1)
        {
          broken = true;
          break;
        }
      }
      if (!broken)
      {
        continue;
      }

      bool warningIssued = false;
      for (pluint i = 0; i < particles_per_cell.at(cellid).size(); i++)
      {
        if (particles_per_cell.at(cellid)[i] == -1)
        {
          continue;
        }

        // issue warning
        if (verbose)
        {
          if (!warningIssued)
          {
            if (isContainedABS(particles[particles_per_cell.at(cellid)[i]].sv.position, localDomain))
            {
              issueWarning(particles[particles_per_cell.at(cellid)[i]]);
              warningIssued = true;
            }
          }
        }

        // actually add to tobedeleted list
        particles[particles_per_cell.at(cellid)[i]].setTag(1);
        deleted++;
      }
    }

    // We have our list, now abuse the removeall function
    removeParticles(1);

    return deleted;
  }

  int HemoCellParticleField::deleteIncompleteCells(const bool verbose)
  {
    int deleted = 0;
    const map<int, vector<int>> &particles_per_cell = get_particles_per_cell();
    // Warning, TODO, high complexity, should be rewritten
    // For now abuse tagging and the remove function
    for (const auto &lpc_it : particles_per_cell)
    {
      int cellid = lpc_it.first;
      bool broken = false;
      for (pluint i = 0; i < particles_per_cell.at(cellid).size(); i++)
      {
        if (particles_per_cell.at(cellid)[i] == -1)
        {
          broken = true;
          break;
        }
      }
      if (!broken)
      {
        continue;
      }

      bool warningIssued = false;
      for (pluint i = 0; i < particles_per_cell.at(cellid).size(); i++)
      {
        if (particles_per_cell.at(cellid)[i] == -1)
        {
          continue;
        }

        // issue warning
        if (verbose)
        {
          if (!warningIssued)
          {
            if (isContainedABS(particles[particles_per_cell.at(cellid)[i]].sv.position, localDomain))
            {
              issueWarning(particles[particles_per_cell.at(cellid)[i]]);
              warningIssued = true;
            }
          }
        }

        // actually add to tobedeleted list
        particles[particles_per_cell.at(cellid)[i]].setTag(1);
        deleted++;
      }
    }

    // We have our list, now abuse the removeall function

    removeParticles(1);

    return deleted;
  }

  void HemoCellParticleField::setlocalDomain(Box3D &localDomain_)
  {
    localDomain = localDomain_;
    localDomain.x0 -= this->getLocation().x;
    localDomain.x1 -= this->getLocation().x;
    localDomain.y0 -= this->getLocation().y;
    localDomain.y1 -= this->getLocation().y;
    localDomain.z0 -= this->getLocation().z;
    localDomain.z1 -= this->getLocation().z;
  }

  void HemoCellParticleField::advanceParticles()
  {
    for (HemoCellParticle &particle : particles)
    {
      particle.advance();
      // By lack of better place, check if it is on a boundary, if so, delete it
      plb::Box3D const box = atomicLattice->getBoundingBox();
      plb::Dot3D const &location = atomicLattice->getLocation();
      plint x = (particle.sv.position[0] - location.x) + 0.5;
      plint y = (particle.sv.position[1] - location.y) + 0.5;
      plint z = (particle.sv.position[2] - location.z) + 0.5;

      if ((x >= box.x0) && (x <= box.x1) &&
          (y >= box.y0) && (y <= box.y1) &&
          (z >= box.z0) && (z <= box.z1))
      {
        if (atomicLattice->get(x, y, z).getDynamics().isBoundary())
        {
          particle.tag = 1;
        }
      }
    }
    removeParticles(1);

    lpc_up_to_date = false;
    pg_up_to_date = false;
  }

  void HemoCellParticleField::separateForceVectors()
  {
    // Also save the total force, therfore recalculate in advance
    applyConstitutiveModel();

    for (HemoCellParticle &sparticle : particles)
    {
      // Save Total Force
      sparticle.force_total = sparticle.sv.force + sparticle.sv.force_repulsion + sparticle.sv.force_adhesion;

      // Just repoint all possible outputs for now //TODO only repoint the ones we
      // want

      sparticle.force_volume = new hemo::Array<T, 3>({0.0, 0.0, 0.0});
      allocated_for_output.push_back(sparticle.force_volume);
      sparticle.force_link = new hemo::Array<T, 3>({0.0, 0.0, 0.0});
      allocated_for_output.push_back(sparticle.force_link);
      sparticle.force_area = new hemo::Array<T, 3>({0.0, 0.0, 0.0});
      allocated_for_output.push_back(sparticle.force_area);
      sparticle.force_bending = new hemo::Array<T, 3>({0.0, 0.0, 0.0});
      allocated_for_output.push_back(sparticle.force_bending);
      sparticle.force_visc = new hemo::Array<T, 3>({0.0, 0.0, 0.0});
      allocated_for_output.push_back(sparticle.force_visc);
      sparticle.force_inner_link = new hemo::Array<T, 3>({0.0, 0.0, 0.0});
      allocated_for_output.push_back(sparticle.force_inner_link);
    }
  }

  void HemoCellParticleField::updateResidenceTime(unsigned int rtime)
  {
    for (HemoCellParticle &sparticle : particles)
    {
      sparticle.sv.restime += rtime;
    }
  }

  void HemoCellParticleField::unifyForceVectors()
  {
    for (const hemo::Array<T, 3> *mem : allocated_for_output)
    {
      delete mem;
    }
    allocated_for_output.clear();
    for (HemoCellParticle &sparticle : particles)
    {
      sparticle.repoint_force_vectors();
    }
  }

  void HemoCellParticleField::applyConstitutiveModel(bool forced)
  {
    map<int, vector<HemoCellParticle *>> *ppc_new = new map<int, vector<HemoCellParticle *>>();
    const map<int, vector<int>> &particles_per_cell = get_particles_per_cell();
    map<int, bool> lpc;
    // Fill it here, probably needs optimization, ah well ...
    for (const auto &pair : particles_per_cell)
    {
      const int &cid = pair.first;
      const vector<int> &cell = pair.second;
      (*ppc_new)[cid].resize(cell.size());
      for (unsigned int i = 0; i < cell.size(); i++)
      {
        if (cell[i] == -1)
        {
          (*ppc_new).erase(cid); // not complete, remove entry
          goto no_add_lpc;
        }
        else
        {
          (*ppc_new)[cid][i] = &particles[cell[i]];
        }
      }
      lpc[cid] = true;
    no_add_lpc:;
    }

    for (pluint ctype = 0; ctype < (*cellFields).size(); ctype++)
    {
      if ((*cellFields).hemocell.iter % (*cellFields)[ctype]->timescale == 0 || forced)
      {
        vector<HemoCellParticle *> found;
        findParticles(getBoundingBox(), found, ctype);
        if (found.size() > 0)
        {
          // only reset forces when the forces actually point at it.
          if (found[0]->force_area == &found[0]->sv.force)
          {
            for (HemoCellParticle *particle : found)
            {
              particle->sv.force = {0., 0., 0.};
#ifdef INTERIOR_VISCOSITY
              particle->normalDirection = {0., 0., 0.};
#endif
            }
          }
        }
        (*cellFields)[ctype]->mechanics->ParticleMechanics(*ppc_new, lpc, ctype);
      }
    }

    delete ppc_new;
  }

#define inner_loop                                                                             \
  const int &l_index = grid_index(x, y, z);                                                    \
  const int &n_index = grid_index(xx, yy, zz);                                                 \
  for (unsigned int i = 0; i < particle_grid_size[l_index]; i++)                               \
  {                                                                                            \
    for (unsigned int j = 0; j < particle_grid_size[n_index]; j++)                             \
    {                                                                                          \
      HemoCellParticle &lParticle = particles[particle_grid[l_index][i]];                      \
      HemoCellParticle &nParticle = particles[particle_grid[n_index][j]];                      \
      if (&nParticle == &lParticle)                                                            \
      {                                                                                        \
        continue;                                                                              \
      }                                                                                        \
      if (lParticle.sv.cellId == nParticle.sv.cellId)                                          \
      {                                                                                        \
        continue;                                                                              \
      }                                                                                        \
      const hemo::Array<T, 3> dv = lParticle.sv.position - nParticle.sv.position;              \
      const T distance = sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);                  \
      if (distance < r_cutoff)                                                                 \
      {                                                                                        \
        const hemo::Array<T, 3> rfm = r_const * (1 / (distance / r_cutoff)) * (dv / distance); \
        lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + rfm;                     \
        nParticle.sv.force_repulsion = nParticle.sv.force_repulsion - rfm;                     \
      }                                                                                        \
    }                                                                                          \
  }

#define inner_loop2                                                                                                                                                       \
  bool stochastic_model = 1;                                                                                                                                              \
  const int &l_index = grid_index(x, y, z);                                                                                                                               \
  const int &n_index = grid_index(xx, yy, zz);                                                                                                                            \
  for (unsigned int i = 0; i < particle_grid_size[l_index]; i++)                                                                                                          \
  {                                                                                                                                                                       \
    for (unsigned int j = 0; j < particle_grid_size[n_index]; j++)                                                                                                        \
    {                                                                                                                                                                     \
      HemoCellParticle &lParticle = particles[particle_grid[l_index][i]];                                                                                                 \
      HemoCellParticle &nParticle = particles[particle_grid[n_index][j]];                                                                                                 \
      \  
      if (&nParticle == &lParticle)                                                                                                                                       \
      {                                                                                                                                                                   \
        continue;                                                                                                                                                         \
      }                                                                                                                                                                   \
      if (lParticle.sv.cellId == nParticle.sv.cellId)                                                                                                                     \
      {                                                                                                                                                                   \
        continue;                                                                                                                                                         \
      }                                                                                                                                                                   \
      const hemo::Array<T, 3> dv = lParticle.sv.position - nParticle.sv.position;                                                                                         \
      const T distance = sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);                                                                                             \
      float dist2 = distance * 5 * pow(10, -7);                                                                                                                           \
      float df2 = 6.40625 * pow(10, -9);                                                                                                                                  \
      pluint cells1 = lParticle.sv.celltype;                                                                                                                              \
      pluint cells2 = nParticle.sv.celltype;                                                                                                                              \
      if (distance < 2)                                                                                                                                                   \
      {                                                                                                                                                                   \
        if (stochastic_model)                                                                                                                                             \
        {                                                                                                                                                                 \
          float k_on0 = 1000;                                                                                                                                             \
          float k_off0 = 1;                                                                                                                                               \
          float sigma_on = 0.0000000005;                                                                                                                                  \
          float sigma_off = 0.00000000005;                                                                                                                                \
          float kbT = 4.28 * (pow(10, -21));                                                                                                                              \
          float deltaT = 0.0000001;                                                                                                                                       \
          float P_on;                                                                                                                                                     \
          float P_off;                                                                                                                                                    \
          float k_on = k_on0 * exp((-1 * sigma_on * dist2 * dist2) / (2 * kbT));                                                                                          \
          float k_off = k_off0 * exp(((sigma_off - sigma_on) * dist2 * dist2) / (2 * kbT));                                                                               \
          P_on = (1 - exp(-1 * k_on * deltaT));                                                                                                                           \
          P_off = (1 - exp(-1 * k_off * deltaT));                                                                                                                         \
          float alpha = 1.5;                                                                                                                                              \
          float de = 15 * kbT; /*original value was 15*kbT*/                                                                                                              \
          const hemo::Array<T, 3> rfm = 2 * alpha * de * 1e6 * (exp(2 * alpha * (0.1 - distance * 0.5)) - exp(alpha * (0.1 - distance * 0.5))) * (dv / distance) / df2;   \
          if (distance < 0.4)                                                                                                                                             \
          {                                                                                                                                                               \
            const hemo::Array<T, 3> rfm2 = 4 / (distance * 0.5e-6) * 0.00001 * kbT * (6 * pow(0.4 / distance, 6) - 12 * pow(0.4 / distance, 12)) * (dv / distance) / df2; \
            lParticle.sv.force_repulsion = lParticle.sv.force_repulsion - 1 * rfm2;                                                                                       \
            nParticle.sv.force_repulsion = nParticle.sv.force_repulsion + 1 * (rfm2);                                                                                     \
          }                                                                                                                                                               \
          if ((cells1 == 2 && cells2 == 1))                                                                                                                               \
          {                                                                                                                                                               \
            lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + 1 * rfm / 1;                                                                                    \
            nParticle.sv.force_repulsion = nParticle.sv.force_repulsion - (rfm / 5);                                                                                      \
          }                                                                                                                                                               \
          else if ((cells1 == 1 && cells2 == 2))                                                                                                                          \
          {                                                                                                                                                               \
            lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + (rfm / 5);                                                                                      \
            nParticle.sv.force_repulsion = nParticle.sv.force_repulsion - 1 * rfm / 1;                                                                                    \
          }                                                                                                                                                               \
          else if ((cells1 == 1 && cells2 == 1))                                                                                                                          \
          {                                                                                                                                                               \
            lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + (rfm / 20);                                                                                     \
            nParticle.sv.force_repulsion = nParticle.sv.force_repulsion - (rfm / 20);                                                                                     \
          }                                                                                                                                                               \
          else if ((cells1 == 2 && cells2 == 2))                                                                                                                          \
          {                                                                                                                                                               \
            lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + (rfm);                                                                                          \
            nParticle.sv.force_repulsion = nParticle.sv.force_repulsion - (rfm);                                                                                          \
          }                                                                                                                                                               \
          else if ((cells1 == 2 && cells2 == 3))                                                                                                                          \
          {                                                                                                                                                               \
            lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + (rfm);                                                                                          \
            nParticle.sv.force_repulsion = nParticle.sv.force_repulsion - (rfm);                                                                                          \
          }                                                                                                                                                               \
          else if ((cells1 == 3 && cells2 == 2))                                                                                                                          \
          {                                                                                                                                                               \
            lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + (rfm);                                                                                          \
            nParticle.sv.force_repulsion = nParticle.sv.force_repulsion - (rfm);                                                                                          \
          }                                                                                                                                                               \
        }                                                                                                                                                                 \
      }                                                                                                                                                                   \
    }                                                                                                                                                                     \
  }
  void HemoCellParticleField::applyRepulsionForce(bool forced)
  {
    const T r_const = cellFields->repulsionConstant;
    const T r_cutoff = cellFields->repulsionCutoff;
    plint numberOfParticlesInContact = 0;
    if (!pg_up_to_date)
    {
      update_pg();
    }
    for (HemoCellParticle &particle : particles)
    {
      particle.sv.force_repulsion = {0., 0., 0.};
    }

    for (int x = 0; x < atomicLattice->getNx() - 1; x++)
    {
      for (int y = 0; y < atomicLattice->getNy() - 1; y++)
      {
        for (int z = 0; z < atomicLattice->getNz() - 1; z++)
        {
          // Manual finding, we could make a map, but for now this should be fast enough
          // 0, 0, 0 //0, 1, 0 //0, 0, 1 //0, 1, 1
          int xx = x;
          for (int yy = y; yy <= y + 1; yy++)
          {
            for (int zz = z; zz <= z + 1; zz++)
            {
              inner_loop2
            }
          }
          // 1, 0, 0
          // 1, 1, 0
          // 1, 0, 1
          // 1, 1, 1
          // 1, 0,-1
          // 1,-1, 0
          // 1,-1,-1
          // 1, 1,-1
          // 1,-1, 1
          xx = x + 1;
          for (int yy = y - 1; yy <= y + 1; yy++)
          {
            for (int zz = z - 1; zz <= z + 1; zz++)
            {
              if (yy < 0)
              {
                continue;
              }
              if (zz < 0)
              {
                continue;
              }
              inner_loop2
            }
          }
          xx = x;
          int yy = y + 1, zz = z - 1;
          if (zz < 0)
          {
            continue;
          }
          // 0, 1,-1
          inner_loop2
        }
      }
    }
  }

#ifdef INTERIOR_VISCOSITY
  void HemoCellParticleField::internalGridPointsMembrane(Box3D domain)
  {
    // This could be done less complex I guess?
    for (const HemoCellParticle &particle : particles)
    { // Go over each particle
      if (!(*cellFields)[particle.sv.celltype]->doInteriorViscosity)
      {
        continue;
      }

      for (unsigned int i = 0; i < particle.kernelCoordinates.size(); i++)
      {
        const hemo::Array<T, 3> latPos = particle.kernelCoordinates[i] - (particle.sv.position - atomicLattice->getLocation());
        const hemo::Array<T, 3> &normalP = particle.normalDirection;

        if (computeLength(latPos) > (*cellFields)[particle.sv.celltype]->mechanics->cellConstants.edge_mean_eq)
        {
          continue;
        }

        T dot1 = hemo::dot(latPos, normalP);

        if (dot1 < 0.)
        { // Node is inside
          InteriorViscosityHelper::get(*cellFields).add(*this, {particle.kernelCoordinates[i][0], particle.kernelCoordinates[i][1], particle.kernelCoordinates[i][2]}, (*cellFields)[particle.sv.celltype]->interiorViscosityTau);
          particle.kernelLocations[i]->attributeDynamics((*cellFields)[particle.sv.celltype]->innerViscosityDynamics);
        }
        else
        { // Node is outside
          InteriorViscosityHelper::get(*cellFields).remove(*this, {particle.kernelCoordinates[i][0], particle.kernelCoordinates[i][1], particle.kernelCoordinates[i][2]});
          particle.kernelLocations[i]->attributeDynamics(&atomicLattice->getBackgroundDynamics());
        }
      }
    }
  }

  // For performance reason, this is only executed once every n iterations to make
  // sure that there are no higher viscosity grid points left after substantial movement
  void HemoCellParticleField::findInternalParticleGridPoints(Box3D domain)
  {
    // Reset all the lattice points to the orignal relaxation parameter
    for (const Dot3D &internalPoint : internalPoints)
    {
      atomicLattice->get(internalPoint.x, internalPoint.y, internalPoint.z);
      atomicLattice->get(internalPoint.x, internalPoint.y, internalPoint.z).attributeDynamics(&atomicLattice->getBackgroundDynamics());
    }
    InteriorViscosityHelper::get(*cellFields).empty(*this);

    for (const auto &pair : get_lpc())
    { // Go over each cell?
      const int &cid = pair.first;
      const vector<int> &cell = get_particles_per_cell().at(cid);
      const pluint ctype = particles[cell[0]].sv.celltype;

      // Plt and Wbc now have normal tau internal, so we don't have
      // to raycast these particles
      if (!(*cellFields)[ctype]->doInteriorViscosity)
      {
        continue;
      }

      hemo::OctreeStructCell octCell(3, 1, 30,
                                     (*cellFields)[ctype]->mechanics->cellConstants.triangle_list,
                                     particles, cell);

      std::set<Array<plint, 3>> innerNodes;
      octCell.findInnerNodes(atomicLattice, particles, cell, innerNodes);
      for (const Array<plint, 3> &node : innerNodes)
      {
        InteriorViscosityHelper::get(*cellFields).add(*this, {node[0], node[1], node[2]}, (*cellFields)[ctype]->interiorViscosityTau);
        atomicLattice->get(node[0], node[1], node[2]).attributeDynamics((*cellFields)[ctype]->innerViscosityDynamics);
      }
    }
  }
#else
  void HemoCellParticleField::findInternalParticleGridPoints(Box3D domain)
  {
    pcout << "(HemoCellParticleField) (Error) findInternalParticleGridPoints called, but INTERIOR_VISCOSITY not defined, exiting..." << endl;
    exit(1);
  }
  void HemoCellParticleField::internalGridPointsMembrane(Box3D domain)
  {
    pcout << "(HemoCellParticleField) (Error) internalGridPointsMembrane called, but INTERIOR_VISCOSITY not defined, exiting..." << endl;
    exit(1);
  }
#endif

  void HemoCellParticleField::interpolateFluidVelocity(Box3D domain)
  {
    // Preallocating
    hemo::Array<T, 3> velocity;
    plb::Array<T, 3> velocity_comp;

    for (HemoCellParticle &particle : particles)
    {

      // Trick to allow for different kernels for different particle types.
      // (*cellFields)[particle.sv.celltype]->kernelMethod(*atomicLattice,particle);

      // We have the kernels, now calculate the velocity of the particles.
      velocity = {0.0, 0.0, 0.0};
      for (pluint j = 0; j < particle.kernelLocations.size(); j++)
      {
        // Direct access
        particle.kernelLocations[j]->computeVelocity(velocity_comp);
        velocity += (velocity_comp * particle.kernelWeights[j]);
      }

      if (particle.sv.celltype >= 3)
      {
        continue;
      }
      particle.sv.v = velocity;
    }
  }

  void HemoCellParticleField::determineApoptosisFromConcentration()
  {

    // map<int,CellInformation> info_per_cell;
    // CellInformationFunctionals::calculateCellInformation(&hemocell,info_per_cell);
    // HemoCellGatheringFunctional<CellInformation>::gather(info_per_cell);
    plb::Box3D const box = sourceLattice->getBoundingBox();
    plb::Dot3D const &location = sourceLattice->getLocation();
    plint n = 0;
    plint total = 0;
    for (HemoCellParticle &particle : particles)
    {
      plint x = (particle.sv.position[0] - location.x) + 0.5;
      plint y = (particle.sv.position[1] - location.y) + 0.5;
      plint z = (particle.sv.position[2] - location.z) + 0.5;

      particle.sv.nearbyConcentration = sourceLattice->get(x, y, z).computeDensity();
      total += particle.sv.nearbyConcentration;
      n++;

      T threshold = 1;
      if (particle.sv.nearbyConcentration < threshold)
      {
        particle.sv.cellState = APOPTOSIS;
        particle.tag = 1;
      }
      else
      {
        particle.sv.cellState = MOVEMENT;
      }
    }

    removeParticles(1);
    lpc_up_to_date = false;
    pg_up_to_date = false;
  }

  void HemoCellParticleField::determineImmuneResponseToCTC(HemoCell *hemocell)
  {
    map<int, CellInformation> info_per_cell;
    CellInformationFunctionals::calculateCellInformation(hemocell, info_per_cell);
    HemoCellGatheringFunctional<CellInformation>::gather(info_per_cell);
    unsigned char typeNKCnumber = (*cellFields)["NKC"]->ctype;
    unsigned char typeCTLnumber = (*cellFields)["CTL"]->ctype;
    unsigned char typeCTCnumber = (*cellFields)["CTC"]->ctype;

    vector<CellInformation> NKCList;
    vector<CellInformation> CTCList;
    vector<CellInformation> nearestCTC;

    bool randomWalk = true; //**a debug variable, set to false to have NKCs move directly towards CTCs
    
    //Loops around all the cell information and adds CTCs NKCsto vectors
    for (auto cell = info_per_cell.begin(); cell != info_per_cell.end(); cell++)
    {
      if (cell->second.cellType == typeCTCnumber)
      {
        CTCList.push_back(cell->second);
      }
      else if (cell->second.cellType == typeNKCnumber)
      {
        NKCList.push_back(cell->second);
      }
    }
  
    if (CTCList.empty())
    {
      return;
    }

    //Checking for nearest CTC to each natural killer cell
    for (int i = 0; i < NKCList.size(); i++)
    {
      hemo::Array<T, 3> NKCPos = NKCList[i].position;
      double minDistance = INT_MAX;
      int cellId = 0;
      for (int i = 0; i < CTCList.size(); i++)
      {
        hemo::Array<T, 3> CTCPos = CTCList[i].position;
        T x = CTCPos[0] - NKCPos[0];
        T y = CTCPos[1] - NKCPos[1];
        T z = CTCPos[2] - NKCPos[2];
        T distance = std::sqrt(x * x + y * y + z * z);

        if (distance < minDistance)
        {
          minDistance = distance;
          cellId = CTCList[i].base_cell_id;
        }
      }
      nearestCTC.push_back(info_per_cell[cellId]);
    }

    /*
      Loops through all the CTCs and NKCs to:
      1. Find the total adhesion contact area
      2. Release cytokines at the position of the NKC if they are touching
      3. Kill the CTC if the contact area is above a percentage
    */
    for (int i = 0; i < CTCList.size(); i++)
    {
      hemo::Array<T, 3> CTCPos = CTCList[i].position;
      plint numberOfNKCContact = 0;
      vector<int> CTCparticles = get_particles_per_cell().at(CTCList[i].base_cell_id);

      for (int j = 0; j < NKCList.size(); j++)
      {
        T totalContactArea = 0.0;
        hemo::Array<T, 3> NKCPos = NKCList[j].position;
        T x = NKCPos[0] - CTCPos[0];
        T y = NKCPos[1] - CTCPos[1];
        T z = NKCPos[2] - CTCPos[2];
        T distance = std::sqrt(x * x + y * y + z * z);

        if (distance <= 20) // Stops search if the center of each NKC/CTC is too far, 20 is just a safe estimate
        {
          vector<int> NKCparticles = get_particles_per_cell().at(NKCList[j].base_cell_id);
          for (int CTC = 0; CTC < CTCparticles.size(); CTC++)
          {
            for (int NKC = 0; NKC < NKCparticles.size(); NKC++)
            {
              if (particles[CTCparticles[CTC]].sv.cellId != CTCList[i].base_cell_id)
              {
                std::cout << "INCORRECT PARTICLE FINDING" << std::endl;
              }
              T particleDistance = std::sqrt(
                  (particles[CTCparticles[CTC]].sv.position[0] - particles[NKCparticles[NKC]].sv.position[0]) * (particles[CTCparticles[CTC]].sv.position[0] - particles[NKCparticles[NKC]].sv.position[0]) +
                  (particles[CTCparticles[CTC]].sv.position[1] - particles[NKCparticles[NKC]].sv.position[1]) * (particles[CTCparticles[CTC]].sv.position[1] - particles[NKCparticles[NKC]].sv.position[1]) +
                  (particles[CTCparticles[CTC]].sv.position[2] - particles[NKCparticles[NKC]].sv.position[2]) * (particles[CTCparticles[CTC]].sv.position[2] - particles[NKCparticles[NKC]].sv.position[2]));
              if (particleDistance < 0.2) //<0.2 is set here because the distance never gets closer than 0.1
              {

                const hemo::Array<T, 3> &v0 = particles[NKCparticles[NKC] - (NKCparticles[NKC]) % 3].sv.position;
                const hemo::Array<T, 3> &v1 = particles[NKCparticles[NKC] + 1 - (NKCparticles[NKC]) % 3].sv.position;
                const hemo::Array<T, 3> &v2 = particles[NKCparticles[NKC] + 2 - (NKCparticles[NKC]) % 3].sv.position;
                totalContactArea += 6 * ((T)1 / (T)3) * computeTriangleArea(v0, v1, v2);
               
              }
            }
          }
          //Release cytokines if they are touching
          if (totalContactArea > 0){
            if(hemocell->iter % 100 == 0)
              std::cout << " TOTAL CONTACT AREA: " << totalContactArea << " TOTAL SURFACE AREA: " << NKCList[j].area << " PERCENTAGE: " << totalContactArea / NKCList[j].area << std::endl;
            T concentration = 1000.0;
            Box3D NKCLocation(NKCPos[0], NKCPos[0], NKCPos[1], NKCPos[1], NKCPos[2], NKCPos[2]);
            hemocell->setConcentration(NKCLocation, concentration, plb::Array<T, 3>((T)0., (T)0., (T)0.));
          }
          //Kill CTC if contact area above a certain percentage
          if (totalContactArea/NKCList[j].area > 0.04)
          {
            std::cout << "KILLING CTC, " << " TOTAL CONTACT AREA: " << totalContactArea << " TOTAL SURFACE AREA: " << NKCList[j].area << " PERCENTAGE: " << totalContactArea / NKCList[j].area << std::endl;
            for (HemoCellParticle &particle : particles)
            {
              if (particle.sv.cellId == CTCList[i].base_cell_id)
              {
                particle.tag = 1;
              }
            }
            removeParticles(1);
            lpc_up_to_date = false;
            pg_up_to_date = false;
          }        
        }

        
      }
    }

    //Moves all NKCs, according to a random walk or gradient NOTE: interpolateFluidVelocity is changed to skip NKCs because their velocity will be computed here
    for (int i = 0; i < NKCList.size(); i++)
    {
      //Change the random vector every 1000 timesteps, random_vector is contained as a global var
      if (hemocell->iter % 1000 == 0)
      {
        T randx = (T)rand() / (T)RAND_MAX;
        T randy = (T)rand() / (T)RAND_MAX;
        T randz = (T)rand() / (T)RAND_MAX;
        random_vector = {randx, randy, randz};
      }
      hemo::Array<T, 3> NKCPos = NKCList[i].position;
      hemo::Array<T, 3> CTCPos = nearestCTC[i].position;
      hemo::Array<T, 3> CTCDirection = {CTCPos[0] - NKCPos[0], CTCPos[1] - NKCPos[1], CTCPos[2] - NKCPos[2]};
      T CTCDistance = std::sqrt(CTCDirection[0] * CTCDirection[0] + CTCDirection[1] * CTCDirection[1] + CTCDirection[2] * CTCDirection[2]);
      plint iX, iY, iZ = 0;
      computeGridPosition(NKCPos, &iX, &iY, &iZ);
      T densityMax = 0;
      bool isGradient = false;
      hemo::Array<T, 3> directionOfGradient = {0, 0, 0};
      //The three loops here find the approximate location of highest concentration
      for (plint x = iX - 2; x < iX + 2; x++)
      {
        if (x >= atomicLattice->getNx() || x < 0)
          continue;

        for (plint y = iY - 2; y < iY + 2; y++)
        {
          if (y >= atomicLattice->getNx() || x < 0)
            continue;

          for (plint z = iZ - 2; z < iZ + 2; z++)
          {
            if (z >= atomicLattice->getNx() || x < 0)
              continue;

            T density = sourceLattice->get(x, y, z).computeDensity();
            if (density > densityMax)
            {
              directionOfGradient = {x - NKCPos[0], y - NKCPos[1], z - NKCPos[2]}; // a concentration difference is found
              densityMax = density;
              isGradient = true;
            }
          }
        }
      }

      //Moving the actual particles:
      hemo::Array<T, 3> velocity;
      plb::Array<T, 3> velocity_comp;
      for (HemoCellParticle &particle : particles)
      {

        if (particle.sv.cellId == NKCList[i].base_cell_id && CTCDistance > 15.5)
        {
          //Random walk can be set to false to debug immune functions
          if (randomWalk)
          {

            velocity = {0.0, 0.0, 0.0};
            //Computes the impact of fluid on each particle, and adds either a random direction or the direction of gradient
            for (pluint j = 0; j < particle.kernelLocations.size(); j++)
            {
              // Direct access
              particle.kernelLocations[j]->computeVelocity(velocity_comp);
              //velocity = velocity + (velocity_comp * particle.kernelWeights[j]);
              if (isGradient)
              {
                velocity = velocity + (velocity_comp * particle.kernelWeights[j]) + particle.kernelWeights[j] * 0.00001*directionOfGradient;
              }
              else
              {
                velocity = velocity + (velocity_comp * particle.kernelWeights[j]) + particle.kernelWeights[j] * 0.001*random_vector;
              }
            }

            particle.sv.v = velocity;
          }
          else
          {
            particle.sv.v = 0.005 * CTCDirection; //For debugging purposes, it will have the NKCs move towards CTCs
          }
        }
      }
    }
  }

  void HemoCellParticleField::spreadParticleForce(Box3D domain)
  {
    for (HemoCellParticle &particle : particles)
    {

      // Trick to allow for different kernels for different particle types.
      (*cellFields)[particle.sv.celltype]->kernelMethod(*atomicLattice, particle);

      // Capping force to ensure stability -> NOTE: this can introduce an error if forces are large!
#ifdef FORCE_LIMIT
      const T force_mag = norm(particle.sv.force);
      if (force_mag > param::f_limit)
        particle.sv.force *= param::f_limit / force_mag;
#endif

      // Directly change the force on a node, quick-and-dirty solution.
      for (pluint j = 0; j < particle.kernelLocations.size(); j++)
      {
        // Direct access
        particle.kernelLocations[j]->external.data[0] += ((particle.sv.force_repulsion[0] + particle.sv.force_adhesion[0] + particle.sv.force[0]) * particle.kernelWeights[j]);
        particle.kernelLocations[j]->external.data[1] += ((particle.sv.force_repulsion[1] + particle.sv.force_adhesion[1] + particle.sv.force[1]) * particle.kernelWeights[j]);
        particle.kernelLocations[j]->external.data[2] += ((particle.sv.force_repulsion[2] + particle.sv.force_adhesion[2] + particle.sv.force[2]) * particle.kernelWeights[j]);
      }
    }
  }

  void HemoCellParticleField::populateBoundaryParticles()
  {

    for (int x = 0; x < this->atomicLattice->getNx() - 1; x++)
    {
      for (int y = 0; y < this->atomicLattice->getNy() - 1; y++)
      {
        for (int z = 0; z < this->atomicLattice->getNz() - 1; z++)
        {
          if (this->atomicLattice->get(x, y, z).getDynamics().isBoundary())
          {
            for (int xx = x - 1; xx <= x + 1; xx++)
            {
              if (xx < 0 || xx > this->atomicLattice->getNx() - 1)
              {
                continue;
              }
              for (int yy = y - 1; yy <= y + 1; yy++)
              {
                if (yy < 0 || yy > this->atomicLattice->getNy() - 1)
                {
                  continue;
                }
                for (int zz = z - 1; zz <= z + 1; zz++)
                {
                  if (zz < 0 || zz > this->atomicLattice->getNz() - 1)
                  {
                    continue;
                  }
                  if (!this->atomicLattice->get(xx, yy, zz).getDynamics().isBoundary())
                  {
                    boundaryParticles.push_back({x, y, z});
                    goto end_inner_loop;
                  }
                }
              }
            }
          }
        end_inner_loop:;
        }
      }
    }
  }

  void HemoCellParticleField::applyBoundaryRepulsionForce()
  {
    if (!pg_up_to_date)
    {
      update_pg();
    }
    const T &br_cutoff = cellFields->boundaryRepulsionCutoff;
    const T &br_const = cellFields->boundaryRepulsionConstant;
    for (Dot3D &b_particle : boundaryParticles)
    {
      for (int x = b_particle.x - 1; x <= b_particle.x + 1; x++)
      {
        if (x < 0 || x > this->atomicLattice->getNx() - 1)
        {
          continue;
        }
        for (int y = b_particle.y - 1; y <= b_particle.y + 1; y++)
        {
          if (y < 0 || y > this->atomicLattice->getNy() - 1)
          {
            continue;
          }
          for (int z = b_particle.z - 1; z <= b_particle.z + 1; z++)
          {
            if (z < 0 || z > this->atomicLattice->getNz() - 1)
            {
              continue;
            }
            const int &index = grid_index(x, y, z);
            for (unsigned int i = 0; i < particle_grid_size[index]; i++)
            {
              HemoCellParticle &lParticle = particles[particle_grid[index][i]];
              const hemo::Array<T, 3> dv = lParticle.sv.position - (b_particle + this->atomicLattice->getLocation());
              const T distance = sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
              if (distance < br_cutoff)
              {
                const hemo::Array<T, 3> rfm = br_const * (1 / (distance / br_cutoff)) * (dv / distance);
                lParticle.sv.force_repulsion = lParticle.sv.force_repulsion + rfm;
              }
            }
          }
        }
      }
    }
  }

  void HemoCellParticleField::applyBoundaryAdhesionForce(HemoCell &hemocell)
  {
    if (!pg_up_to_date)
    {
      update_pg();
    }
    for (HemoCellParticle &particle : particles)
    {
      particle.sv.force_adhesion = {0., 0., 0.};
    }
    map<int, CellInformation> info_per_cell;
    CellInformationFunctionals::calculateCellInformation(&hemocell, info_per_cell);
    HemoCellGatheringFunctional<CellInformation>::gather(info_per_cell);
    vector<vector<float>> distancesvectors;
    distancesvectors.clear();
    // vector<vector<float>> bonds;
    //  cout<<"size of bonds" <<bonds.size()<<endl;
    const T &br_cutoff = cellFields->boundaryAdhesionCutoff;
    const T &br_const = cellFields->boundaryAdhesionConstant;
    // Developed adhesion force based on stochastic model by  Nahid Rahmati
    for (auto b_particle = boundaryParticles.begin(); b_particle != boundaryParticles.end(); ++b_particle)
    {
      int __index = b_particle - boundaryParticles.begin();
      int y1 = b_particle->y - 1;
      int y2 = b_particle->y + 1;
      int x0 = b_particle->x;
      int y0 = b_particle->y;
      int z0 = b_particle->z;
      int z1 = b_particle->z - 1;
      int z2 = b_particle->z + 1;

      if (y1 < 0 || y2 > this->atomicLattice->getNy() - 1)
      { // cout<<"boundary:"<<__index<<endl;
        continue;
      }
      if (z1 < 0 || z2 > this->atomicLattice->getNz() - 1)
      { // cout<<"boundary:"<<__index<<endl;
        continue;
      }

      if (this->atomicLattice->get(x0, y1, z0).getDynamics().isBoundary() && atomicLattice->get(x0, y2, z0).getDynamics().isBoundary() && atomicLattice->get(x0, y0, z1).getDynamics().isBoundary() && atomicLattice->get(x0, y0, z2).getDynamics().isBoundary())
      {
        // cout<<"boundary:"<<__index<<";" <<x0 <<";"<<y0<<";"<< z0 <<endl;
        continue;
      }
      for (int x = b_particle->x - 2; x <= b_particle->x + 2; x++)
      {
        if (x < 0 || x > this->atomicLattice->getNx() - 1)
        {
          continue;
        }

        for (int y = b_particle->y - 2; y <= b_particle->y + 2; y++)
        {
          if (y < 0 || y > this->atomicLattice->getNy() - 1)
          {
            continue;
          }
          for (int z = b_particle->z - 2; z <= b_particle->z + 1; z++)
          {
            if (z < 0 || z > this->atomicLattice->getNz() - 2)
            {
              continue;
            }
            if (y > 7 && y < 30 && z > 7 && z < 30)
            {
              continue;
            }
            const int &index = grid_index(x, y, z);

            for (unsigned int i = 0; i < particle_grid_size[index]; i++)
            {
              int in = particle_grid[index][i];
              HemoCellParticle &lParticle = particles[particle_grid[index][i]];
              const hemo::Array<T, 3> dv = lParticle.sv.position - (*b_particle + this->atomicLattice->getLocation());
              int basecellId = info_per_cell[lParticle.sv.cellId].base_cell_id;
              pluint ct = lParticle.sv.celltype;
              if (ct == 0)
              {
                continue;
              }

              const T distance = sqrt(dv[0] * dv[0] + dv[1] * dv[1] + (dv[2]) * (dv[2]));
              const T distance_normal = sqrt((dv[2] - 0.5) * (dv[2] - 0.5));
              if (distance > 1)
              {
                continue;
              }

              float dist2 = distance * 5e-7;
              float dist2_n = distance_normal * 5e-07;
              float df2 = 6.4e-9;
              bool existed_bond = false;
              int a = bonds.size();

              if (bonds.size() != 0)
              {
                for (int j = 0; j < bonds.size(); j++)
                {
                  if ((bonds[j][0] == __index) && (bonds[j][2] == ct) && (bonds[j][3] == lParticle.sv.vertexId) && (bonds[j][4] == basecellId))
                  {
                    bonds[j][5] = distance; // distance
                    bonds[j][6] = distance_normal;
                    existed_bond = true;
                    float P_off;
                    float P_random_break = (float)rand() / (float)(RAND_MAX);

                    if (distance >= 3)
                    {
                      // cout<<"break happen due to oversize distance"<<endl;
                      P_off = 1.0;
                    }
                    else
                    {
                      // calculate P_off for existed bond

                      float k_off0 = 100; // unstressed off rate      10   350                    			   (will be validated in future) xiao 2017

                      float sigma_off = 5e-6;            // off strength       3e-7     0.0000003215625                   			   (will be validated in future) xiao 2017
                      float kbT = 4.14 * (pow(10, -21)); // boltmann constant * temperature                  	(will be validated in future) xiao 2017
                      float deltaT = 1e-7;               // time step of the simulation (time interval)        		   (will be validated in future) xiao 2017
                      float k_off = k_off0 * exp(((sigma_off) * (dist2 - 5e-7) * (dist2 - 5e-7)) / (2 *
                                                                                                    kbT)); // simplification: we used distance instead of (distance - equilibrium spring length )
                      P_off = (1 - exp(-1 * k_off * deltaT));
                      ;
                      // cout<<"K_off: "<<k_off<<"; P_off: "<<P_off<<" dist2: "<< dist2 << endl;
                    }
                    if (P_off > P_random_break)
                    {

                      //   cout<<"break happens"<<"; P_off: "<<P_off<<" dist2: "<< dist2 << " dist2_n: "<< dist2_n<<endl;
                      existed_bond = false;
                      bonds.erase(bonds.begin() + j);
                      j = j - 1;
                    }
                    else
                    {
                      existed_bond = true;
                    }
                  }
                  else if (/*(!(bonds[j][0]==__index) && ((bonds[j][2]==ct) &&(bonds[j][3]==lParticle.sv.vertexId) && (bonds[j][4]==basecellId))) || */ (!(bonds[j][0] == __index) && ((bonds[j][2] == ct) && (bonds[j][3] == lParticle.sv.vertexId) && (bonds[j][4] == basecellId))))
                  {
                    existed_bond = true;
                  }
                }
              } // checking bond breaking finished

              // probabilistic adhesion function (modified by Nahid n3rahmat@uwaterloo.ca)
              if (distance < 3)
              { // br_cutoff

                if (existed_bond == true)
                {
                  continue;
                }
                if (lParticle.sv.celltype != 0)
                {
                  // non-existed bonds that may be formed
                  vector<float> local_a = {__index, in, ct, lParticle.sv.vertexId, basecellId, distance, distance_normal, b_particle->x, b_particle->y, b_particle->z, lParticle.sv.position[0], lParticle.sv.position[1], lParticle.sv.position[2]};
                  distancesvectors.insert(distancesvectors.end(), local_a);
                }
              }
            }
          }
        }
      }
    }

    std::sort(distancesvectors.begin(), distancesvectors.end(), [=](const vector<float> &v1, const vector<float> &v2)
              { return v1[5] < v2[5]; });
    int m = distancesvectors.size();
    int n = distancesvectors[0].size();

    for (int i = 0; i < distancesvectors.size(); i++)
    {
      for (int j = i; j < distancesvectors.size(); j++)
      {
        if (i == j)
        {
          continue;
        }
        if (distancesvectors[i][2] == distancesvectors[j][2] && distancesvectors[i][3] == distancesvectors[j][3] && distancesvectors[i][4] == distancesvectors[j][4])
        {
          distancesvectors.erase(distancesvectors.begin() + j);
          i = 0;
          j = 0;
        }
        // if ((distancesvectors[i][1]==distancesvectors[j][1]) && (distancesvectors[i][2]==distancesvectors[j][2]) && (distancesvectors[i][3]==distancesvectors[j][3])){distancesvectors.erase(distancesvectors.begin()+j); i=0;
        //  j = 0;}
      }
    }
    /* std::sort(distancesvectors.begin(), distancesvectors.end(), [=] (const vector<float>& v1, const vector<float>& v2) {
         return v1[6] < v2[6];
     });*/
    int a = distancesvectors.size();
    // cout<<"distance size"<<a<< endl;
    for (int i = 0; i < distancesvectors.size(); i++)
    {

      float k_on0 = 10000;   // unstressed on rate    12050000     350000                     	       (will be validated in future) xiao 2017
      float sigma_on = 5e-7; // on strength        0.00000065625                         			   (will be validated in future) xiao 2017
      float kbT = 4.14 * (pow(10,
                              -21)); // boltmann constant * temperature                  			   (will be validated in future) xiao 2017
      float deltaT = 1e-7;           // time step of the simulation (time interval)        		   (will be validated in future) xiao 2017
      float P_on;
      float dist2 = distancesvectors[i][5] * 5e-07;
      float dist2_n = distancesvectors[i][6] * 5e-07;
      float k_on = k_on0 * exp((-1 * sigma_on * (dist2 - 5e-7) * (dist2 - 5e-7)) / (2 *
                                                                                    kbT)); // simplification: we used distance instead of (distance - equilibrium spring length )
      P_on = (1 - exp(-1 * k_on * deltaT));
      if (dist2_n < 50e-9)
      { /*P_on =1; */
      }
      // cout<<"K_on: "<<k_on<<"; P_on: "<<P_on<<" dist2: "<< dist2 << endl;
      float P_random_form = (float)rand() / (float)(RAND_MAX);
      if (P_on > P_random_form)
      {
        auto temp = vector<float>{distancesvectors.at(i).begin(), distancesvectors.at(i).begin() + 7};
        bonds.insert(bonds.end(), temp);
      }
    }
    // cout<< "size bond:"<< bonds.size() <<endl;
    for (int i = 0; i < bonds.size(); i++)
    {

      int x1 = boundaryParticles[bonds[i][0]].x;
      int y1 = boundaryParticles[bonds[i][0]].y;
      int z1 = boundaryParticles[bonds[i][0]].z;

      // cout<<"x1:"<< x1 <<"; y1: " << y1 << "; z1:" <<z1 <<endl;
      for (int x = x1 - 2; x <= x1 + 2; x++)
      {
        if (x < 0 || x > this->atomicLattice->getNx() - 1)
        {
          continue;
        }
        for (int y = y1 - 2; y <= y1 + 2; y++)
        {
          if (y < 0 || y > this->atomicLattice->getNy() - 1)
          {
            continue;
          }
          for (int z = z1 - 2; z <= z1 + 2; z++)
          {
            if (z < 0 || z > this->atomicLattice->getNz() - 1)
            {
              continue;
            }
            const int &index = grid_index(x, y, z);
            for (unsigned int j = 0; j < particle_grid_size[index]; j++)
            {
              HemoCellParticle &lParticle = particles[particle_grid[index][j]];
              auto b_particle = boundaryParticles.begin() + bonds[i][0];
              const hemo::Array<T, 3> dv =
                  lParticle.sv.position - (*b_particle + this->atomicLattice->getLocation());
              int basecellId = info_per_cell[lParticle.sv.cellId].base_cell_id;
              pluint ct = lParticle.sv.celltype;
              if (ct == 0)
              {
                continue;
              }
              hemo::Array<T, 3> dv2 = {dv[0], dv[1], dv[2] - 0.5};
              const T distance = sqrt(dv[0] * dv[0] + dv[1] * dv[1] + (dv[2]) * (dv[2]));
              float dist2 = distance * 5e-07;
              float df2 = 6.4e-9;
              // cout<<"x1:"<< dv[0] <<"; y1: " << dv[1] << "; z1:" <<dv[2] <<endl;
              if ((bonds[i][2] == ct) && (bonds[i][3] == lParticle.sv.vertexId) &&
                  (bonds[i][4] == basecellId))
              {
                const hemo::Array<T, 3> rfm =
                    -1e-3 *
                    (dist2 - 5e-7) * (dv2 / distance) / df2;

                if (ct == 1)
                {
                  lParticle.sv.force_adhesion = lParticle.sv.force_adhesion + 1 * (rfm / 1); /// 5
                }
                else
                {
                  lParticle.sv.force_adhesion = lParticle.sv.force_adhesion + 1 * rfm;
                }

                // cout << "rfm: " << rfm[0] << ";" << rfm[1] << ";" << rfm[2] << "; lparticle_vertexid: " << lParticle.sv.vertexId << endl;
              }
            }
          }
        }
      }
    }
    //   cout<<"size of bonds end of iteration: " <<bonds.size()<<endl;
  }

  void HemoCellParticleField::populateBindingSites(plb::Box3D &domain_)
  {

    // Translate domain to make sense in the lattice domain.
    Dot3D shift = this->getLocation() - this->atomicLattice->getLocation();
    Box3D domain = domain_.shift(shift.x, shift.y, shift.z);

    for (int x = domain.x0; x <= domain.x1; x++)
    {
      for (int y = domain.y0; y <= domain.y1; y++)
      {
        for (int z = domain.z0; z <= domain.z1; z++)
        {
          if (this->atomicLattice->get(x, y, z).getDynamics().isBoundary())
          {
            for (int xx = x - 1; xx <= x + 1; xx++)
            {
              if (xx < 0 || xx > this->atomicLattice->getNx() - 1)
              {
                continue;
              }
              for (int yy = y - 1; yy <= y + 1; yy++)
              {
                if (yy < 0 || yy > this->atomicLattice->getNy() - 1)
                {
                  continue;
                }
                for (int zz = z - 1; zz <= z + 1; zz++)
                {
                  if (zz < 0 || zz > this->atomicLattice->getNz() - 1)
                  {
                    continue;
                  }
                  if (!this->atomicLattice->get(xx, yy, zz).getDynamics().isBoundary())
                  {
                    bindingFieldHelper::get(*cellFields).add(*this, {x, y, z});
                    goto end_inner_loop;
                  }
                }
              }
            }
          }
        end_inner_loop:;
        }
      }
    }
  }

  T HemoCellParticleField::eigenValueFromCell(plb::Cell<T, DESCRIPTOR> &cell)
  {
    plb::Array<T, SymmetricTensor<T, DESCRIPTOR>::n> element;
    cell.computePiNeq(element);
    T omega = cell.getDynamics().getOmega();
    T rhoBar = cell.getDynamics().computeRhoBar(cell);
    T prefactor = -omega * DESCRIPTOR<T>::invCs2 *
                  DESCRIPTOR<T>::invRho(rhoBar) / (T)2;
    for (int iTensor = 0; iTensor < SymmetricTensor<T, DESCRIPTOR>::n; ++iTensor)
    {
      element[iTensor] *= prefactor;
    }

    Array<Array<T, 3>, 3> S; // Strain-rate tensor (symmetric).
    S[0][0] = element[0];    // s[xx];
    S[0][1] = element[1];    // s[xy];
    S[0][2] = element[2];    // s[xz];

    S[1][0] = element[1]; // s[yx];
    S[1][1] = element[3]; // s[yy];
    S[1][2] = element[4]; // s[yz];

    S[2][0] = element[2]; // s[zx];
    S[2][1] = element[4]; // s[zy];
    S[2][2] = element[5]; // s[zz];

    Eigen::Matrix<T, 3, 3> A;
    for (plint i = 0; i < 3; i++)
    {
      for (plint j = 0; j < 3; j++)
      {
        A(i, j) = S[i][j];
      }
    }

    bool computeEigenvectors = false;
    Eigen::EigenSolver<Eigen::Matrix<T, 3, 3>> es(A, computeEigenvectors);
    std::vector<T> lambda(3);
    lambda[0] = std::real(es.eigenvalues()[0]);
    lambda[1] = std::real(es.eigenvalues()[1]);
    lambda[2] = std::real(es.eigenvalues()[2]);
    std::sort(lambda.begin(), lambda.end());

    //    Array<Array<T,3>,3> x;  // Eigenvectors of S.
    //    Array<T,3> d;           // Eigenvalues of S.
    //    Eigen::eigenDecomposition(S, x, d);
    //    std::vector<T> lambda(3);
    //    lambda[0] = d[0];
    //    lambda[1] = d[1];
    //    lambda[2] = d[2];
    //   std::sort(lambda.begin(), lambda.end());
    T tresca = (lambda[2] - lambda[0]) / 2;
    return tresca;
  }

  void HemoCellParticleField::prepareSolidification()
  {
#ifdef SOLIDIFY_MECHANICS
    for (HemoCellField *type : cellFields->cellFields)
    {
      ppc_up_to_date = false;
      if (type->doSolidifyMechanics)
      {

        // Invoke cell-type specific solidification mechanics implementation.
        type->mechanics->solidifyMechanics(get_particles_per_cell(), particles, this->atomicLattice, this->CEPAClattice, type->ctype, *this);
      }
    }
#else
    hlog << "(HemoCellParticleField) prepareSolidification called but SOLIDIFY_MECHANICS not enabled" << endl;
    exit(1);
#endif
  }

  void HemoCellParticleField::solidifyCells()
  {
#ifdef SOLIDIFY_MECHANICS

    // Remove any to be removed particles (tagged with `tag == 1`).
    removeParticles(1);
    if (!pg_up_to_date)
    {
      update_pg();
    }

    // Detect particles to be solidified by looping over a block of 3x3 LBM cells
    // around each binding site. When a cell statisfies both:
    // - close enough in space to a binding site,
    // - shows a minimum tresca stress,
    // the particle is labelled to be solified.
    for (const Dot3D &b_particle : bindingSites)
    {
      for (int x = b_particle.x - 1; x <= b_particle.x + 1; x++)
      {
        if (x < 0 || x > this->atomicLattice->getNx() - 1)
        {
          continue;
        }

        for (int y = b_particle.y - 1; y <= b_particle.y + 1; y++)
        {
          if (y < 0 || y > this->atomicLattice->getNy() - 1)
          {
            continue;
          }

          for (int z = b_particle.z - 1; z <= b_particle.z + 1; z++)
          {
            if (z < 0 || z > this->atomicLattice->getNz() - 1)
            {
              continue;
            }

            const int &index = grid_index(x, y, z);

            for (unsigned int i = 0; i < particle_grid_size[index]; i++)
            {
              HemoCellParticle &lParticle = particles[particle_grid[index][i]];
              const hemo::Array<T, 3> dv = lParticle.sv.position - (b_particle + this->atomicLattice->getLocation());
              const T distance = sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
              T tresca = eigenValueFromCell(this->atomicLattice->get(x, y, z));

              // FIXME: both user-defined constants could be extracted outside the loop.
              if ((distance <= (*cellFields)[lParticle.sv.celltype]->mechanics->cfg["MaterialModel"]["distanceThreshold"].read<T>()) && (abs(tresca / 1e-7) > (*cellFields)[lParticle.sv.celltype]->mechanics->cfg["MaterialModel"]["shearThreshold"].read<T>()))
              {
                lParticle.sv.solidify = true;
              }
            }
          }
        }
      }
    }
#else
    hlog << "(HemoCellParticleField) SolidifyCells called but SOLIDIFY_MECHANICS not enabled" << endl;
    exit(1);
#endif
  }

  HemoCellParticleDataTransfer &HemoCellParticleField::getDataTransfer()
  {
    return particleDataTransfer;
  }
  HemoCellParticleDataTransfer const &HemoCellParticleField::getDataTransfer() const
  {
    return particleDataTransfer;
  }

  std::string HemoCellParticleField::getBlockName()
  {
    return std::string("HemoParticleField3D");
  }

  HemoCellFields *HemoCellParticleField::cellFields = 0;
}
