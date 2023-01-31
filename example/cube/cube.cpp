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
#include <hemocell.h>
#include <helper/voxelizeDomain.h>
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "wbcHighOrderModel.h"
#include "helper/hemocellInit.hh"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include "writeCellInfoCSV.h"
#include <fenv.h>

#include "palabos3D.h"
#include "palabos3D.hh"

typedef double T;

using namespace hemo;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }

  HemoCell hemocell(argv[1], argc, argv);
  Config *cfg = hemocell.cfg;

  // number of cells along each axis
  int nx, ny, nz;
  nx = (*cfg)["domain"]["nx"].read<int>();
  ny = (*cfg)["domain"]["ny"].read<int>();
  nz = (*cfg)["domain"]["nz"].read<int>();

  hlog << "(unbounded) (Parameters) calculating flow parameters" << endl;
  param::lbm_shear_parameters((*cfg), nz);
  param::printParameters();

  hlog << "(unbounded) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.initializeLattice(defaultMultiBlockPolicy3D().getMultiBlockManagement(nx, ny, nz, (*cfg)["domain"]["fluidEnvelope"].read<int>()));

  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
                = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

/*  hemocell.lattice->toggleInternalStatistics(false);

  // extract sides of the rectangular domain for assignment of boundary
  // conditions along the outer planes of the domain
  Box3D top = Box3D(0, nx-1, 0, ny-1, nz-1, nz-1);
  Box3D bottom = Box3D(0, nx-1, 0, ny-1, 0, 0);
//  Box3D front = Box3D(0, nx-1, 0, 0, 0, nz-1);
//  Box3D back  = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);

  // all directions have periodicity
  lattice.periodicity().toggleAll(true);

  // bounce back conditions along the front and back of the domain
  defineDynamics(*hemocell.lattice, front, new BounceBack<T, DESCRIPTOR> );
  defineDynamics(*hemocell.lattice, back, new BounceBack<T, DESCRIPTOR> );

  // assign velocity at top/bottom
  (*boundaryCondition).setVelocityConditionOnBlockBoundaries ( *hemocell.lattice, top );
  (*boundaryCondition).setVelocityConditionOnBlockBoundaries ( *hemocell.lattice, bottom );

  // define shear velocity along top/bottom planes (z axis)
  // shear velocity given by `height * shear rate / 2`
  T vtop = (nz-1)*param::shearrate_lbm;
  setBoundaryVelocity(*hemocell.lattice, top, plb::Array<T,3>(vtop,0.0,0.0));
  setBoundaryVelocity(*hemocell.lattice, bottom, plb::Array<T,3>(0.0,0.0,0.0));

  hemocell.latticeEquilibrium(1., plb::Array<T, 3>(0.0, 0.0, 0.0));

  // report basic information regarding the current multi-block configuration
  hlog << getMultiBlockInfo(*hemocell.lattice) << endl;

  // assign force vector equal in x, y, z directions
  setExternalVector(*hemocell.lattice, (*hemocell.lattice).getBoundingBox(),
                    DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                    plb::Array<T, DESCRIPTOR<T>::d>(
                        0.0, 0.0, 0.0)); */
  hemocell.lattice->toggleInternalStatistics(false);

  iniLatticeSquareCouette(*hemocell.lattice, nx, ny, nz, *boundaryCondition, param::shearrate_lbm);

  hemocell.lattice->initialize();

  // initialise the lattice
  hemocell.lattice->initialize();

  // initialise the cells
  hemocell.initializeCellfield();

  // the desired RBC type
  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);

  // define update increments
  hemocell.setMaterialTimeScaleSeparation(
      "RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setParticleVelocityUpdateTimeScaleSeparation(
      (*cfg)["ibm"]["stepParticleEvery"].read<int>());

  hemocell.addCellType<WbcHighOrderModel>("WBC", WBC_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("WBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());


    hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

    hemocell.setRepulsion((*cfg)["domain"]["kRep"].read<T>(), (*cfg)["domain"]["RepCutoff"].read<T>());
    hemocell.setRepulsionTimeScaleSeperation((*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  // hemocell output fields
  vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC,OUTPUT_FORCE_REPULSION,OUTPUT_FORCE_ADHESION,OUTPUT_CELL_ID, OUTPUT_VERTEX_ID};

  hemocell.setOutputs("RBC", outputs);
  hemocell.setOutputs("WBC", outputs);
  // LBM fluid output fields
  outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE};
  hemocell.setFluidOutputs(outputs);




  // Turn on periodicity in the `x`-direction: along shear
  hemocell.setSystemPeriodicity(0, true);

  // Enable boundary particles
  hemocell.enableBoundaryParticles((*cfg)["domain"]["BkRep"].read<T>(), (*cfg)["domain"]["BRepCutoff"].read<T>(),(*cfg)["domain"]["BkAdh"].read<T>(), (*cfg)["domain"]["BAdhCutoff"].read<T>(),(*cfg)["ibm"]["stepMaterialEvery"].read<int>());

  // loading the cellfield
  if (not cfg->checkpointed) {
    hemocell.loadParticles();
    hemocell.writeOutput();
  } else {
    hemocell.loadCheckPoint();
  }

  if (hemocell.iter == 0) {
    hlog << "(unbounded) fresh start: warming up cell-free fluid domain for "
         << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..."
         << endl;
    for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>();
         ++itrt) {
      hemocell.lattice->collideAndStream();
    }
  }

  unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
  unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();

  hlog << "(unbounded) Starting simulation..." << endl;
  int ncells =
      CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC");
  hlog << " | RBC Volume ratio [x100%]: "
       << ncells * 77.0 * 100 / (nx * ny * nz) << endl;
  hlog << "(main)   nCells (global) = " << ncells << endl;

  while (hemocell.iter < tmax) {
    hemocell.iterate();

    if (hemocell.iter % tmeas == 0) {
      hlog << "(main) Stats. @ " << hemocell.iter << " ("
           << hemocell.iter * param::dt << " s):" << endl;
      hlog << "\t # of cells: "
           << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
      hlog << " | # of RBC: "
           << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell,
                                                                   "RBC");
      hlog << " | # of WBC: "
          << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell,
                            "WBC");
      FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell);
      double toMpS = param::dx / param::dt;
      hlog << "\t Velocity  -  max.: " << finfo.max * toMpS
           << " m/s, mean: " << finfo.avg * toMpS
           << " m/s, rel. app. viscosity: "
           << (param::u_lbm_max * 0.5) / finfo.avg << endl;

      hemocell.writeOutput();
    }
  }

  hemo::global.statistics.printStatistics();
  hemo::global.statistics.outputStatistics();

  hlog << "(main) Simulation finished :) " << endl;
  return 0;
}
