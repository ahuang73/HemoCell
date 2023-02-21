// Include most of the interface offered by the HemoCell library
#include <hemocell.h>
// This is the mechanical model for the cells that we want to use later on,
// alternatives can be found in the mechanics folder
#include <rbcHighOrderModel.h>
// These are functions found in the helpers folder, they are not in the core of
// HemoCell but can be handy nonetheless
#include <cellInfo.h>
#include <fluidInfo.h>
#include "pltSimpleModel.h"
#include "particleInfo.h"
#include <fenv.h>
#include "palabos3D.h"
#include "palabos3D.hh"


using namespace std;
using namespace hemo;

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
        return -1;
    }
    // The first argument is the config.xml location, the second and third argument
    // are necessary as a passthrough for the Palabos initialization
    HemoCell hemocell(argv[1], argc, argv);
    Config * cfg = hemocell.cfg;

    // Calculate and load in the lattice Boltzmann parameters from the config file
    // that will be used later on. Pretend that we are calculating the parameters
    // for a pipe, to get an acceptable maximum velocity.

    int refDirN = (*cfg)["domain"]["refDirN"].read<int>();
    param::lbm_pipe_parameters((*hemocell.cfg), refDirN);
    // Also print the parameters so we have visual confirmation.
    param::printParameters();

    // Although we are not creating a pipe, we still must define a driving force,
    // We pretend that this is a pipe, therefore the resulting velocity will be higher,
    // but acceptable. It is possible to analytically solve this correctly if you
    // want.
    T poiseuilleForce = 8 * param::nu_lbm * (param::u_lbm_max * 0.5) / param::pipe_radius / param::pipe_radius;
    // First we create a Palabos management object
    // The first three arguments are the number of fluid cells in x,y and z
    // direction, so this is a 50x50x50 block, the fourth argument is the fluid
    // envelope size and must be two

    // Initialize the fluid lattice within hemocell
    hemocell.lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(
        defaultMultiBlockPolicy3D().getMultiBlockManagement(refDirN, refDirN, refDirN, (*cfg)["domain"]["fluidEnvelope"].read<int>()),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
        defaultMultiBlockPolicy3D().getCombinedStatistics(),
        defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
        new GuoExternalForceBGKdynamics<T, DESCRIPTOR>(1.0 / param::tau));

    // Just to be sure disable all periodicity. Afterwards enable it in the
    // x-direction
    hemocell.lattice->periodicity().toggleAll(false);
    hemocell.lattice->periodicity().toggle(0, true);
    // Set up bounceback boundaries in the other directions
    Box3D topChannel(0, refDirN - 1, 0, refDirN - 1, refDirN - 1, refDirN - 1);
    Box3D bottomChannel(0, refDirN - 1, 0, refDirN - 1, 0, 0);
    Box3D backChannel(0, refDirN - 1, refDirN - 1, refDirN - 1, 0, refDirN - 1);
    Box3D frontChannel(0, refDirN - 1, 0, 0, 0, refDirN - 1);
    
  
    
    defineDynamics(*hemocell.lattice, topChannel, new BounceBack<T, DESCRIPTOR>);
    defineDynamics(*hemocell.lattice, bottomChannel, new BounceBack<T, DESCRIPTOR>);
    defineDynamics(*hemocell.lattice, backChannel, new BounceBack<T, DESCRIPTOR>);
    defineDynamics(*hemocell.lattice, frontChannel, new BounceBack<T, DESCRIPTOR>);
    // Disable statistics to run faster
    hemocell.lattice->toggleInternalStatistics(false);
    // Equilibrate everything
    hemocell.latticeEquilibrium(1., plb::Array<double, 3>(0., 0., 0.));
    Box3D source ((*cfg)["domain"]["sourcex1"].read<int>(),(*cfg)["domain"]["sourcex2"].read<int>(),
    (*cfg)["domain"]["sourcey1"].read<int>(),
    (*cfg)["domain"]["sourcey2"].read<int>(),
    (*cfg)["domain"]["sourcez1"].read<int>(),(*cfg)["domain"]["sourcez2"].read<int>());
   
    

    // After we set up the fluid, it is time to set up the particles in the
    // simulation

    hemocell.initializeCellfield();

    // Add a particleType to the simulation, the template argument refers to the
    // corresponding mechanics in the mechanics/ folder
    // The first argument must correspond with the CELL.xml and CELL.pos present in
    // the directory (where CELL is the string input).
    // The second argument defines how a cell is build up. see
    // config/constant_defaults.h for options.
    hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);

    // Only update the forces resulting from the mechanical deformation every X
    // timesteps, recalculating this is the most costly step and since our
    // timestep is so small it can be done intermittently
    hemocell.setMaterialTimeScaleSeparation("RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
    
    //double check what this does
    hemocell.setInitialMinimumDistanceFromSolid("RBC", 0.5); 
    
    // Only update the integrated velocity (from the fluid field to the particles)
    // every X timesteps.
    hemocell.setParticleVelocityUpdateTimeScaleSeparation(5);

    
    hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
    hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());

    vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC};
    hemocell.setOutputs("RBC", outputs);
    hemocell.setOutputs("PLT", outputs);
    hemocell.setFluidOutputs({OUTPUT_VELOCITY, OUTPUT_DENSITY, OUTPUT_FORCE,
                              OUTPUT_SHEAR_RATE, OUTPUT_STRAIN_RATE,
                              OUTPUT_SHEAR_STRESS});

    hemocell.setSourceOutputs({OUTPUT_DENSITY});
    plb::initializeAtEquilibrium(*hemocell.cellfields->sourceLattice, (*hemocell.cellfields->sourceLattice).getBoundingBox(), 0.0, {0.0,0.0,0.0});

    OnLatticeAdvectionDiffusionBoundaryCondition3D<T,CEPAC_DESCRIPTOR>* diffusionBoundary = createLocalAdvectionDiffusionBoundaryCondition3D<T, CEPAC_DESCRIPTOR>();
    diffusionBoundary->addTemperatureBoundary2N(source,*hemocell.cellfields->sourceLattice, boundary::density);
    plb::setBoundaryDensity(*hemocell.cellfields->sourceLattice,source,0.5);
    
    hemocell.lattice->initialize();
    hemocell.cellfields->sourceLattice->initialize();

    hemocell.enableBoundaryParticles((*cfg)["domain"]["kRep"].read<T>(), (*cfg)["domain"]["BRepCutoff"].read<T>(),(*cfg)["ibm"]["stepMaterialEvery"].read<int>());

    // Turn on periodicity in the X direction
    hemocell.setSystemPeriodicity(0, true);

    // Load the particles from all the *.pos files
    hemocell.loadParticles();

    // Load some basic values from the config.xml file that define how long the
    // simulation must run and when we want to save output
    unsigned int tmax = (*hemocell.cfg)["sim"]["tmax"].read<unsigned int>();
    unsigned int tmeas = (*hemocell.cfg)["sim"]["tmeas"].read<unsigned int>();
    // This is the main running loop, run for tmax iterations.
    while (hemocell.iter < tmax)
    {
        // Advance the fluid field and cellfields one tick.
        hemocell.iterate();

        // Set driving force as required after each iteration
        setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                          DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                          plb::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));

        // When we want to save
        if (hemocell.iter % tmeas == 0)
        {
            hemocell.writeOutput();
        }
    }
    return 0;
}