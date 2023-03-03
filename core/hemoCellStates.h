#include "hemoCellParticle.h"
#include "hemocell.h"
#include "hemoCellFields.h"
#include "palabos3D.h"
#include "palabos3D.hh"

namespace hemo {
    class HemoCellState{
        public:
        unsigned char state = 0;
        T oxygenConcentration = 0.0;
        HemoCellState();
        void determineApoptosis(Box3D domain);
        void calculateNearbyOxygenConcentration(Box3D domain);


    };
    
}