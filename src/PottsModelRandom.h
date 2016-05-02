/*
 *  PottsModelRandom.h
 *  POTTS
 *
 *  Implements Potts model with two states distributed randomly and a
 *  harmonic constraint on the density of one of the states.
 *
 *  Created by Kelsey Schuster
 *
 */

#ifndef PottsModelRandom_H
#define PottsModelRandom_H

#include "PottsModelInterface.h"

namespace potts {
    
    class PottsModelRandom : public PottsModelInterface
    {
    public:
        PottsModelRandom(int=1000, int=2, double=1, int=2, double=1, int=0, int=1, double=1.0);
        ~PottsModelRandom();
        
    protected:
        void initLattice();
    };
}

#endif
