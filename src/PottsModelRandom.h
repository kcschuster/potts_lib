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

#include "PottsModel.h"
#include "PottsModelInterface.h"

namespace potts {
    
    class PottsModelRandom : public PottsModelInterface
    {
    public:
        PottsModelRandom(int n=1000, int d=2, double beta=1, int nstate=2, double J=1, int s1=0, int s2=1, double kappa=1.0);
        
    protected:
        void initLattice();
    };
}

#endif
