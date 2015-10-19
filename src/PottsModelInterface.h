/*
 *  PottsModelInterface.h
 *  POTTS
 *
 *  Implements Potts model with two states (forming an interface) and a
 *  harmonic constraint on the density of one of the states.
 *
 *  Created by Kelsey Schuster
 *
 */

#ifndef PottsModelInterface_H
#define PottsModelInterface_H

#include "PottsModel.h"

namespace potts {
    
    class PottsModelInterface : public PottsModel
    {
    public:
        PottsModelInterface(int n=1000, int d=2, double beta=1, int nstate=2, double J=1, int s1=0, int s2=1, double kappa=1.0);
        
        void run(double);
        void equilibrate();
        double getStateFractions(int);
        
    protected:
        void initLattice();
        int attemptMove(int);
        
        int _state1;
        int _state2;
        double _kappa;
        
        double _fractionS1;
        double _fractionS2;
    };
}

#endif