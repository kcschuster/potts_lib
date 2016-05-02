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
        PottsModelInterface(int=1000, int=2, double=1, int=2, double=1, int=0, int=1, double=1.0);
        virtual ~PottsModelInterface();
        
        void run(double);
        void equilibrate();
        double getStateFractions(int);
        
    protected:
        virtual void initLattice();
        int attemptMove(int);
        
        int _state1;
        int _state2;
        double _kappa;
        
        double _fractionS1;
        double _fractionS2;
    };
}

#endif