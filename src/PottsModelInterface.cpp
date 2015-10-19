/*
 *  PottsModelInterface.cpp
 *  POTTS
 *
 *  Implements Potts model with two states (forming an interface) and a
 *  harmonic constraint on the density of one of the states.
 *
 *  Created by Kelsey Schuster
 *
 */

#include "PottsModelInterface.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <cassert>
#include <string>
#include <map>
#include <list>

using namespace std;

namespace potts
{
    PottsModelInterface::PottsModelInterface(int n, int d, double beta, int nstate, double J, int s1, int s2, double kappa)
    : PottsModel(n, d, beta, nstate, J)
    {
        //nstate must equal 2
        if (nstate != 2) {
            std::cerr << "Error: must have nstate=2 for interface initialization" << std::endl;
            exit(1);
        }
        
        //set states for desired interface - state1 has umbrella potential
        _state1 = s1;
        _state2 = s2;
        _kappa = kappa;
        
        //initialize lattice
        initLattice();
        
        _time = 0;
        _equilibrate = false;
    }
    
    
    //initializes lattice with interface
    void PottsModelInterface::initLattice()
    {
        _lattice.clear();
        for (unsigned int i=0; i<_n; i++) {
            Site c(_d);
            _lattice.push_back(c);
        }
        
        //d=2 square lattice
        if (_d == 2) {
            for (int i=0; i<_nsite; i++) {
                for (int j=0; j<_nsite; j++) {
                    
                    //get site index
                    int c = j + i*_nsite;
                    _lattice[c].nbr.clear();
                    
                    //set x, y position
                    _lattice[c].x[0] = i;
                    _lattice[c].x[1] = j;
                    
                    //build neighbor list
                    for (int x=-1; x<=1; x++) {
                        for (int y=-1; y<= 1; y++) {
                            
                            //only one dimension can change at a time
                            if (abs(x+y) != 1) {
                                continue;
                            }
                            int neighbor =  (i+x+_nsite)%_nsite*_nsite +
                            (j+y+_nsite)%_nsite;
                            if (neighbor != c) {
                                _lattice[c].nbr.push_back(neighbor);
                            }
                        }
                    }
                }
            }
        }
        
        //d=3 square lattice
        if (_d == 3) {
            for (int i=0; i<_nsite; i++) {
                for (int j=0; j<_nsite; j++) {
                    for (int k=0; k<_nsite; k++) {
                        
                        //get site index
                        int c = k + j*_nsite + i*_nsite*_nsite;
                        _lattice[c].nbr.clear();
                        
                        //set x, y, z position
                        _lattice[c].x[0] = i;
                        _lattice[c].x[1] = j;
                        _lattice[c].x[2] = k;
                        
                        //build neighbor list
                        for (int x=-1; x<=1; x++) {
                            for (int y=-1; y<=1; y++) {
                                for (int z=-1; z<=1; z++) {
                                    
                                    //only one dim can change at a time
                                    if (!(x==0 && z==0) && !(x==0 && y==0) && !(y==0 && z==0)) {
                                        continue;
                                    }
                                    int neighbor = (i+x+_nsite)%_nsite*_nsite*_nsite
                                    + (j+y+_nsite)%_nsite*_nsite + (k+z+_nsite)%_nsite;
                                    if (neighbor != c) {
                                        _lattice[c].nbr.push_back(neighbor);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        //set number of neighbors as variable
        _z = _lattice[0].nbr.size();
        
        //check to make sure we have correct number of neighbors
        if (_z != 2*_d) {
            std::cout << "error: num neighbors not equal to 2*d" << std::endl;
        }
        
        //add Potts states to make interface 50/50
        for (unsigned int i=0; i<_n; i++) {
            if (_lattice[i].x[0] < _nsite*0.25 || _lattice[i].x[0] > _nsite*0.75) {
                _lattice[i].s = _state1;
            } else {
                _lattice[i].s = _state2;
            }
        }
        
        //print system composition
        std::cout << "fraction state = " << _state1 << ": " << getFractCrystal(_state1) << "\nfraction state = " << _state2 << ": " << getFractCrystal(_state2) << std::endl;
        
        //calculate energy of entire lattice
        _energy = calculateTotalEnergy();
        
        //get fraction/density of each state
        _fractionS1 = getFractCrystal(_state1);
        _fractionS2 = getFractCrystal(_state2);
        
        _initial_lattice = _lattice;
        _time = 0.0;
        _equilibrate = false;
    }
    
    //runs model with interface
    void PottsModelInterface::run(double nstep)
    {
        for (int t=0; t<nstep; t++) {
            for (int ii=0; ii<_n; ii++) {
                
                //choose site and attempt move
                int i = nrand48(_seed) % _n;
                attemptMove(i);
                _time++;
            }
            
            //print data for movie
            if (_makeMovie) {
                printSpinCoords(_spinOut);
            }
            
            //std::cout << "energy: " << _energy << "\t" << calculateTotalEnergy() << std::endl;
        }
    }
    
    //attempts move and uses umbrella potential on density to stabilize interface
    int PottsModelInterface::attemptMove(int i)
    {
        int currS,newS;
        double currE,finalE,deltaE,dev;
        int s1Flag = 0;
        int s2Flag = 0;
        
        //selected site for which to attempt move
        currS = _lattice[i].s;
        
        //randomly choose another orientation to attempt move to
        newS = currS;
        while (newS == currS) {
            newS = nrand48(_seed)%(_nstate);
        }
            
        //calculate change in energy resulting from move
        currE = calcSiteEnergy(currS,i);
        finalE = calcSiteEnergy(newS,i);
        deltaE = finalE - currE;
            
        //calculate change in chosen state fraction
        if (newS == _state1) {
            s1Flag = 1;
            _fractionS1 += 1.0*_invN;
        }
        else if (currS == _state1) {
            s1Flag = -1;
            _fractionS1 -= 1.0*_invN;
        }
        //deviation from 50% (where harmonic potential is centered)
        dev = (_fractionS1 - 0.5);
            
        // update fract of other state in system
        if (newS == _state2) {
            s2Flag = 1;
            _fractionS2 += 1.0*_invN;
        }
        else if (currS == 1) {
            s2Flag = -1;
            _fractionS2 -= 1.0*_invN;
        }
            
        //accept or reject move by metropolis criterion with harmonic constraint
        if ((deltaE <= 0.0) || (erand48(_seed) < exp(-_beta*(deltaE + _kappa*dev*dev)))) {
                
            //make move and update energy
            _lattice[i].s = newS;
            _energy += deltaE;
                
        } else {
                
            //reset fractions to original values if move rejected
            _fractionS1 -= ((double)s1Flag)*_invN;
            _fractionS2 -= ((double)s2Flag)*_invN;
        }
    
        //std::cout << "energy: " << _energy << "\t" << calculateTotalEnergy() << std::endl;
        
        return 1;
    }
    
    //equilibrate system
    void PottsModelInterface::equilibrate()
    {
        bool old_value = _equilibrate;
        setEquilibrate(false);
        
        run(1000);
        setEquilibrate(old_value);
        
        _time = 0.0;
    }
    
    //return state fraction variables
    double PottsModelInterface::getStateFractions(int s)
    {
        if (s == 1) {
            return _fractionS1;
        }
        else if (s == 2) {
            return _fractionS2;
        }
        return 0;
    }
}


