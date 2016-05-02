/*
 *  PottsModelRandom.cpp
 *  POTTS
 *
 *  Implements Potts model with two states randomly initialized and a
 *  harmonic constraint on the density of one of the states.
 *
 *  Created by Kelsey Schuster
 *
 */

#include "PottsModelRandom.h"

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
    PottsModelRandom::PottsModelRandom(int n, int d, double beta, int nstate, double J, int s1, int s2, double kappa)
    : PottsModelInterface(n, d, beta, nstate, J, s1, s2, kappa)
    {
        
        //initialize lattice
        initLattice();
        
        _time = 0;
        _equilibrate = false;
    }
    
    PottsModelRandom::~PottsModelRandom()
    {
        cerr << "~PottsModelRandom" << endl;
    }
    
    
    //initializes lattice with interface
    void PottsModelRandom::initLattice()
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
        
        //add Potts states to make interface
        for (unsigned int i=0; i<_n; i++) {
            if (erand48(_seed) < 0.5) {
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
}
