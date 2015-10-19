/*
 *  PottsModel.cpp
 *  POTTS
 *
 *  Implements basic Potts model.
 *
 *  Created by Kelsey Schuster
 *
 */

#include "PottsModel.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <cassert>
#include <string>
#include <map>
#include <list>


namespace potts
{
    PottsModel::Site::Site()
    {
    }
    
    PottsModel::Site::~Site()
    {
    }
    
    PottsModel::Site::Site(int d)
    {
        s = 0;
        x = std::vector<int>(d,0);
    }
    
    PottsModel::Site::Site(const Site& c)
    {
        s = c.s;
        x = c.x;
        nbr = c.nbr;
    }
    
    const PottsModel::Site& PottsModel::Site::operator=(const Site& c)
    {
        s = c.s;
        x = c.x;
        nbr = c.nbr;
        return *this;
    }
    
    //initialize model
    PottsModel::PottsModel(int n, int d, double beta, int nstate, double J)
    : _d(d), _beta(beta), _nstate(nstate), _J(J)
    {
        _nsite = int(pow(n, 1.0/d));
		_n = pow(_nsite, d);
		_lattice.resize(n);
        
        _energy = 0.0;
        
        _seed[0] = 1; _seed[1] = 2; _seed[2] = 3;
        
        //initialize lattice
        initLattice(0);
        
        //when dividing by system size, avoid division later on
        _invN = 1.0/((double)_n);
        
        _time = 0;
        _equilibrate = false;
    }
    
    PottsModel::~PottsModel()
    {
        std::cerr << "here" << std::endl;
    }
    
    //initializes lattice in desired state/states
    void PottsModel::initLattice(int state)
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
        
        //check to make sure valid orientation is chosen for initialization
        if (state > (_nstate-1)) {
            std::cout << "error: invalid state for initialization" << std::endl;
        }
        
        //populate lattice with single desired state
        if (state >= 0) {
            for (unsigned int i=0; i<_n; i++) {
                _lattice[i].s = state;
            }
        } else {
            for (unsigned int i=0; i<_n; i++) {
                _lattice[i].s = 0;
            }
        }
        
        //calculate energy of entire lattice
        _energy = calculateTotalEnergy();
        
        _initial_lattice = _lattice;
        _time = 0.0;
        _equilibrate = false;
    }
    
    //calculate energy of lattice according to effective Hamiltonian
    double PottsModel::calculateTotalEnergy()
    {
        double value = 0.0;
        
        //loop through lattice sites
        for (unsigned int i=0; i<_n; i++) {
            
            int site = _lattice[i].s;
            
            //loop through neighbors
            for (unsigned int j=0; j<_z; j++) {
                
                //add half the value so we don't double-count
                value += 0.5*_J*deltaFxn(site, _lattice[_lattice[i].nbr[j]].s);
            }
        }
        return value;
    }
    
    //evaluates dirac delta functions - returns 1 if same potts state
    int PottsModel::deltaFxn(int i,int j)
    {
        if (i == j) {
            return 1;
        } else {
            return 0;
        }
    }
    
    //sets seed values
    void PottsModel::setSeed(unsigned short* seed)
    {
        _seed[0] = seed[0];
        _seed[1] = seed[1];
        _seed[2] = seed[2];
    }
    
    //sets lattice and saves initial lattice
    void PottsModel::setLattice(const std::vector<PottsModel::Site>& lattice)
    {
        _lattice = lattice;
        _initial_lattice = _lattice;
    }

    //set value of beta
    void PottsModel::setBeta(double beta)
    {
        _beta = beta;
        _energy = calculateTotalEnergy();
    }
    
    //set energetic interaction between two sites of the same state
    void PottsModel::setJ(double J)
    {
        _J = J;
        _energy = calculateTotalEnergy();
    }

    //set make movie parameter
    void PottsModel::setMovie(int m, std::string sos)
    {
        _makeMovie = m;
        _spinOut = sos;
    }
    
    //returns official time
    double PottsModel::getTime()
    {
        return _time;
    }
    
    //equilibrate system
    void PottsModel::equilibrate()
    {
        bool old_value = _equilibrate;
        setEquilibrate(false);
        
        run(10000);
        setEquilibrate(old_value);
        
        _time = 0.0;
    }

    //set equilibrate variable
    void PottsModel::setEquilibrate(bool b)
    {
        _equilibrate = b;
    }
    
    //run simulation
    void PottsModel::run(double nstep)
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
    
    //attempts move
    int PottsModel::attemptMove(int i)
    {
        int currS,newS;
        double currE,finalE,deltaE;
        
        currS = _lattice[i].s;
        
        //randomly choose another orientation to attempt move to
        newS = currS;
        while (currS == newS) {
            newS = nrand48(_seed)%(_nstate);
        }
        
        //calculate change in energy if move accepted
        currE = calcSiteEnergy(currS,i);
        finalE = calcSiteEnergy(newS,i);
        deltaE = finalE - currE;
        
        //accept if satisfied
        if ((deltaE <= 0.0) || (erand48(_seed) < exp(-_beta*deltaE))) {
            _lattice[i].s = newS;
            _energy += deltaE;
        }

        //std::cout << "energy: " << _energy << "\t" << calculateTotalEnergy() << std::endl;
        
        return 1;
    }

    //calculates energy associated with single lattice site (eff. hamiltonian)
    double PottsModel::calcSiteEnergy(int s,int i)
    {
        double value = 0.0;
            
        int site = s;
        
        for (unsigned int j=0; j<_z; j++) {
            value += _J*deltaFxn(site, _lattice[_lattice[i].nbr[j]].s);
        }
        return value;
    }
    
    //print spin coordinates
    void PottsModel::printSpinCoords(std::string s)
    {
        std::ostringstream sfn;
        sfn << s;
        std::ofstream sfile(sfn.str().c_str());
        
        for (unsigned int i=0; i<_n; i++) {
            for (unsigned int j=0; j<_d; j++) {
                sfile << _lattice[i].x[j] << "\t";
            }
            sfile << _lattice[i].s << "\n";
        }
    }
    
    //returns fraction of sites that are specific phase of crystal
    double PottsModel::getFractCrystal(int state)
    {
        double s = 0.0;
        for (unsigned int i=0; i<_n; i++) {
            if (_lattice[i].s == state) {
                s += 1.0;
            }
        }
        return s*_invN;
    }
    
    //returns potts state of site i
    int PottsModel::getPottsState(int i)
    {
        return _lattice[i].s;
    }
}

