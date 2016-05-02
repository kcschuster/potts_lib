/*
 *  PottsModel.h
 *  POTTS
 *
 *  Created by kelsey
 *
 *
 */

#ifndef PottsModel_H
#define PottsModel_H

#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <map>
#include <list>

namespace potts {
    
    class PottsModel
    {
        
    public:
        PottsModel(int=1000, int=2, double=1.0, int=2, double=1.0);
        virtual ~PottsModel();
        
        void setBeta(double);
        virtual void initLattice(int);
        double calculateTotalEnergy();
        void setSeed(unsigned short*);
        void setJ(double);
        void setMovie(int,std::string);
        double getTime();
        virtual void equilibrate();
        void setEquilibrate(bool);
        virtual void run(double);
        double getFractCrystal(int);
        int getPottsState(int);
        void printSpinCoords(std::string);
        
        
        struct Site
        {
            Site();
            ~Site();
            Site(int);
            Site(const Site&);
			const Site& operator=(const Site&);
			std::vector<int> x;     //vector of x,y,(z) positions
			int s;                  //potts state with values (0, nstate-1)
            std::vector<int> nbr;   //list of nearest-neighbor sites
        };
        
        
        
    protected:
        
        bool _equilibrate;          //true if equilibrating system
        double _time;               //time elapsed
        double _beta;               //inverse temperature
        double _J;                  //energetic interaction between sites of same state
        double _energy;             //total energy of system
        int _z;                     //number of neighbors for each site (=2*d)
        int _nstate;                //total number potts states
        unsigned int _d;            //dimension
        unsigned int _n;            //total number sites
        unsigned int _nsite;        //edge length of simulation box (_n = _nsite^_d)
        unsigned int _makeMovie;    //flag for printing out frames after time step
        double _invN;               //to avoid division, _invN = 1/_n
        
        std::string _spinOut;       //file name for spin output
        
        std::vector<Site> _lattice;
		std::vector<Site> _initial_lattice;
        
        unsigned short _seed[3];    //seed for generating distinct realizations

        
        int deltaFxn(int,int);
        void setLattice(const std::vector<PottsModel::Site>&);
        virtual int attemptMove(int);
        double calcSiteEnergy(int,int);
        
    };
}

#endif