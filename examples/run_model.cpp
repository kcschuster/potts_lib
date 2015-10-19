/*
 *  run_model.cpp
 *
 *  Runs Potts Model simulation and prints configurations
 *  of states.
 *
 *  Created by Kelsey Schuster
 *
 */

#include <iostream>
#include <sstream>
#include <potts/potts.h>
#include <fstream>
#include <stdio.h>
#include <math.h>

using namespace std;
using namespace potts;


// main function
int main()
{
    //==================================================================
    
    //set system size and dimension
    int systemSize = 1000;
    int dimension = 3;
    int nstate = 2;
    
    //set model parameters
    double beta = 2.0;
    double J = 1.0;
    
    //==================================================================
    
    
    //initialize model
    PottsModel A(systemSize, dimension, beta, nstate, J);
    
    //initialize in state 1 (defaults to state 0)
    A.initLattice(1);
    
    //equilibrate
    A.equilibrate();
    
    
    //run model
    for (int i=0; i<10; i++) {
        
        //run for 100 sweeps
        A.run(100);
        
        if (i % 10 == 0) {
            cout << i << endl;
        }
        
        //print current spin coordinates
        ostringstream file;
        file << "spins" << i << ".txt";
        A.printSpinCoords(file.str());
    }    
}


