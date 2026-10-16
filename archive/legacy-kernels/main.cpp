#include "ITvariables.hpp"
#include "ITlattice.hpp"
#include "ITupdates.hpp"			
#include "ITobservables.hpp"
#include "ITparalleltemp.hpp"
#include "writeresults.hpp"
#include "mpi.h"
#include <iomanip>
#include <cmath>
#include <math.h>
#include <iostream>
#include <set>

static ran rann;

int main()
{
	ITvariables vars;	
	ITupdates update;
	ITparalleltemp partemp;
	ITobservables obs;
	writeresults write;

	MPI::Status status;	
	MPI::Init();
	myrank = MPI::COMM_WORLD.Get_rank();
	numprocs = MPI::COMM_WORLD.Get_size();	

	vars.declare_variables(myrank);
	vars.set_temperatures(myrank);

  	ITlattice *** lattice;
  	
  	//initialize spins randomly
  	lattice = new ITlattice ** [M];
  	for (int xc = 0; xc < M; ++xc) {
  	  lattice[xc] = new ITlattice * [Ns];
  	  for (int yc = 0; yc < Ns; ++yc) {
  	    lattice[xc][yc] = new ITlattice[Nl];
  	    for (int zc = 0; zc < Nl; ++zc) { 
  	      lattice[xc][yc][zc].set_S((rann() <= 0.5) ? -1 : 1);
  	    }
  	    /*for (int zc = 1; zc < Nl; ++zc) { 
  	      lattice[xc][yc][zc].set_S(lattice[xc][yc][0].S());
  	    }  */	    
  	  }
  	}
  	/*
  	for (int xc = 1; xc < M; ++xc) {
  	  lattice[xc] = new ITlattice * [Ns];
  	  for (int yc = 0; yc < Ns; ++yc) {
  	    lattice[xc][yc] = new ITlattice[Nl];
  	    for (int zc = 0; zc < Nl; ++zc) { 
  	      lattice[xc][yc][zc].set_S( lattice[0][yc][zc].S() );
  	    } 	    
  	  }
  	}
  	*/
  	ITlattice * plaquette;
  	plaquette = new ITlattice[Np];
  	for (int i = 0; i < Np; i++) {plaquette[i].set_S(-1);}

	double t1, t2;
	if (myrank == 0)
	{
		t1 = MPI::Wtime();
	}
	
	for (int i = 0; i < numprocs; ++i){
	partemp.ID[i]=i;
	}
	exc_index=0;	
	if(bool_J3){J3=J3_temp;}
	//cout << "Equilibration done    " << h << "    "  << Kx << endl;	
	for (int i = 0; i < isteps; ++i)
	{
		//Equilibration      
		update.swendsen_flip(lattice, plaquette, T);
		//update.E=update.enrg(lattice, 1, 1);
		//partemp.parallel_tempering(i, update);		
		//partemp.paralleltemp(i, update);			
	}
	

	//cout << "Equilibration done    " << myrank << endl;	
	//obs.pbc();
	obs.Initiate_observables();
			
	
	for (int j = 0; j < nbins; ++j)
	{
		// Measurements
		for (int i=0; i < iter; ++i)
		{
			update.swendsen_flip(lattice, plaquette, T);
			//update.E=update.enrg(lattice, 1, 1);			
			//partemp.parallel_tempering(i, update);			
			obs.observables(lattice, update);
		}
		obs.binning_data(nbins, iter);
	}
	
	//cout << "observables done    " << myrank << "  "<< Kx << endl;	
	if (versus == "TF" or versus == "T") {
	write.output(myrank, T, obs.wdata1[0], (obs.wdata1[1]-pow(obs.wdata1[0],2))/T, obs.wdata1[3],  (obs.wdata1[4]-pow(obs.wdata1[3],2))/T, 
					       (obs.wdata1[6]-pow(obs.wdata1[5],2))/T,  1 / 2 - ( pow(obs.wdata1[1], 2) / (6 * obs.wdata1[2] )), obs.wdata1[7], numprocs);
	}
	else if (versus == "F") {
	write.output(myrank, h, obs.wdata1[0], (obs.wdata1[1]-pow(obs.wdata1[0],2))/T, obs.wdata1[3],  (obs.wdata1[4]-pow(obs.wdata1[3],2))/T, 
				(obs.wdata1[6]-pow(obs.wdata1[5],2))/T,  1 / 2 - ( pow(obs.wdata1[1], 2) / (6 * obs.wdata1[2] )), obs.wdata1[7], numprocs);
	}
	//}
	if (myrank == 0)
	{       
		t2 = MPI::Wtime();
		std::cout << "Total time taken by the process is about " << (t2 - t1) / 60 << " minutes approx" << std::endl;
	}
	for(int i = 0; i < M; ++i)
		{
		    for(int j = 0; j < Ns; ++j)
		    {
		        delete[] lattice[i][j];
		    }
		    delete[] lattice[i];
		
		}
	delete[] lattice;
	delete[] plaquette;
	MPI::Finalize();
	return 0;
}





















	/*
	int n_anneal=100+1;
	J3=0;
	for (int i = 0; i < isteps/2; ++i)
	{
		//Equilibration      
		update.swendsen_flip(lattice, plaquette, T);
		update.E=update.enrg(lattice, 1);
		partemp.paralleltemp(i, update);
	if (i % ((isteps/2)/(n_anneal-1)) ==0 ) J3 += J3_temp/(n_anneal-1);			
	}
	//Turn on J3;
  	*/
  	/*
	const int lentag=0;
   	const int datatag=1;
   	if (myrank == 0) {
	
   	     // Send length, then data
   	     int len = M*Ns*Nl;
   	     for (int i=1; i<numprocs; ++i){
   	     MPI_Send( &len, 1, MPI::INTEGER, i, lentag, MPI::COMM_WORLD );
   	     }
	     
             std::vector<int> send_vec;   	     
   	     for (int xc = 0; xc < M; ++xc) {
  	      for (int yc = 0; yc < Ns; ++yc) {
  	       for (int zc = 0; zc < Nl; ++zc) { 
 		   send_vec.push_back(lattice[xc][yc][zc].S());  	              
  	       }
  	      }
  	     }  
   	        
   	     for (int i=1; i<numprocs; ++i){	
   	     MPI_Send( send_vec.data(), len, MPI::INTEGER, i, datatag, MPI::COMM_WORLD );
   	     }
	     send_vec.clear();
	     
   	 } else {
	     MPI_Status status;		
   	     int len;
   	     MPI_Recv( &len, 1, MPI::INTEGER, 0, lentag, MPI::COMM_WORLD, &status);
	     
   	     int *recv_data =  new int[len];
   	     MPI_Recv( recv_data, len, MPI::INTEGER, 0, datatag, MPI::COMM_WORLD, &status);
	
	     int j=0;
   	     for (int xc = 0; xc < M; ++xc) {
  	      for (int yc = 0; yc < Ns; ++yc) {
  	       for (int zc = 0; zc < Nl; ++zc) { 
 		   lattice[xc][yc][zc].set_S(recv_data[j]);  	
 		   j+=1;              
  	       }
  	      }
  	     } 
  	     delete[] recv_data;  
  	       	         
   	 }  
   	*/ 
