#include "fn.hpp"
//#include "sw.hpp"
//#include <random>
#include <iostream>
#include <stack>

//static std::random_device rd;
//static std::mt19937 gen(rd());




static ran rann;
using namespace std;



inline int fn::PBC(int a, int b, int L){return (a+b+L) % L;}
inline int fn::PBYU(int a, int b, int L){return (a+b+L) % L;}
inline int fn::PBYD(int a, int b, int L){return (a%Np < lx ) ? a+b+L: a+b;}
inline int fn::PBXR(int a, int b, int L){return ((a+b) % L==0)? a+b-L: a+b;}
inline int fn::PBXL(int a, int b, int L){return (a % L == 0) ? a+b+L: a+b;}
/*
int fn::PBUR(int a, int b, int L){
// Along n1 unit vector.
if ((a+1) % lx == 0 && a < Ns){return (a+1)%L;}
else if (a > Ns-lx){return (a+b)%L;}
else {return (a+b+L)%L;}}
int fn::PBUL(int a, int b, int L){
if (a % lx == 0 && a > lx-1){return a-1;}
else if (a < lx && a != 0){return (L+a-lx-1);}
else if (a == 0){return (L-1);}
else {return (a+b+L)%L;}}
*/

void fn::merge_clusters(sw** lattice, int cluster_id1, int cluster_id2 )
{
  int merged_cluster_id = std::min(cluster_id1, cluster_id2);
  int max_cluster_id = std::max(cluster_id1, cluster_id2); 
  for(int xc = 0; xc<M; xc++)
  {
    for(int yc=0; yc<Ns; yc++)
    {
        //for(int zc=0; zc<ly; zc++)
        //{
        if(lattice[xc][yc].cl()== max_cluster_id) lattice[xc][yc].set_cl(merged_cluster_id); 
        //}
    }
  }
}

void fn::make_bonds(sw** lattice, int xc, int yc, double temperature)
{
  double beta = 1/temperature;
  sw central = lattice[xc][yc]; 
  int ycR=yc; int nbr=0,y=0; int n1,n2,n3,n4,n5,n6;
  stack<int> cluster;
  
  while(1) {
  y+=1;
  if (yc < Np){  
  n1=PBYU(yc,lx,Np);
  n2=PBYD(yc+Np,-lx,Np);
  n3=PBXR(yc+Np,1,lx);
  n4=PBYD(yc,-lx,Np);
  n5=PBXR(PBYD(yc+Np,-lx,Np),1,lx);
  n6=yc+Np;
  
  //int G1 = lattice[xc][yc].S() * lattice[xc][n1].S() * lattice[xc][n2].S() * lattice[xc][n3].S();
  //int G2 = lattice[xc][yc].S() * lattice[xc][n4].S() * lattice[xc][n5].S() * lattice[xc][n6].S();
  } // If
  
  else {  
  n1=PBXR(yc,1,lx);
  n2=yc%Np;
  n3=PBYU(yc%Np,lx,Np);
  n4=PBXL(yc,-1,lx);
  n5=PBXL(yc%Np,-1,lx);
  n6=PBXL(PBYU(yc%Np,lx,Np),-1,lx);
 
  }
  
  int G1 = lattice[xc][yc].S() * lattice[xc][n1].S() * lattice[xc][n2].S() * lattice[xc][n3].S();
  int G2 = lattice[xc][yc].S() * lattice[xc][n4].S() * lattice[xc][n5].S() * lattice[xc][n6].S();
  
  if(rann() < central.ice_rule(lattice[xc][n1], beta, J0/M) and G1 > 0 ) 
  {
    if(!lattice[xc][n1].bonded()) 
    {
      lattice[xc][n1].set_cl(central.cl()); yc=n1; nbr=1; central = lattice[xc][yc]; cluster.push(yc);
    }
    else
    {
      merge_clusters(lattice, lattice[xc][n1].cl(), central.cl()) ;
    }
  }
    
  else if(rann() < central.ice_rule(lattice[xc][n4], beta, J0/M) and G2 > 0 and y == 1) 
  {
    if(!lattice[xc][n4].bonded()) 
    {
      lattice[xc][n4].set_cl(central.cl()); yc=n4; nbr=1; central = lattice[xc][yc];  cluster.push(yc);
    }
    else
    {
      merge_clusters(lattice, lattice[xc][n4].cl(), central.cl()) ;
    }
  }
  
  else if(rann() < central.ice_rule(lattice[xc][n6], beta, J1/M) and G1 > 0 and y == 1) 
  {
    if(!lattice[xc][n6].bonded()) 
    {
      lattice[xc][n6].set_cl(central.cl()); nbr=1;yc=n6;central = lattice[xc][yc]; cluster.push(yc);
    }
    else
    {
      merge_clusters(lattice, lattice[xc][n6].cl(), central.cl()) ;
    } 
  }
    
  else if(rann() < central.ice_rule(lattice[xc][n3], beta, J1/M) and G1 > 0) 
  {
    if(!lattice[xc][n3].bonded()) 
    {
      lattice[xc][n3].set_cl(central.cl()); nbr=1;yc=n3; central = lattice[xc][yc]; cluster.push(yc);
    }
    else
    {
      merge_clusters(lattice, lattice[xc][n3].cl(), central.cl()) ;   
    } 
  }
  // Spatial clusters along y   
  else if(rann() < central.ice_rule(lattice[xc][n5], beta, J1/M) and G2 > 0 and y == 1) 
  {
    if(!lattice[xc][n5].bonded()) 
    {
      lattice[xc][n5].set_cl(central.cl()); nbr=1;yc=n5;central = lattice[xc][yc];  cluster.push(yc);
    }
    else
    {
      merge_clusters(lattice, lattice[xc][n5].cl(), central.cl()) ;
    }
  }
    
  else if(rann() < central.ice_rule(lattice[xc][n2], beta, J1/M) and G2 > 0) 
  {
    if(!lattice[xc][n2].bonded()) 
    {
      lattice[xc][n2].set_cl(central.cl()); nbr=1;yc=n2;central = lattice[xc][yc]; cluster.push(yc);
    }
    else
    {
      merge_clusters(lattice, lattice[xc][n2].cl(), central.cl()) ;
    } 
  
  //**************************************************************************************************************************************
 
  if (ycR==yc or nbr==0){break;}}
  } // while loop
}


int fn::make_cluster(sw** lattice, double temperature)
{
  int cluster_id = 0;
  for(int xc=0; xc<M; xc++)
  {      		
    for(int ycR=0; ycR<Ns; ycR++)
    {
        //for(int zc=0; zc<ly; zc++)
        //{
        
           int yc = rann()*Ns;
           if(!lattice[xc][yc].bonded()) 
           {
             lattice[xc][yc].set_cl(cluster_id);
           }
	   // Spatial update
           make_bonds(lattice, xc, yc, temperature);
      
           if(lattice[xc][yc].cl()==cluster_id)
           {
             cluster_id++;
           }
           //std::cout << yc << "  " << lattice[xc][yc].cl() << std::endl;
        //}
   }
  }
           //std::cout <<  "  " << std::endl;
  return cluster_id;
}


void fn::flip_cluster(sw** lattice, bool* decision_list)
{
  for(int xc=0; xc<M; xc++)
  {
    for(int yc=0; yc<Ns; yc++)
    {
        //for(int zc=0; zc<ly; zc++)
        //{
          if(decision_list[lattice[xc][yc].cl()]) lattice[xc][yc].flip();
        //}
    }
  }  
}

void fn::swendsen_flip(sw** lattice, double temperature)
{

  /*int no_clusters = make_cluster(lattice, temperature);
  bool* flip_decision = new bool [no_clusters];
  for(int i=0; i<no_clusters; i++) flip_decision[i] = (rann()<=0.5);
  flip_cluster(lattice, flip_decision);
  delete[] flip_decision;
  */
  double dE,p;
  
  for(int xc=0; xc<M; xc++)
  {
    for(int yc=0; yc<Ns; yc++)
    {pl(lattice, xc, yc, 0, & p);
        //for(int zc=0; zc<ly; zc++)
        //{
           // Imaginary-time update
          //std::cout << "vikas" << p << std::endl; 
          
	  if (yc < Np){
	   dE=2*lattice[xc][yc].S()*(
	  				  -(J0/M)*(lattice[xc][PBYU(yc,lx,Np)].S() * lattice[xc][yc+Np].S() * lattice[xc][PBXR(yc+Np,1,lx)].S()+
	  				       lattice[xc][PBYD(yc,-lx,Np)].S() * lattice[xc][PBYD(yc+Np,-lx,Np)].S() * lattice[xc][PBXR(PBYD(yc+Np,-lx,Np),1,lx)].S())
	  	   			  +(J1/M)*(lattice[xc][PBYU(yc,lx,Np)].S() + lattice[xc][PBYD(yc,-lx,Np)].S())
	  	   			  +(J2/M)*p
	  	   			  );
	  	      }		
	  	      
	  else 	{
	   dE=2*lattice[xc][yc].S()*(
	  				  -(J0/M)*(lattice[xc][PBXR(yc,1,lx)].S() * lattice[xc][yc%Np].S() * lattice[xc][PBYU(yc%Np,lx,Np)].S()+
	  				       lattice[xc][PBXL(yc,-1,lx)].S() * lattice[xc][PBXL(yc%Np,-1,lx)].S() * lattice[xc][PBXL(PBYU(yc%Np,lx,Np),-1,lx)].S())
	  	   			  +(J1/M)*(lattice[xc][PBXR(yc,1,lx)].S() + lattice[xc][PBXL(yc,-1,lx)].S())
	  	   			  +(J2/M)*p
	  	   			  );
	  	      }	      	  
	  if (exp(dE/temperature) > rann()){lattice[xc][yc].flip();}
          lattice[xc][yc].reset_bond();
          
        //}
     }   
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////


double fn::magnetization(sw** lattice)
{
  double mag = 0.0;
  for(int xc=0; xc<M; xc++)
  {
    for(int yc=0; yc<Np; yc++)
    {
        //for(int zc=0; zc<ly; zc++)
        //{
        //std::cout << yc << " " << lattice[xc][yc].S() << std::endl;
        //if (pow(-1,((yc%lx)+(yc)/lx)) == 1){
          mag += double(lattice[xc][yc].S() * lattice[xc][PBYU(yc,lx,Ns)].S() * lattice[xc][yc+Np].S() * lattice[xc][PBXR(yc+Np,1,lx)].S());//*pow(-1,((yc%lx)+(yc)/lx));}
         // mag += double(lattice[xc][yc].S())*pow(-1,((yc%lx)+(yc)/lx));
     //     std::cout << "no of clusters is " << double(lattice[xc][yc].S()*lattice[xc][PBXR(yc,1,lx)].S()*lattice[xc][PBYU(yc,lx,Ns)].S()*lattice[xc][PBYU(yc,lx+1,Ns)].S()) << "  "<< std::endl; 
        //}
   // std::cout << PBYU(yc,lx,Np) <<  "  " << yc+Np << "  " << PBXR(yc+Np,1,lx) << "  " << PBYD(yc,-lx,Np) << "  " << PBYD(yc+Np,-lx,Np) << "  " << PBXR(PBYD(yc+Np,-lx,Np),1,lx) << std::endl;    
    }    
  }
    	/*for(int xc=0; xc<M; xc++){
  		  for(int yc=0; yc<Ns; yc++){
		          std::cout << yc << "  " << PBXR(yc,1,lx) << "  " << PBYU(yc,lx,Ns) << "  " << PBYU(yc,lx+1,Ns) << " " << std::endl;
		  }}*/
  return abs(2*mag)/(double(M*Ns));
}
//***************************************************************************************


void
fn::phys_quant(sw** lattice, double * op, double * op_sq, double * op_four, double * Ic, double * Ic_sq) {

    // Order parameter P,rhp, and etc, to detect antiferroelectric phase.
    double k_x = 0.f, k_y = 3.14159265359f;
    double ms1 = 0.f, ms2 = 0.f, mc1 = 0.f, mc2 = 0.f, cs = 1.f;
    * op = 0.f;* Ic = 0.f;
    
    for (int xc = 0; xc < M; ++xc) {
	for (int yc = 0; yc < Ns; ++yc) {
	    //for (int zc = 0; zc < ly; ++zc) {
	        // int i = get_index(yc, zc);
	        //dip(lattice, xc, yc, 0, & P1, & P2);

	        // *op += lattice[xc][yc].S()*lattice[xc][PBYU(yc,lx,Np)].S() * lattice[xc][yc+Np].S() * lattice[xc][PBXR(yc+Np,1,lx)].S();
	       	// if (yc >= Np){
    		// std::cout << PBXR(yc,1,lx) << "  " << yc%Np << " " << PBYU(yc%Np,lx,Np) << " " << PBXL(yc,-1,lx) << " " << PBXL(yc%Np,-1,lx) << " " << PBXL(PBYU(yc%Np,lx,Np),-1,lx) <<  std::endl;   } 
                //!------------------ stripe order-----------------------!
                double w = cs * ((k_x) * psite[yc][0] + (k_y) * psite[yc][1]);
                ms1 += lattice[xc][yc].S() * sin(w);
                mc1 += lattice[xc][yc].S() * cos(w);
                w = cs * ((k_y) * psite[yc][0] + (k_x) * psite[yc][1]);
                ms2 += lattice[xc][yc].S() * sin(w);
                mc2 += lattice[xc][yc].S() * cos(w);
                //!------------------------------------------------------!
                int v1 = lattice[xc][yc].S() - lattice[xc][PBC(yc,lx,Np)].S() + lattice[xc][yc+Np].S() - lattice[xc][PBC(yc,Np+1,lx)].S();
                int v2 = lattice[xc][yc].S() - lattice[xc][PBC(yc,lx,Np)].S() - lattice[xc][yc+Np].S() + lattice[xc][PBC(yc,Np+1,lx)].S();
                if (abs(v1) == 4) {
                    if (yc < lx*ly) {
                        * Ic = * Ic + 1;
                    }
                } else if (abs(v2) == 4) {
                    if (yc < lx*ly) {
                        * Ic = * Ic + 1;
                    }
                } else if (abs(v1) != 4 or abs(v2) != 4) {
                    if (yc < lx*ly) {
                        * Ic = * Ic - 1 / 3;
                    }
                }
            //}
        }
    }

    * op = pow(pow(ms1, 2) + pow(mc1, 2), 2);
    * op += pow(pow(ms2, 2) + pow(mc2, 2), 2);
    * op = pow( * op, 0.5);
    * op = * op / (pow(Ns * 1 * M, 2));
    * op = * op / ( 1 * 1 );
    * Ic = * Ic / (lx * ly * 1 * M);
/*
    // Correlation function in imaginary-time, to detect deconfined phase.
    if (J2 == 0) {
        for (int xc = 10; xc < 32; ++xc) 
        {
            for (int yc = 0; yc < lx*ly; ++yc) 
            {
		//for (int zc = 0; zc < ly; ++zc) 
		//{
                     double P1 = lattice[0][yc].S();//-(1/4)*(lattice[0][yc].S() - lattice[0][PBC(yc,1,lx)].S() + lattice[0][yc+lx].S() - lattice[0][yc+lx][PBC(zc,1,lx)].S());
                     double P2 = lattice[xc][yc].S();// -(1/4)*(lattice[xc][yc].S() - lattice[xc][PBC(yc,1,lx)].S() + lattice[xc][yc+lx].S() - lattice[xc][yc+lx][PBC(zc,1,lx)].S());
                     c[xc - 10] += P1 * P2 / (lx*ly);
                //}
            }
        }
              // std::cout << "vikas" <<std::endl;     
    * enrg = * enrg / (2 * lx * ly * Nl);
    * enrg_sq = pow( * enrg, 2);
    * m1 = abs( * m1 / (M * 2 * lx * ly * Nl));
    * m1_sq = pow( * m1, 2);
    */
    * op_sq = pow( * op, 2);
    * op_four = pow( * op, 4);
    * Ic_sq = ( * Ic) * ( * Ic);
    
    
}

void
fn::pl(sw** lattice, int xc, int yc, int zc, double * pl_t) {

        int l, k, dp[8];
	double Px[8], Py[8];
	if (yc<Np){
        l = yc;//bsp[0][i];
        k = PBYD(l,-lx,Np);
        
                dp[0] = PBXR(l,1,lx);//bpn[0][l];
                dp[1] = PBYU(l,lx,Np);//bpn[1][l];
                dp[2] = PBXL(l,-1,lx);//bpn[2][l];
                dp[3] = PBYD(l,-lx,Np);//bpn[3][l];
                dp[4] = PBXR(k,1,lx);
                dp[5] = PBYU(k,lx,Np);
                dp[6] = PBXL(k,-1,lx);
                dp[7] = PBYD(k,-lx,Np);

		for (int b = 0; b < 8; b++){dip(lattice, xc, dp[b], zc, & Px[b], & Py[b]);}
		 
                * pl_t =  1*( Px[0] + Px[1] + Px[2] + Px[3] ) + 1*( Py[0] + Py[1] + Py[2] + Py[3] ) +
                          -1*( Px[4] + Px[5] + Px[6] + Px[7] ) - 1*( Py[4] + Py[5] + Py[6] + Py[7] );
        	* pl_t /= 4;//pow(-1, ((yc%Np) % lx) + ((yc%Np) / lx)) 
        }//PBC()//bsp[1][i];
        else {
        l = yc%Np;//bsp[0][i];
        k = PBXL(l,-1,lx);
        
                dp[0] = PBYD(l,-lx,Np);//bpn[0][l];
                dp[1] = PBXR(l,1,lx); //bpn[1][l];
                dp[2] = PBYU(l,lx,Np);//bpn[2][l];
                dp[3] = PBXL(l,-1,lx);  //bpn[3][l];
                dp[4] = PBYD(k,-lx,Np);
                dp[5] = PBXR(k,1,lx); 
                dp[6] = PBYU(k,lx,Np);
                dp[7] = PBXL(k,-1,lx);
                //std::cout << dp[4] << " " << dp[5] << " " << dp[6] << " " << dp[7] << " " << std::endl;
		for (int b = 0; b < 8; b++){dip(lattice, xc, dp[b], zc, & Px[b], & Py[b]);}
		 
                * pl_t =  1*( Px[0] + Px[1] + Px[2] + Px[3] ) - 1*( Py[0] + Py[1] + Py[2] + Py[3] ) +
                          -1*( Px[4] + Px[5] + Px[6] + Px[7] ) + 1*( Py[4] + Py[5] + Py[6] + Py[7] );  
                * pl_t /= 4;                 
        	//* pl_t *= pow(-1, ((yc%Np) % lx) + ((yc%Np) / lx))         
        }//PBC()//bsp[1][i];            
}



void
fn::dip(sw** lattice, int xc, int yc, int zc, double * px, double * py) {
 * px = 0.25 * //pow(-1, (yc % lx) + (yc / lx)) * 
       (lattice[xc][yc + (zc * Ns)].S() + lattice[xc][yc + Np + (zc * Ns)].S() - lattice[xc][PBXR(yc+Np,1,lx)].S() - lattice[xc][PBYU(yc,lx,Np)].S());
 * py = 0.25 * //pow(-1, (yc % lx) + (yc / lx)) * 
       (lattice[xc][yc + (zc * Ns)].S() + lattice[xc][PBXR(yc+Np,1,lx)].S() - lattice[xc][yc + Np + (zc * Ns)].S() - lattice[xc][PBYU(yc,lx,Np)].S());
 	        //std::cout << yc << " " << * px << " " << * py << std::endl;      
}
/*

void fn::make_bonds(sw** lattice, int xc, int yc, double temperature)
{
  double beta = 1/temperature;
  sw central = lattice[xc][yc];

  if(rann() < central.bond_prob(lattice[PBC(xc,1,M)][yc], beta, Kx)) 
  {
    if(!lattice[PBC(xc,1,M)][yc].bonded()) 
    {
      lattice[PBC(xc,1,M)][yc].set_cl(central.cl());
    }
    else
    {
      merge_clusters(lattice, lattice[PBC(xc,1,M)][yc].cl(), central.cl()) ;
    }
  }
  
  if(rann() < central.bond_prob(lattice[PBC(xc,-1,M)][yc], beta, Kx)) 
  {
    if(!lattice[PBC(xc,-1,M)][yc].bonded()) 
    {
      lattice[PBC(xc,-1,M)][yc].set_cl(central.cl());
    }
    else
    {
      merge_clusters(lattice, lattice[PBC(xc,-1,M)][yc].cl(), central.cl()) ;
    }
  }
  // Spatial clusters along x 
  if(rann() < ice_rule(lattice[xc][PBXR(yc,Np+1,lx)], beta, J0/M)) 
  {
    if(!lattice[xc][PBXR(yc,Np+1,lx)].bonded()) 
    {
      lattice[xc][PBXR(yc,Np+1,lx)].set_cl(central.cl());
    }
    else
    {
      merge_clusters(lattice, lattice[xc][PBXR(yc,Np+1,lx)].cl(), central.cl()) ;
    }
  }
    
  if(rann() < ice_rule(lattice[xc][PBXL(yc,-(Np+1),lx)], beta, J0/M)) 
  {
    if(!lattice[xc][PBXL(yc,-(Np+1),lx)].bonded()) 
    {
      lattice[xc][PBXL(yc,-(Np+1),lx)].set_cl(central.cl());
    }
    else
    {
      merge_clusters(lattice, lattice[xc][PBXL(yc,-(Np+1),lx)].cl(), central.cl()) ;
    }
  }
  
  if(rann() < ice_rule(lattice[xc][PBYU(yc,Np+2,Ns)], beta, J0/M)) 
  {
    if(!lattice[xc][PBYU(yc,Np+2,Ns)].bonded()) 
    {
      lattice[xc][PBYU(yc,Np+2,Ns)].set_cl(central.cl());
    }
    else
    {
      merge_clusters(lattice, lattice[xc][PBYU(yc,Np+2,Ns)].cl(), central.cl()) ;
    }
  }
    
  if(rann() < ice_rule(lattice[xc][PBYD(yc,-(Np+2),Ns)], beta, J0/M)) 
  {
    if(!lattice[xc][PBYD(yc,-(Np+2),Ns)].bonded()) 
    {
      lattice[xc][PBYD(yc,-(Np+2),Ns)].set_cl(central.cl());
    }
    else
    {
      merge_clusters(lattice, lattice[xc][PBYD(yc,-(Np+2),Ns)].cl(), central.cl()) ;
    }
  }
  
}

*/

