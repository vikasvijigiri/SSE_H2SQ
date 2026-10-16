#include "fn.hpp"
#include <iostream>
#include <deque>


static ran rann;
using namespace std;

inline int fn::PBC(int a, int b, int L) {
  return (a + b + L) % L;
}
inline int fn::PBYU(int a, int b, int L) {
  return (a + b + L) % L;
}
inline int fn::PBYD(int a, int b, int L) {
  return (a % Np < lx) ? a + b + L : a + b;
}
inline int fn::PBXR(int a, int b, int L) {
  return ((a + b) % L == 0) ? a + b - L : a + b;
}
inline int fn::PBXL(int a, int b, int L) {
  return (a % L == 0) ? a + b + L : a + b;
}

void fn::merge_clusters_intl(sw ** * lattice, int cluster_id1, int cluster_id2) {
  int merged_cluster_id = std::min(cluster_id1, cluster_id2);
  int max_cluster_id = std::max(cluster_id1, cluster_id2);
  for (int xc = 0; xc < M; xc++) {
    for (int yc = 0; yc < Ns; yc++) {
      for (int zc = 0; zc < Nl; zc++) {
        if (lattice[xc][yc][zc].cl() == max_cluster_id) lattice[xc][yc][zc].set_cl(merged_cluster_id);
      }
    }
  }
}

void fn::make_bonds(sw ** * lattice, sw * plaquette, deque < int > * cluster, deque < int > * plaq_cluster, int xc, int yc, int zc, double temperature) {
  double beta = 1 / temperature;
  sw central = lattice[xc][yc][zc];
  sw plaq = plaquette[yc % Np];
  int ycR = yc;
  int n1, n2, n3, n4, n5, n6;

  for (int i = 0; i < 2; ++i) {
    yc = ycR;
    while (1) {
      int nbr = 0;
      //std::cout << "clusters   " << yc << std::endl;
      if (yc < Np) {
        n1 = PBYU(yc, lx, Np);
        n2 = yc + Np;
        n3 = PBXR(yc + Np, 1, lx);
        n4 = PBYD(yc, -lx, Np);
        n5 = PBXR(PBYD(yc + Np, -lx, Np), 1, lx);
        n6 = yc + Np;
      } // If
      else {
        n1 = PBXR(yc, 1, lx);
        n2 = yc % Np;
        n3 = PBYU(yc % Np, lx, Np);
        n4 = PBXL(yc, -1, lx);
        n5 = PBXL(yc % Np, -1, lx);
        n6 = PBXL(PBYU(yc % Np, lx, Np), -1, lx);
      } // else

      int G1 = lattice[xc][yc][zc].S() * lattice[xc][n1][zc].S() * lattice[xc][n2][zc].S() * lattice[xc][n3][zc].S();
      int G2 = lattice[xc][yc][zc].S() * lattice[xc][n4][zc].S() * lattice[xc][n5][zc].S() * lattice[xc][n6][zc].S();
      if (rann() < central.ice_rule(lattice[xc][n1][zc], beta, J0 / M) and G1 > 0) {
        if (!lattice[xc][n1][zc].bonded() and!plaquette[yc % Np].bonded()) {
          lattice[xc][n1][zc].set_cl(central.cl());
          plaquette[yc % Np].set_cl(plaq.cl());
          plaq = plaquette[yc % Np];
          nbr = 1;
          if (n1 != ycR) cluster -> push_back(n1);
          plaq_cluster -> push_back(yc % Np);
          yc = n1;
          central = lattice[xc][yc][zc];
        }
      } else if (rann() < central.ice_rule(lattice[xc][n4][zc], beta, J0 / M) and G2 > 0) {
        if (!lattice[xc][n4][zc].bonded() and!plaquette[n4 % Np].bonded()) {
          lattice[xc][n4][zc].set_cl(central.cl());
          plaquette[n4 % Np].set_cl(plaq.cl());
          plaq = plaquette[n4 % Np];
          nbr = 1;
          if (n4 != ycR) cluster -> push_back(n4);
          plaq_cluster -> push_back(n4 % Np);
          yc = n4;
          central = lattice[xc][yc][zc];
        }
      } else if (rann() < central.ice_rule(lattice[xc][n6][zc], beta, J1 / M) and G1 > 0) {
        if (!lattice[xc][n6][zc].bonded() and!plaquette[n4 % Np].bonded()) {
          lattice[xc][n6][zc].set_cl(central.cl());
          plaquette[n4 % Np].set_cl(plaq.cl());
          plaq = plaquette[n4 % Np];
          nbr = 1;
          if (n6 != ycR) cluster -> push_back(n6);
          plaq_cluster -> push_back(n4 % Np);
          yc = n6;
          central = lattice[xc][yc][zc];
        }
      } else if (rann() < central.ice_rule(lattice[xc][n3][zc], beta, J1 / M) and G1 > 0) {
        if (!lattice[xc][n3][zc].bonded() and!plaquette[yc % Np].bonded()) {
          lattice[xc][n3][zc].set_cl(central.cl());
          plaquette[yc % Np].set_cl(plaq.cl());
          plaq = plaquette[yc % Np];
          nbr = 1;
          if (n3 != ycR) cluster -> push_back(n3);
          plaq_cluster -> push_back(yc % Np);
          yc = n3;
          central = lattice[xc][yc][zc];
        }
      }
      // Spatial clusters along y   
      else if (rann() < central.ice_rule(lattice[xc][n5][zc], beta, J1 / M) and G2 > 0) {
        if (!lattice[xc][n5][zc].bonded() and!plaquette[n4 % Np].bonded()) {
          lattice[xc][n5][zc].set_cl(central.cl());
          plaquette[n4 % Np].set_cl(plaq.cl());
          plaq = plaquette[n4 % Np];
          nbr = 1;
          if (n5 != ycR) cluster -> push_back(n5);
          plaq_cluster -> push_back(n4 % Np);
          yc = n5;
          central = lattice[xc][yc][zc];
        }
      } else if (rann() < central.ice_rule(lattice[xc][n2][zc], beta, J1 / M) and G2 > 0) {
        if (!lattice[xc][n2][zc].bonded() and!plaquette[yc % Np].bonded()) {
          lattice[xc][n2][zc].set_cl(central.cl());
          plaquette[yc % Np].set_cl(plaq.cl());
          plaq = plaquette[yc % Np];
          nbr = 1;
          if (n2 != ycR) cluster -> push_back(n2);
          plaq_cluster -> push_back(yc % Np);
          yc = n2;
          central = lattice[xc][yc][zc];
        }
      }
      //**************************************************************************************************************************************

      if (ycR == yc or nbr == 0) {
        break;
      }
    }
  } // while loop
}

void fn::flip_cluster(sw*** lattice, deque < int > cluster, bool decision, int xc, int zc) {
  if (decision) {
    for (std::deque < int > ::iterator it = cluster.begin(); it != cluster.end(); ++it) {
      lattice[xc][ * it][zc].flip();
      lattice[xc][ * it][zc].reset_bond(); //Here bond is like site.
    }
  }
}

void fn::make_bonds_intl(sw*** lattice, int xc, int yc, int zc, double temperature) {
  double beta = 1 / temperature;
  sw central = lattice[xc][yc][zc];

  if (rann() < central.bond_prob(lattice[xc][yc][PBC(zc, 1, Nl)], beta, J3 / M)) {
    if (!lattice[xc][yc][PBC(zc, 1, Nl)].bonded()) {
      lattice[xc][yc][PBC(zc, 1, Nl)].set_cl(central.cl());
    } else {
      merge_clusters_intl(lattice, lattice[xc][yc][PBC(zc, 1, Nl)].cl(), central.cl());
    }
  }

  if (rann() < central.bond_prob(lattice[xc][yc][PBC(zc, -1, Nl)], beta, J3 / M)) {
    if (!lattice[xc][yc][PBC(zc, -1, Nl)].bonded()) {
      lattice[xc][yc][PBC(zc, -1, Nl)].set_cl(central.cl());
    } else {
      merge_clusters_intl(lattice, lattice[xc][yc][PBC(zc, -1, Nl)].cl(), central.cl());
    }
  }
  
  if (rann() < central.bond_prob(lattice[PBC(xc,1,M)][yc][zc], beta, Kx)) {
    if (!lattice[PBC(xc,1,M)][yc][zc].bonded()) {
      lattice[PBC(xc,1,M)][yc][zc].set_cl(central.cl());
    } else {
      merge_clusters_intl(lattice, lattice[PBC(xc,1,M)][yc][zc].cl(), central.cl());
    }
  }

  if (rann() < central.bond_prob(lattice[PBC(xc,-1,M)][yc][zc], beta, Kx)) {
    if (!lattice[PBC(xc,-1,M)][yc][zc].bonded()) {
      lattice[PBC(xc,-1,M)][yc][zc].set_cl(central.cl());
    } else {
      merge_clusters_intl(lattice, lattice[PBC(xc,-1,M)][yc][zc].cl(), central.cl());
    }
  } 
}

int fn::make_cluster_intl(sw*** lattice, int yc, double temperature) {
  int cluster_id = 0;
  for (int xc = 0; xc < M; xc++) {
  for (int zc = 0; zc < Nl; zc++) {
    if (!lattice[xc][yc][zc].bonded()) {
      lattice[xc][yc][zc].set_cl(cluster_id);
    }
    make_bonds_intl(lattice, xc, yc, zc, temperature);

    if (lattice[xc][yc][zc].cl() == cluster_id) {
      cluster_id++;
    }
  }}
  return cluster_id;
}

void fn::flip_cluster_intl(sw *** lattice, int yc, bool * decision_list) {
  for (int xc = 0; xc < M; xc++) {
  for (int zc = 0; zc < Nl; zc++) {
    if (decision_list[lattice[xc][yc][zc].cl()]) lattice[xc][yc][zc].flip();
  }}
}

double * fn::cluster_field_intl(sw ** * lattice, int no_clusters, int yc) {
  double * cluster_field = new double[no_clusters];
  for (int xc=0; xc < M; xc++){
  for (int zc = 0; zc < Nl; zc++) {
    cluster_field[lattice[xc][yc][zc].cl()] += enrg_diff_local(lattice, xc, yc, zc, 0);
  }}
  return cluster_field;
}

void fn::swendsen_flip(sw *** lattice, sw * plaquette, double temperature) {
  // SW Update both along the interlayers and trotter axis.
  
    for (int yc = 0; yc < Ns; ++yc) {
      int no_clusters = make_cluster_intl(lattice, yc, temperature);
      bool * flip_decision = new bool[no_clusters];
      double * cluster_field = new double[no_clusters];
      cluster_field = cluster_field_intl(lattice, no_clusters, yc);
      for (int i = 0; i < no_clusters; i++) flip_decision[i] = (rann() <= exp(cluster_field[i] / temperature));
      flip_cluster_intl(lattice, yc, flip_decision);
      delete[] flip_decision, cluster_field;
    }

  // Loop Update in 2D along with the single flip metropolis in all directions.  
  int cluster_id = 0;
  double dE;
  deque < int > cluster, plaq_cluster;

  for (int xc = 0; xc < M; ++xc) {
    for (int zcR = 0; zcR < Nl; ++zcR) {
      for (int ycR = 0; ycR < Ns; ++ycR) {
        int yc = rann() * Ns; int zc = rann() * Nl;
        //spatial-update.
        cluster.push_back(yc);
        make_bonds(lattice, plaquette, & cluster, & plaq_cluster, xc, yc, zc, temperature); //filled cluster
        plaq_cluster = dip_enrg_diff_cluster(lattice, plaquette, plaq_cluster, cluster, & dE, xc, zc);
        flip_cluster(lattice, cluster, (rann() >= exp(-dE / temperature)), xc, zc); // perform operations and empty the cluster and return back.
        cluster.clear();
        //lattice[xc][ycR][zcR].reset_bond();
        //single-flip update.
        dE = enrg_diff_local(lattice, xc, yc, zc, 1);
        if (exp(dE / temperature) > rann()) {
          lattice[xc][yc][zc].flip();
        }

      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////

double fn::enrg_diff_local(sw *** lattice, int xc, int yc, int zc, int cl) {

  double dE, p;
  pl(lattice, xc, yc, zc, & p);
  if (yc < Np) {
    dE = 2 * lattice[xc][yc][zc].S() * (
      -(J0 / M) * (lattice[xc][PBYU(yc, lx, Np)][zc].S() * lattice[xc][yc + Np][zc].S() * lattice[xc][PBXR(yc + Np, 1, lx)][zc].S() +
        lattice[xc][PBYD(yc, -lx, Np)][zc].S() * lattice[xc][PBYD(yc + Np, -lx, Np)][zc].S() * lattice[xc][PBXR(PBYD(yc + Np, -lx, Np), 1, lx)][zc].S()) +
      (J1 / M) * (lattice[xc][PBYU(yc, lx, Np)][zc].S() + lattice[xc][PBYD(yc, -lx, Np)][zc].S()) +
      (J2 / M) * p +
      cl * (J3 / M) * (lattice[xc][yc][PBC(zc, 1, Nl)].S() + lattice[xc][yc][PBC(zc, -1, Nl)].S())
    );
  } else {
    dE = 2 * lattice[xc][yc][zc].S() * (
      -(J0 / M) * (lattice[xc][PBXR(yc, 1, lx)][zc].S() * lattice[xc][yc % Np][zc].S() * lattice[xc][PBYU(yc % Np, lx, Np)][zc].S() +
        lattice[xc][PBXL(yc, -1, lx)][zc].S() * lattice[xc][PBXL(yc % Np, -1, lx)][zc].S() * lattice[xc][PBXL(PBYU(yc % Np, lx, Np), -1, lx)][zc].S()) +
      (J1 / M) * (lattice[xc][PBXR(yc, 1, lx)][zc].S() + lattice[xc][PBXL(yc, -1, lx)][zc].S()) +
      (J2 / M) * p +
      cl * (J3 / M) * (lattice[xc][yc][PBC(zc, 1, Nl)].S() + lattice[xc][yc][PBC(zc, -1, Nl)].S())
    );
  }

  return dE;
}

deque < int > fn::dip_enrg_diff_cluster(sw ** * lattice, sw * plaquette, deque < int > cluster, deque < int > clust, double * dE, int xc, int zc) {
  double E_old, E_new, E_J3, px1, py1, px[4], py[4];
  int n[4];

  E_old = 0;
  for (std::deque < int > ::iterator it = cluster.begin(); it != cluster.end(); ++it) {
    int yc = * it;
    plaquette[yc].set_S(0);
    n[0] = PBYU(yc, lx, Np);
    n[1] = PBXL(yc, -1, lx);
    n[2] = PBYD(yc, -lx, Np);
    n[3] = PBXR(yc, 1, lx);
    dip(lattice, xc, yc, zc, & px1, & py1);
    for (int i = 0; i < 4; ++i) {
      if (plaquette[n[i]].S() != 0) {
        dip(lattice, xc, n[i], zc, & px[i], & py[i]);
        E_old += (J2 / M) * (px[i] * px1 + py[i] * py1);
        plaquette[yc].reset_bond();
      }
    }
  }
  E_new = 0;
  E_J3 = 0;
  for (std::deque < int > ::iterator it = clust.begin(); it != clust.end(); ++it) {
    E_J3 += (J3 / M) * 2 * lattice[xc][ * it][zc].S() * (lattice[xc][ * it][PBC(zc, 1, Nl)].S() + lattice[xc][ * it][PBC(zc, -1, Nl)].S());
    lattice[xc][ * it][zc].flip();
    lattice[xc][ * it][zc].reset_bond();
  }

  for (std::deque < int > ::iterator it = cluster.begin(); it != cluster.end(); ++it) {
    int yc = * it;
    plaquette[yc].set_S(1);
    n[0] = PBYU(yc, lx, Np);
    n[1] = PBXL(yc, -1, lx);
    n[2] = PBYD(yc, -lx, Np);
    n[3] = PBXR(yc, 1, lx);
    dip(lattice, xc, yc, zc, & px1, & py1);
    for (int i = 0; i < 4; ++i) {
      if (plaquette[n[i]].S() != 1) {
        dip(lattice, xc, n[i], zc, & px[i], & py[i]);
        E_new += (J2 / M) * (px[i] * px1 + py[i] * py1);
      }
    }

  }
  * dE = E_new - E_old ;
  cluster.clear();
  //std::cout << *dE << std::endl; 
  return cluster;
}
//************************************ Physical quantities ***************************************************

void
fn::phys_quant(sw ** * lattice, double * op, double * op_sq, double * op_four, double * Ic, double * Ic_sq) {

  // Order parameter P,rhp, and other physical quantities.
  double k_x = 0.f, k_y = 3.14159265359f;
  double ms1 = 0.f, ms2 = 0.f, mc1 = 0.f, mc2 = 0.f, cs = 1.f;
  * op = 0.f;* Ic = 0.f;
  int zc = 0;

  for (int xc = 0; xc < M; ++xc) {
    for (int yc = 0; yc < Ns; ++yc) {
      //Energy and specific heat.
    
    
      //!------------------ stripe order P, rho (X_P, X_rho)-----------------------!
      double w = cs * ((k_x) * psite[yc][0] + (k_y) * psite[yc][1]);
      ms1 += lattice[xc][yc][zc].S() * sin(w);
      mc1 += lattice[xc][yc][zc].S() * cos(w);
      w = cs * ((k_y) * psite[yc][0] + (k_x) * psite[yc][1]);
      ms2 += lattice[xc][yc][zc].S() * sin(w);
      mc2 += lattice[xc][yc][zc].S() * cos(w);
      //!------------------------------------------------------!
      if (yc < Np) {
        int v1 = lattice[xc][yc][zc].S() - lattice[xc][PBC(yc, lx, Np)][zc].S() + lattice[xc][yc + Np][zc].S() - lattice[xc][PBC(yc, Np + 1, lx)][zc].S();
        int v2 = lattice[xc][yc][zc].S() - lattice[xc][PBC(yc, lx, Np)][zc].S() - lattice[xc][yc + Np][zc].S() + lattice[xc][PBC(yc, Np + 1, lx)][zc].S();
        if (abs(v1) == 4) {
          if (yc < lx * ly) {
            * Ic = * Ic + 1;
          }
        } else if (abs(v2) == 4) {
          if (yc < lx * ly) {
            * Ic = * Ic + 1;
          }
        } else if (abs(v1) != 4 or abs(v2) != 4) {
          if (yc < lx * ly) {
            * Ic = * Ic - 1 / 3;
          }
        }
      }
    }
  }

  * op = pow(pow(ms1, 2) + pow(mc1, 2), 2);
  * op += pow(pow(ms2, 2) + pow(mc2, 2), 2);
  * op = pow( * op, 0.5);
  * op = * op / (pow(Ns * 1 * M, 2));
  * op = * op / (1 * 1);
  * Ic = * Ic / (lx * ly * 1 * M);
  /*
      // Correlation function in imaginary-time, to detect deconfined phase.
      if (J2 == 0 and M >=32) {
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
fn::pl(sw ** * lattice, int xc, int yc, int zc, double * pl_t) {

  int l, k, dp[8];
  double Px[8], Py[8];
  if (yc < Np) {
    l = yc;
    k = PBYD(l, -lx, Np);

    dp[0] = PBXR(l, 1, lx); //bpn[0][l];
    dp[1] = PBYU(l, lx, Np); //bpn[1][l];
    dp[2] = PBXL(l, -1, lx); //bpn[2][l];
    dp[3] = PBYD(l, -lx, Np); //bpn[3][l];
    dp[4] = PBXR(k, 1, lx);
    dp[5] = PBYU(k, lx, Np);
    dp[6] = PBXL(k, -1, lx);
    dp[7] = PBYD(k, -lx, Np);

    for (int b = 0; b < 8; b++) {
      dip(lattice, xc, dp[b], zc, & Px[b], & Py[b]);
    }

    * pl_t = 1 * (Px[0] + Px[1] + Px[2] + Px[3]) + 1 * (Py[0] + Py[1] + Py[2] + Py[3]) +
      -1 * (Px[4] + Px[5] + Px[6] + Px[7]) - 1 * (Py[4] + Py[5] + Py[6] + Py[7]);
    * pl_t /= 4; //pow(-1, ((yc%Np) % lx) + ((yc%Np) / lx)) 
  } //PBC()//bsp[1][i];
  else {
    l = yc % Np; //bsp[0][i];
    k = PBXL(l, -1, lx);

    dp[0] = PBYD(l, -lx, Np); //bpn[0][l];
    dp[1] = PBXR(l, 1, lx); //bpn[1][l];
    dp[2] = PBYU(l, lx, Np); //bpn[2][l];
    dp[3] = PBXL(l, -1, lx); //bpn[3][l];
    dp[4] = PBYD(k, -lx, Np);
    dp[5] = PBXR(k, 1, lx);
    dp[6] = PBYU(k, lx, Np);
    dp[7] = PBXL(k, -1, lx);
    for (int b = 0; b < 8; b++) {
      dip(lattice, xc, dp[b], zc, & Px[b], & Py[b]);
    }

    * pl_t = 1 * (Px[0] + Px[1] + Px[2] + Px[3]) - 1 * (Py[0] + Py[1] + Py[2] + Py[3]) +
      -1 * (Px[4] + Px[5] + Px[6] + Px[7]) + 1 * (Py[4] + Py[5] + Py[6] + Py[7]);
    * pl_t /= 4;
    //* pl_t *= pow(-1, ((yc%Np) % lx) + ((yc%Np) / lx))         
  }
}

void
fn::dip(sw ** * lattice, int xc, int yc, int zc, double * px, double * py) {
  * px = 0.25 * //pow(-1, (yc % lx) + (yc / lx)) * 
    (lattice[xc][yc][zc].S() + lattice[xc][yc + Np][zc].S() - lattice[xc][PBXR(yc + Np, 1, lx)][zc].S() - lattice[xc][PBYU(yc, lx, Np)][zc].S());
  * py = 0.25 * //pow(-1, (yc % lx) + (yc / lx)) * 
    (lattice[xc][yc][zc].S() + lattice[xc][PBXR(yc + Np, 1, lx)][zc].S() - lattice[xc][yc + Np][zc].S() - lattice[xc][PBYU(yc, lx, Np)][zc].S());
}
