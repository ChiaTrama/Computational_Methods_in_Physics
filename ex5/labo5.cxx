#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <math.h>

// funzione gauss 
 double gauss (double x, double x0, double sigma){
        double g=0.;
        g=(1./(sigma*sqrt(2.*M_PI)))*exp(-pow((x-x0)/sigma,2.)/2.);
        return g;}



using namespace std;
int main(int argc, const char * argv[]){
  
// inserisco il numero di atomi e il lato della cella di simulazione dal file argonlast.dat
    int n;
    double L;
    ifstream fin;
   
    fin.open ("argonlast.dat", std::ios::in);
    fin >> n;
    fin >> L;

// chiedo il raggio massimo, sigma e il numero di snapshot da fare
    double rmax, sigma; 
    int ns;
   
  cout << " inserisci rmax (in unità di Angstrom): \n";
  cin >> rmax;
  cout << "inserisci il parametro di broadening (in unità di Angstrom):\n";
  cin >> sigma;
  cout <<" inserisci il numero di snapshot che vuoi fare(numero intero!): \n";
  cin >> ns;

// ora leggo il file argonlast e metto le posizioni in una matrice che chiamo matrice r
    double r[1000][3];
    for (int i=0;i<n;i++){
      for (int j=0;j<3;j++){ 
          fin >> r[i][j];
      }
    }

   fin.close();
 
// faccio una griglia per il grafico
   double *x,*g;
   g=new double[ns];
   x= new double [ns];
   for (int i=0;i<ns;i++){
      g[i]=0.;
   }
   for (int i=0;i<ns;i++){ 
      x[i]=rmax/((double)ns)*((double)i);
   }


// trovo la distanza minima come in lj.cxx
   double dist, dist_min;
   
   for (int i=0; i<n; i++){

      for (int j=0; j<n; j++){

         if (j!=i){

          dist_min=1e10;
          for(int kx=-1;kx<2;kx++){
               for(int ky=-1;ky<2;ky++){
                    for(int kz=-1;kz<2;kz++){
                        dist=pow(r[i][0]-r[j][0]+L*kx,2.)+pow(r[i][1]-r[j][1]+L*ky,2.)+pow(r[i][2]-r[j][2]+L*kz,2.);
                        dist=sqrt(dist);
                if (dist<dist_min) dist_min=dist;
                                              }  
                                        } 
                                   }    
//e qui chiudo l'if, allora per i e j fissati ora calcolo la gaussiana per ogni t con la dist_min trovata e faccio la sommatoria 
    for(int i=0;i<ns;i++) g[i]+= gauss(x[i],dist_min,sigma);   
            
                            }
                          }
                        }

// qui chiudo i for su i e j

    for (int i=0;i<ns;i++) 
     g[i]*=(1./(4.*M_PI*pow(x[i]*n,2.)*pow(n,2.)))*pow(L,3.);




// ora scrivo i risultati per plottare la funzione in un file chiamato corrfunction.txt
   ofstream fout;
   fout.open("corrfunction.txt", std::ios::out);
     for (int i=0;i<ns;i++){
       fout<< x[i] <<" " << g[i]<< '\n';
     }
   fout.close();
       

   delete []g;
   delete []x;
     
return 0;}








