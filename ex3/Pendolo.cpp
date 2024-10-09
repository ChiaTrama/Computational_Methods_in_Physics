#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <vector>


int main(int argc, const char * argv[])
{ 
    using namespace std;
    double l=1.;
    string filename;
    ofstream fout;
    double g=9.806;
    int i=0; 
    int n=0.;
    double h=0.;
    int a =0.;
    int z=0.;

    cout << "Inserici nome file su cui vuoi trascrivere \n";
    cin >> filename;


    // inseriesco le condizioni iniziali
    double theta0=0.;  
    double dtheta0=0.;  

    cout << "Inserici l'angolo iniziale theta0 \n";
    cin >> theta0;  
    cout<<"inserisci  la spaziatura tra i theta";
    cin>> h;

    // calcolo  il numero di step
     n = (theta0)/h;

    
    double *Y1,*Y2,*Y3,*Y4,*Y;
    double *W1,*W2,*W3,*W4,*W;
    double *y,*Tvero,*Tapprox; 
  
    Y1 = new double[2];
    Y2 = new double[2];
    Y3 = new double[2];
    Y4 = new double[2];
    Y  = new double[2];
    W1 = new double[2];
    W2 = new double[2];
    W3 = new double[2];
    W4 = new double[2];
    W  = new double[2];
    y = new double[n];
    Tvero = new double[n];
    Tapprox = new double[n];
 
 
   for (int i=0;i<n;i++) {
        y[i]=theta0-i*h;}

    for (int i=0;i<n;i++) {

  //setto le condizioni iniziali 
        z=0;
        Y1[0]=y[i];
        Y1[1]=0.;

    do{
            Y2[0]=Y1[0]+(h/2.)*Y1[1];
    	    Y2[1]=Y1[1]-(g/l)*(h/2.)*sin(Y1[0]);
    	
    	    Y3[0]=Y1[0]+(h/2.)*(Y1[1]-(g/l)*(h/2.)*sin(Y1[0]));
    	    Y3[1]=Y1[1]-(g/l)*(h/2.)*sin(Y2[0]);
    	
    	    Y4[0]=Y1[0]+(h)*(Y1[1]-(g/l)*(h/2.)*sin(Y2[0]));
    	    Y4[1]=Y1[1]-(g/l)*(h)*sin(Y3[0]);
    	
    	    Y[0]=Y1[0]+(Y1[1]+2.*(Y1[1]-(g/l)*(h/2.)*sin(Y1[0]))+2.*(Y1[1]-(g/l)*(h/2.)*sin(Y2[0]))+(Y1[1]-(g/l)*(h)*sin(Y3[0])))*(h/6.);
    	    Y[1]=Y1[1]+(-g/l)*(sin(Y1[0]) + 2.*sin(Y2[0]) + 2.*sin(Y3[0]) + sin(Y4[0]))*(h/6.);
          if ((Y1[0]*Y[0])>=0.) {
    	        Y1[0]=Y[0];
    	        Y1[1]=Y[1];
    	        Y2[0]=0  ;
    	        Y2[1]=0.;
    	        Y3[0]=0.;
    	        Y3[1]=0.;
    	        Y4[0]=0.;
    	        Y4[1]=0.;
    	        Y[0]=0.;
    	        Y[1]=0.;
    	        z++;
    	        continue;
    	    } else {   
                    break;}
    	}  while ((Y1[0]*Y[0])>=0.);
  

 //interpolazione per trovare periodo(un quarto)
 Tvero[i]=(z*h +((h/(Y1[0]-Y[0]))*Y1[0]));
}

 // mi calcolo il periodo con seno approssimat0

    for (int i=0;i<n;i++) {
        a=0;
        W1[0]=y[i];
        W1[1]=0.;
        do {
            W2[0]=W1[0]+(h/2.)*W1[1];
            W2[1]=W1[1]-(g/l)*(h/2.)*(W1[0]);

            W3[0]=W1[0]+(h/2.)*(W1[1]-(g/l)*(h/2.)*(W1[0]));
            W3[1]=W1[1]-(g/l)*(h/2.)*(W2[0]);

            W4[0]=W1[0]+(h)*(W1[1]-(g/l)*(h/2.)*(W2[0]));
            W4[1]=W1[1]-(g/l)*(h)*(W3[0]);

            W[0]=W1[0]+(W1[1]+2.*(W1[1]-(g/l)*(h/2.)*(W1[0]))+2.*(W1[1]-(g/l)*(h/2.)*(W2[0]))+(W1[1]-(g/l)*(h)*(W3[0])))*(h/6.);
            W[1]=W1[1]+(-g/l)*((W1[0])+2.*(W2[0])+2.*(W3[0])+(W4[0]))*(h/6.);
    
            if ((W1[0]*W[0])>=0.) {
                W1[0]=W[0];
                W1[1]=W[1];
                W2[0]=0.;
                W2[1]=0.;
                W3[0]=0.;
                W3[1]=0.;
                W4[0]=0.;
                W4[1]=0.;
                W[0]=0.;
                W[1]=0.;
                a++;
                continue;
            } else {
                break;
  	        }
        } while ((W1[0]*W[0])>=0.);

 // Calcolo un quarto di perido
       Tapprox[i]=(a*h+((h/(W1[0]-W[0]))*W1[0]));
    }
    

 fout.open(filename,ios::out);
 fout.precision(10);
 for (i=0;i<n;i++) {
        fout << y[i] << " " << Tvero[i]/Tapprox[i] << '\n';}
 fout.close();

    delete [] Y1;
    delete [] Y2;
    delete [] Y3;
    delete [] Y4;
    delete [] Y;
    delete [] W1; 
    delete [] W2;
    delete [] W3;
    delete [] W4;
    delete [] W;
    delete [] Tvero;
    delete [] Tapprox;
    delete [] y;


    return 0;
}

