//
//  main.cpp
//  Ising
//
//  Created by Paolo Umari on 13/12/19.
//  Copyright (c) 2019 unipd. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;


double energia(double Jconst, int nx, int ny, double* spins)
{
    double ene=0.;
    int icount=0;
    for(int j=0;j<ny;j++){
        for( int i=0; i<nx;i++){
            
            int pv_x,pv_y,pv_count;
            
            if(i<nx-1){
                pv_x=i+1;
            }else{
                pv_x=0;
            }
            pv_y=j;
            pv_count=pv_y*nx+pv_x;
            ene+=-0.5*Jconst*spins[icount]*spins[pv_count];
            
            if(i>0){
                pv_x=i-1;
            }else{
                pv_x=nx-1;
            }
            pv_y=j;
            pv_count=pv_y*nx+pv_x;
            ene+=-0.5*Jconst*spins[icount]*spins[pv_count];
            
            
            if(j<ny-1){
                pv_y=j+1;
            }else{
                pv_y=0;
            }
            pv_x=i;
            pv_count=pv_y*nx+pv_x;
            ene+=-0.5*Jconst*spins[icount]*spins[pv_count];
            
            if(j>0){
                pv_y=j-1;
            }else{
                pv_y=ny-1;
            }
            pv_x=i;
            pv_count=pv_y*nx+pv_x;
            ene+=-0.5*Jconst*spins[icount]*spins[pv_count];

            
           
            icount++;
       }
        
    }
    
    return ene;
}


double numero_random( void){
    const long int a=16807;
    const long int c=0;
    const long mt=2147483647;
    static long int x0=1;
    long int x1;
    double r;
    
    x1=(a*x0+c)%mt;
    r=((double) x1)/((double) mt);
    x0=x1;
    
    return r;
    
}

void genera_configurazione(int nx, int ny, double *spins)
{
    double r;
    int icount=0;
    for(int j=0;j<ny;j++){
        for( int i=0; i<nx;i++){
            r =sqrt(pow(numero_random(),2));
            if(r<0.5){
                spins[icount]=-1;
            }else{
                spins[icount]=1;
            }
            icount++;
        }
    }
    
}


double cambia_configurazione(double Jconst,int nx, int ny, double *spins, double* spins1)
{
    for(int i=0;i< nx*ny;i++){
        spins1[i]=spins[i];
    }
    double r;
    int i,j;
    int icount;
    double ene=0.;
    
    r=sqrt(pow(numero_random(),2));
    icount=floor(r*(nx*ny));
  //  std::cout << " Cambia" << icount << '\n' ;
    spins1[icount]*=-1.;
    j=floor((double)icount/(double)nx);
    i=icount-j*nx;
    
    int pv_x,pv_y,pv_count;
    
    if(i<nx-1){
        pv_x=i+1;
    }else{
        pv_x=0;
    }
    pv_y=j;
    pv_count=pv_y*nx+pv_x;
    ene+=-Jconst*spins1[icount]*spins1[pv_count];
    ene-=-Jconst*spins[icount]*spins[pv_count];
    
    if(i>0){
        pv_x=i-1;
    }else{
        pv_x=nx-1;
    }
    pv_y=j;
    pv_count=pv_y*nx+pv_x;
    ene+=-Jconst*spins1[icount]*spins1[pv_count];
    ene-=-Jconst*spins[icount]*spins[pv_count];
    
    
    if(j<ny-1){
        pv_y=j+1;
    }else{
        pv_y=0;
    }
    pv_x=i;
    pv_count=pv_y*nx+pv_x;
    ene+=-Jconst*spins1[icount]*spins1[pv_count];
    ene-=-Jconst*spins[icount]*spins[pv_count];
    
    if(j>0){
        pv_y=j-1;
    }else{
        pv_y=ny-1;
    }
    pv_x=i;
    pv_count=pv_y*nx+pv_x;
    ene+=-Jconst*spins1[icount]*spins1[pv_count];
    ene-=-Jconst*spins[icount]*spins[pv_count];
    

    return ene;
    
}

double magnetizzazione(int nx, int ny, double *spins)
{   int icount=0;
    double ma=0;
    for(int j=0;j<ny;j++){
        for( int i=0; i<nx;i++){
            ma+=spins[icount];
            icount++;
        }
    }
    return ma;
}

int main(int argc, const char * argv[])
{
    int nx,ny;
    long int nsteps;
    double J,kbT;
    double ene_ave=0., ene_var=0.;
    double ma_ave=0., ma_var=0.;
    double masq_ave=0., enesq_ave=0.;
    double E=0.;
    double M=0.;
    double de=0.;
    double test=0.;
    int w=1;
    double m=0.;
    double Etot=0.;
    double Mtot=0.;
    double Etot2=0.;
    double Mtot2=0.;
    int sw=0;
    
    double ene0,ene1;
    double ratio,r;
    int metodo;
    int icount=0;
    double diff_ene;
    
    std::cout << "Nx :\n";
    std::cin >> nx;
    std::cout << "Ny :\n";
    std::cin >> ny;
    std::cout << "J :\n";
    std::cin >> J;
    std::cout << "kbT :\n";
    std::cin >> kbT;
//    std::cout << "N passi :\n";
//    std::cin >> nsteps;
    
    
    std::cout.precision(10);
    
    
//    double* en = new double[nsteps];
//    double* mag = new double[nsteps];
    double* spins = new double[nx*ny];
    double* spins1 = new double[nx*ny];
    
    ofstream fout1("energie.txt");
    ofstream fout2("magnet.txt");
    
    genera_configurazione(nx, ny, spins);
    E = energia(J, nx, ny, spins);
    M = magnetizzazione(nx,ny,spins);
    fout1 << 0 << " " << E/(nx*ny) << endl;
    fout2 << 0 << " " << M/(nx*ny) << endl;	
    
    
    for(int i=0;i<pow(10, 8);i++){
    	E = energia(J, nx, ny, spins);
    	M = magnetizzazione(nx,ny,spins);
    	de = cambia_configurazione(J,nx,ny,spins,spins1);
    	test = numero_random();
    	if (pow(M_E, -(de/kbT))>sqrt(pow(test,2))) {
    		m = magnetizzazione(nx,ny,spins1);
    		if (sw==99) {
    			fout1 << w << " " << (E + de)/(nx*ny) << endl;
    			fout2 << w << " " << m/(nx*ny) << endl;
    			sw=0;
    		}
    		if(i>(2*pow(10,7))) {
    			Etot += (E + de)/pow(10, 8);
    			Mtot += m/(pow(10, 8));
    			Etot2 += pow((E+de),2)/pow(10, 8);
    			Mtot2 += pow(m,2)/pow(10, 8);
    		}
    		for(int i=0;i< nx*ny;i++){
        		spins[i]=spins1[i];
    		}
		} else {
			if (sw==99) {
				fout1 << w << " " << E/(nx*ny) << endl;
    			fout2 << w << " " << M/(nx*ny) << endl;	
    			sw=0;
    		}
    		if(i>(2*pow(10,7))) {
				Etot += (E)/pow(10, 8);
    			Mtot += M/(pow(10, 8));	
    			Etot2 += pow(E,2)/pow(10, 8);
    			Mtot2 += pow(M,2)/pow(10, 8);
    		}
		}
		w++;
		sw++;
		
        
    }
    cout << "Media delle energie per sito: " << Etot/(nx*ny) << endl;
    cout << "Media della magnetizzazione per sito: " << Mtot/(nx*ny) << endl;
    cout << "Errore associato alle energie per sito: " << (Etot2/pow((nx*ny),2)) - pow(Etot/(nx*ny),2) << endl;
    cout << "Errore associato alla magnetizzazione per sito: " << (Mtot2/pow((nx*ny),2)) - pow(Mtot/(nx*ny),2) << endl;
    
    
   
//    delete []  en;
//    delete [] mag;
    delete [] spins;
    delete [] spins1;
   
    
    return 0;
}

