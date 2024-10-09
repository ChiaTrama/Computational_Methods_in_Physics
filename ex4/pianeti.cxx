#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <cstdlib>
#include <array>
#include <cctype>

using namespace std;

int main()
{
    double N;

//inserisco numeri intervallo pianeti 
    cout << "inserire N numero di step da far fare ai pianeti" << endl;
    cin >> N;
// per pianeta 5 ne servono 1000000
// inizializzo i due array che mi serviranno per il velocity verlet
// nel primo inserisco i dati del pianeti dal file sistema.dat 

        string file;
	file = "sistema.dat";
 	ifstream input(file);
 	if (!input) cout <<"errore nell'apertura del file"<<endl;

     double p1[10][7];
     for (int i=0;i<10;i++){
        for (int j=0;j<7;j++){
	  input >> p1[i][j];
	   }
      }
// il secondo:
     double  p2[10][7]; 
    
      for (int i=0;i<10;i++){
         for (int j=0;j<7;j++){
	 p2[i][j]=0.;
	 }
  }
// ora trovo la massa totale del sistema solare
double Msistema;
Msistema=0.0;
      for(int k=0;k<10;k++){
		Msistema+=p1[k][0];
		}

//trovo velocità del cm 
double vcm[3];
	for (int i=0;i<3;i++) {
		for(int k=0;k<10;k++){
			vcm[i]+=p1[k][4+i]*p1[k][0]/Msistema;} }

// la tolgo alle varie velocità iniziali
for (int i=0;i<10;i++) {
		for(int k=0;k<3;k++){
			p1[i][4+k]=p1[i][4+k]-vcm[k];}
}
// in questo programma calcolo le interazioni tra i vari pianeti e il Sole  (è una buona prima approssimazione)
   ofstream a;
    cout << "i dati del pianeta 2 verranno salvati nel file object2.dat "<< endl;
    a.open("object2.dat");
    a.precision(5);
   ofstream b;
    cout << "i dati del pianeta 3 verranno salvati nel file object3.dat "<< endl;
    b.open("object3.dat");
    b.precision(5);
   ofstream c;
    cout << "i dati del pianeta 4 verranno salvati nel file object4.dat "<< endl;
    c.open("object4.dat");
    c.precision(5);
    ofstream d;
    cout << "i dati del pianeta 5 verranno salvati nel file object5.dat "<< endl;
    d.open("object5.dat");
    d.precision(5);
 
   

long double G = 6.67259E-20;
double dt;
dt=10000; 
// printo nei file i valori INIZIALI di posizione:
a<< p1[1][1]<<" "<<p1[1][2]<<" "<<p1[1][3]<<endl;
b<< p1[2][1]<<" "<<p1[2][2]<<" "<<p1[2][3]<<endl;
c<< p1[3][1]<<" "<<p1[3][2]<<" "<<p1[3][3]<<endl;
d<< p1[4][1]<<" "<<p1[4][2]<<" "<<p1[4][3]<<endl;


//pianeta 2
for (int q=0; q<N; q=q+100)
{
 
//MI TROVO TUTTE LE COSE CHE MI SERVONO PER IL VELOCITY VERLET:
// definisco raggio tra sole e pianeta 2
    double r1,r2;
    r1=0.0;
    r2=0.0;
   r1=sqrt(pow(p1[0][1]-p1[1][1], 2) + pow(p1[0][2]-p1[1][2], 2)+ pow (p1[0][3]-p1[1][3],2));
  
// calcolo la forza gravitazionale tra Sole e pianeta 2 
double f1[3];
f1[0]=0.0;
f1[1]=0.0;
f1[2]=0.0;
for (int j=0;j<3;j++){
f1[j] +=  -G * p1[0][0] * p1[1][0]/pow(r1,3)*(p1[1][1+j]-p1[0][1+j]);
}

//faccio il velocity verlet
p2[0][0]=p1[0][0];// massa sole rimane invariata
p2[1][0]=p1[1][0];
// calcolo le tre nuove coordinate spaziali 

for (int j=0;j<3;j++){
// per il sole 
 p2[0][j+1]=p1[0][j+1]+p1[0][j+4]*dt+(1.0/(2.0*p1[0][0]))*pow(dt,2)*f1[j];
// per il pianeta 2
p2[1][j+1]=p1[1][j+1]+p1[1][j+4]*dt+(1.0/(2.0*p1[1][0]))*pow(dt,2)*f1[j];
}
// a questo punto ricalcolo il raggio 
r2=sqrt(pow(p2[0][1]-p2[1][1], 2) + pow(p2[0][2]-p2[1][2], 2)+ pow (p2[0][3]-p2[1][3],2));
// a questo punto rifaccio la forza gravitazionale con le nuove coordinate spaziali 
double f2[3];
f2[0]=0.0;
f2[1]=0.0;
f2[2]=0.0;
for (int j=0;j<3;j++){
f2[j] +=  -G * p2[0][0] * p2[1][0]/pow(r2,3)*(p2[1][1+j]-p2[0][1+j]);}
// allora ho tutto per fare il verlet per la velocità per le tre componenti e trovo le velocità del pianeta 2 all'istante successivo  
for (int j=0;j<3;j++){
p2[1][j+4]=p1[1][j+4]+dt*(1.0/(2.0*p1[1][0]))*(f1[j]+f2[j]);
p2[0][j+4]=p1[0][j+4]+dt*(1.0/(2.0*p1[0][0]))*(f1[j]+f2[j]);} 
// printo 
a <<" "<< p2[1][1]<<" "<<p2[1][2]<<" "<<p2[1][3]<<endl;
// ripulisco i vettori
  for (int j=0;j<7;j++){
                       p1[1][j]=p2[1][j];
			p1[0][j]=p2[0][j];}
}




// ora rifaccio lo stesso per tutti gli altri 
//pianeta 3
for (int q=0; q<N; q=q+100)
{
    double r1,r2;
    r1=0.0;
    r2=0.0;
   r1=sqrt(pow(p1[0][1]-p1[2][1], 2) + pow(p1[0][2]-p1[2][2], 2)+ pow (p1[0][3]-p1[2][3],2));
double f1[3];
f1[0]=0.0;
f1[1]=0.0;
f1[2]=0.0;
for (int j=0;j<3;j++){
f1[j] +=  -G * p1[0][0] * p1[2][0]/pow(r1,3)*(p1[2][1+j]-p1[0][1+j]);
}
p2[0][0]=p1[0][0];
p2[2][0]=p1[2][0]; 
for (int j=0;j<3;j++){
 p2[0][j+1]=p1[0][j+1]+p1[0][j+4]*dt+(1.0/(2.0*p1[0][0]))*pow(dt,2)*f1[j];
p2[2][j+1]=p1[2][j+1]+p1[2][j+4]*dt+(1.0/(2.0*p1[2][0]))*pow(dt,2)*f1[j];
}
r2=sqrt(pow(p2[0][1]-p2[2][1], 2) + pow(p2[0][2]-p2[2][2], 2)+ pow (p2[0][3]-p2[2][3],2));
double f2[3];
f2[0]=0.0;
f2[1]=0.0;
f2[2]=0.0;
for (int j=0;j<3;j++){
f2[j] +=  -G * p2[0][0] * p2[2][0]/pow(r2,3)*(p2[2][1+j]-p2[0][1+j]);}
for (int j=0;j<3;j++){
p2[2][j+4]=p1[2][j+4]+dt*(1.0/(2.0*p1[2][0]))*(f1[j]+f2[j]);
p2[0][j+4]=p1[0][j+4]+dt*(1.0/(2.0*p1[0][0]))*(f1[j]+f2[j]);} 
b <<" "<< p2[2][1]<<" "<<p2[2][2]<<" "<<p2[2][3]<<endl;
  for (int j=0;j<7;j++){
                       p1[2][j]=p2[2][j];
			p1[0][j]=p2[0][j];}
}


//pianeta 4
for (int q=0; q<N; q=q+100)
{
    double r1,r2;
    r1=0.0;
    r2=0.0;
   r1=sqrt(pow(p1[0][1]-p1[3][1], 2) + pow(p1[0][2]-p1[3][2], 2)+ pow (p1[0][3]-p1[3][3],2));
double f1[3];
f1[0]=0.0;
f1[1]=0.0;
f1[2]=0.0;
for (int j=0;j<3;j++){
f1[j] +=  -G * p1[0][0] * p1[3][0]/pow(r1,3)*(p1[3][1+j]-p1[0][1+j]);
}
p2[0][0]=p1[0][0];
p2[3][0]=p1[3][0];
for (int j=0;j<3;j++){
 p2[0][j+1]=p1[0][j+1]+p1[0][j+4]*dt+(1.0/(2.0*p1[0][0]))*pow(dt,2)*f1[j];
p2[3][j+1]=p1[3][j+1]+p1[3][j+4]*dt+(1.0/(2.0*p1[3][0]))*pow(dt,2)*f1[j];
}
r2=sqrt(pow(p2[0][1]-p2[3][1], 2) + pow(p2[0][2]-p2[3][2], 2)+ pow (p2[0][3]-p2[3][3],2));
double f2[3];
f2[0]=0.0;
f2[1]=0.0;
f2[2]=0.0;
for (int j=0;j<3;j++){
f2[j] +=  -G * p2[0][0] * p2[3][0]/pow(r2,3)*(p2[3][1+j]-p2[0][1+j]);}
for (int j=0;j<3;j++){
p2[3][j+4]=p1[3][j+4]+dt*(1.0/(2.0*p1[3][0]))*(f1[j]+f2[j]);
p2[0][j+4]=p1[0][j+4]+dt*(1.0/(2.0*p1[0][0]))*(f1[j]+f2[j]);} 
c <<" "<< p2[3][1]<<" "<<p2[3][2]<<" "<<p2[3][3]<<endl;
  for (int j=0;j<7;j++){
                       p1[3][j]=p2[3][j];
			p1[0][j]=p2[0][j];}
}


//pianeta 5
for (int q=0; q<N; q=q+100)
{
    double r1,r2;
    r1=0.0;
    r2=0.0;
   r1=sqrt(pow(p1[0][1]-p1[4][1], 2) + pow(p1[0][2]-p1[4][2], 2)+ pow (p1[0][3]-p1[4][3],2));
double f1[3];
f1[0]=0.0;
f1[1]=0.0;
f1[2]=0.0;
for (int j=0;j<3;j++){
f1[j] +=  -G * p1[0][0] * p1[4][0]/pow(r1,3)*(p1[4][1+j]-p1[0][1+j]);
}
p2[0][0]=p1[0][0];
p2[4][0]=p1[4][0];
for (int j=0;j<3;j++){ 
 p2[0][j+1]=p1[0][j+1]+p1[0][j+4]*dt+(1.0/(2.0*p1[0][0]))*pow(dt,2)*f1[j];
p2[4][j+1]=p1[4][j+1]+p1[4][j+4]*dt+(1.0/(2.0*p1[4][0]))*pow(dt,2)*f1[j];
} 
r2=sqrt(pow(p2[0][1]-p2[4][1], 2) + pow(p2[0][2]-p2[4][2], 2)+ pow (p2[0][4]-p2[3][3],2));
double f2[3];
f2[0]=0.0;
f2[1]=0.0;
f2[2]=0.0;
for (int j=0;j<3;j++){
f2[j] +=  -G * p2[0][0] * p2[4][0]/pow(r2,3)*(p2[4][1+j]-p2[0][1+j]);}
for (int j=0;j<3;j++){
p2[4][j+4]=p1[4][j+4]+dt*(1.0/(2.0*p1[4][0]))*(f1[j]+f2[j]);
p2[0][j+4]=p1[0][j+4]+dt*(1.0/(2.0*p1[0][0]))*(f1[j]+f2[j]);} 
d <<" "<< p2[4][1]<<" "<<p2[4][2]<<" "<<p2[4][3]<<endl;
  for (int j=0;j<7;j++){
                       p1[4][j]=p2[4][j];
			p1[0][j]=p2[0][j];}
}


a.close();
b.close();
c.close();
d.close();

    return 0;}

