#include <iostream>
#include <string>
#include <array>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>

/*PARAMETRI CONSIGLIATI
// parametri energie
    int N = 1000;
    double L = 10.;
    double E_min = 0.;
    double E_max = -2. * 10.;
    double dE = -2. * 0.001;
    double s = 1E-08;

// BUCA PIATTA  
// POTENZIALE A SCALINO
    double a = 3.;
    double b = 7.;
    double Vo = 1.;
// POTENZIALE ARMONICO
    double k = 1.;
*/
   using namespace std;

double bisezione(double Eneg, double Epos,double yneg,double ypos,double *c, double s,double h);
void autovalori(int N, double L,double c[], double E_min, double E_max,double dE, double s, vector <double> &En);
void potenziale1(double L,double N, double c[]){
    for(int i=0; i<N+1; i++){
        c[i] = 0.0;
    }
}
void potenziale2(double L,double N,double a,double b,double V0,double c[]){
	double h=L/N;
	int ia=(int)(a/L*N);
	int ib=(int)(b/L*N);
	for(int i=0; i<ia; i++){
        c[i] = 0.0;
    }
    for(int i=ia; i<ib+1; i++){
        c[i] = -2;
    }
    for(int i=ib+1; i<N+1; i++){
        c[i] = 0.0;
    }
}
void potenziale3(double L,int N,double k,double c[]){
	double x[N+1];
    for(int i=0; i<N+1; i++){
        x[i] = i*(L/N);
    }
    for(int i=0; i<N+1; i++){
        c[i] = -2.*k*(x[i]-0.5*L)*(x[i]-0.5*L);
    }
	
}
void stampanumerov (int N, double L,double c[],double Em,std::ofstream& fout);
int main(){
    

// parametri energie
    int N = 1000;
    double L = 10.;
    double E_min = 0.;
    double E_max = -2. * 10.;
    double dE = -2. * 0.001;
    double s = 1E-08;

// BUCA PIATTA
    
// POTENZIALE A SCALINO
    double a = 3.;
    double b = 7.;
    double Vo = 1.;
// POTENZIALE ARMONICO
    double k = 1.;
// SCELTA DEL TIPO DI POTENZIALE
    int caso;
	cout<<"Programma che calcola gli autovalori e stampa i file per le autofunzioni di un e- in una buca di infinita di potenziale di larghezza L"<<endl;
	cout<<endl;
    cout << "Digita 1 per buca biatta, 2 per potenziale a gradino, 3 per potenziale armonico" << endl;
    cin >> caso;
    
    if ((caso!=1) && (caso!=2) && (caso!=3)){
        while((caso!=1) && (caso!=2) && (caso!=3)){
            cout << "Errore: digita 1 per buca biatta, 2 per potenziale a gradino, 3 per potenziale armonico" << endl;
            cin >> caso;
        }
    }
//RICHIESTA PARAMETRI
	 if(caso==1){
		cout<<"Le lunghezze sono in Bohr, Le energie in Hartree"<<endl;
		cout << "Inserire un numero intero di x-step  N"<<endl;
		cin >> N;
   	    cout << "Inserire Larghezza della buca L"<<endl;
		cin >> L;
		cout << "Inserire in ordine Emin,Emax, il passo dE:"<<endl;
		cin >>E_min>>E_max>>dE;
		cout<< "Inserire la soglia per la bisezione s"<< endl;
		cin >> s;
		E_min = -2.*E_min;
   		E_max = -2. * E_max;
   		dE = -2. * dE;
    }
    
    else if(caso==2){
    	cout<<"Le lunghezze sono in Bohr, Le energie in Hartree"<<endl;
		cout << "Inserire un numero intero di x-step  N"<<endl;
		cin >> N;
   	    cout << "Inserire Larghezza della buca L"<<endl;
		cin >> L;
		cout << "Inserire in ordine Emin,Emax, il passo dE:"<<endl;
		cin >>E_min>>E_max>>dE;
		cout<< "Inserire la soglia per la bisezione s"<< endl;
		cin >> s;
		E_min = -2.*E_min;
   		E_max = -2. * E_max;
   		dE = -2. * dE;
		cout <<"Inserire in ordine l'estremo sinistro e destro del gradino"<<endl;
		cin >> a >>b;
		cout <<"Inserire altezza del gradino"<<endl;
		cin >> Vo;
    }
    
    else if(caso==3){
       	cout<<"Le lunghezze sono in Bohr, Le energie in Hartree"<<endl;
		cout << "Inserire un numero intero di x-step  N"<<endl;
		cin >> N;
   	    cout << "Inserire Larghezza della buca L"<<endl;
		cin >> L;
		cout << "Inserire in ordine Emin,Emax, il passo dE:"<<endl;
		cin >>E_min>>E_max>>dE;
		cout<< "Inserire la soglia per la bisezione s"<< endl;
		cin >> s;
		E_min = -2.*E_min;
   		E_max = -2. * E_max;
   		dE = -2. * dE;
		cout <<"Inserire il parametro k del potenziale armonico"<<endl;
		cin >> k;
    }



// COSTRUZIONE POTENZIALI
    double c[N+1];
     if(caso==1){
       potenziale1(L,N,c);
    }
    
    else if(caso==2){
        potenziale2(L,N,a,b,Vo,c);
    }
    
    else if(caso==3){
       potenziale3(L,N,k,c);
    }
// CALCOLO AUTOVALORI
   vector <double> En;
   autovalori(N,L,c,E_min,E_max,dE,s,En);
// SCRITTURA SU FILE x psi(x)
	double e;
    for (int i = 1; i < En.size()+1; i++) {
		e=En[i-1];
		std::ofstream file("autofunction" + std::to_string(i) + ".txt");
		stampanumerov (N,L,c,e,file);
		file.close();
    }
    return 0;
}

// funzioni

double bisezione(double Eneg, double Epos,double yneg,double ypos,double c[], double s,double h){
    //Energie
   	int n=1000;
   	double a=Eneg;
   	double b=Epos;
   	double ya=yneg;
   	double yb=ypos;
   	double Em=(a+b)/2.;
    long double y[n+1];
	 //Bisezione
	double errore=(a-b)*(a-b);
	double sdue=s*s;
	while(errore >= sdue){
		y[0]=0.;
		y[1]=0.1;
        for(int j=1; j<n; j++){    
            y[j+1] = (1./(12. - (Em - c[j+1])*pow(h, 2.))) * ( (10.*pow(h, 2.)*(Em - c[j]) + 24.)*y[j] + ((Em-c[j-1])*pow(h, 2.) -12.)*y[j-1]) ;
        }
        if (ya * y[n] < 0.){
            yb=y[n];
            b=Em;
            Em=(a+b)*0.5;
        }
		else if (yb*y[n]<0.)
		{
			ya=y[n];
            a=Em;
            Em=(a+b)*0.5;
		} 
		else{
			break;
		}
		errore=(a-b)*(a-b);   
    }
	return Em;		
}

void stampanumerov (int N, double L,double c[],double Em,std::ofstream& fout){
	double x[N+1];
	double y[N+1];	
	for(int i=0; i<N+1; i++){
        x[i] = i*(L/N);
    }
	double h = L/N;
	y[0]=0.;
	y[1]=0.1;
    for(int j=1; j<N; j++){    
    y[j+1] = (1./(12.-(Em - c[j+1])*pow(h, 2.)))*( (10.*pow(h,2.)*(Em - c[j]) + 24.)*y[j] + ((Em-c[j-1])*pow(h, 2.)-12.)*y[j-1]) ;
    }
	for(int i=0; i<N+1; i++){
	fout<<x[i]<<" "<<y[i]<<endl;
	}
}



void autovalori(int N, double L,double c[], double E_min, double E_max,double dE, double s, vector <double> &En){
    using namespace std;
    s=-2*1E-8;
    std::cout.precision(10);
    int n = (E_max - E_min)/dE + 1;      
    double Energie[n];
    
    for(int i=0; i<n; i++){
        Energie[i] = E_min + i*dE;
    }
    
    double x[N+1];
    
    for(int i=0; i<N+1; i++){
        x[i] = i*(L/N);
    }
    double h = L/N;
    
    double y[N+1];
    y[0] = 0.;
    y[1] = 1.;
    
    vector <double> Emaggiore;
    vector <double> Eminore;
    vector <double> ymaggiore;
    vector <double> yminore;
    
    double y_finale = 0.;
    
    for(int i=1; i<n; i++){
        
        for(int j=1; j<N; j++){   
            
            y[j+1] = (1./(12. - (Energie[i] - c[j+1])*pow(h, 2.))) * ( (10.*pow(h, 2.)*(Energie[i] - c[j]) + 24.)*y[j] + ((Energie[i]-c[j-1])*pow(h, 2.) -12.)*y[j-1]) ;
        }
        
        if (y_finale * y[N] < 0.){
            Emaggiore.push_back(Energie[i]);
            Eminore.push_back(Energie[i-1]);
            ymaggiore.push_back(y[N]);
            yminore.push_back(y_finale);
        }
        
        else{double dummy= 0;};
        
        y_finale = y[N];
        
    }
    
    int count = Emaggiore.size();      
    double Epos[count];
    double Eneg[count];
    double ypos[count];
    double yneg[count];
    
    for(int i=0; i<count; i++){
        Epos[i] = Emaggiore[i];
        Eneg[i] = Eminore[i];
        ypos[i] = ymaggiore[i];
        yneg[i] = yminore[i];
    }
	double b;
	cout<<"Autovalori"<<endl;
 	 for(int i=0; i<count; i++){
	b=-0.5*bisezione(Eneg[i],Epos[i],yneg[i],ypos[i],c,s,h);
   	cout <<"E"<<i<<"="<<b<<endl;
	En.push_back(-b*2.0);
	}   
}
