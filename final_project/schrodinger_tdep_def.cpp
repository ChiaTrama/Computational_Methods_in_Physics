#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>

using namespace std;

void solve_tridiagonal(int n, complex<double> *d, complex<double> *u, complex<double> *l,complex<double> *b, complex<double> *a ){
    
    complex<double>* alfa = new complex<double>[n];
    complex<double>* beta = new complex<double>[n];
    
    alfa[0]=-d[0]/u[0];
    beta[0]=b[0]/u[0];
    
    for(int i=1;i < n; i++){
        alfa[i]=(-l[i-1]/(u[i]*alfa[i-1])-d[i]/u[i]);
        beta[i]=b[i]/u[i]+l[i-1]*beta[i-1]/(u[i]*alfa[i-1]);
                 
    }
    
    
    a[n-1]=(b[n-1]/l[n-2]+beta[n-2]/alfa[n-2])/(1./alfa[n-2]+d[n-1]/l[n-2]);
    
    for(int i=n-1;i>0;i--){
        a[i-1]=a[i]/alfa[i-1]-beta[i-1]/alfa[i-1];
    }
    
    
    delete [] alfa;
    delete [] beta;
    
    
}


int main(int argc, const char * argv[])
{
    double L,x0,q,sigma;
    double a,b, V0;
    double dt;
    int Nx,Nsteps,Nprint,s;
    complex<double> pic (0.,1.0);
    
    cout << "Lunghezza L\n";
    cin >> L;
    cout << "Posizione x0\n";
    cin >> x0;
    cout << "Momento q\n";
    cin >> q;
    cout << "Sigma \n";
    cin >> sigma;
    cout << "Limite potenziale a\n";
    cin >> a;
    cout << "Limite potenziale b\n";
    cin >> b;
    cout << "Valore potenziale V0\n";
    cin >> V0;
    cout << "Numero punti asse x Nx\n";
    cin >> Nx;
    cout << "Time step dt \n";
    cin >> dt;
    cout << "Numero time steps\n";
    cin >> Nsteps;
    cout << "Scrivi ogni numero passi\n";
    cin >> Nprint;
    cout << "Che metodo si vuole utilizzare?\n";
    cout << "1)Metodo Crank Nicolson\n";
    cout << "2)Metodo Eulero esplicito\n";
    cout << "3)Metodo Stormer Verlet\n";
    cin >> s;

    
    double h=L/(Nx-1);
    int Nmat=Nx-1;
    double norm, norma, norma2, normaz, normaz2;
    
    complex<double> *psi0 = new complex<double>[Nmat];
    complex<double> *psi1 = new complex<double>[Nmat];
    complex<double> *psi1a = new complex<double>[Nmat-1];
    complex<double> *psi1b = new complex<double>[Nmat-1];
    complex<double> *psi2 = new complex<double>[Nmat];
    complex<double> *d = new complex<double>[Nmat];
    complex<double> *u = new complex<double>[Nmat];
    complex<double> *l = new complex<double>[Nmat];
    complex<double> *f = new complex<double>[Nmat];
    complex<double> *f2 = new complex<double>[Nmat-1];

    
    switch (s) {
        case 1:{
        	
//METODO DI CRANK NICOLSON
        	
        	int w=0;
    		cout << "Senza (1) o con (2) normalizzazione?\n";
    		cin >> w;
    		switch (w) {
    			case 1:{	
							
//Set di psi0 e relativa normalizzazione

    				norm=0.;
    				for(int i=0;i<Nmat;i++){
        				double x=h*i;     
        				psi0[i]=exp(pic*q*x)*exp(-pow(x-x0,2)/(2*pow(sigma,2)));
        				norm+=pow(abs(psi0[i]),2);       
    				}
    
    				norm=norm*h;
    				norm=sqrt(norm);
    				cout << "Norma" <<norm <<'\n';
    				
    				for(int i=0;i<Nmat;i++){
        				psi0[i]=psi0[i]/norm;
    				} 
    				
//Set della matrice M
    
   					for(int i=0;i<Nmat;i++){
        				u[i]=1.;
        				l[i]=1.;
    				}
    
   
      
    				for(int i=0;i<Nmat;i++){
        				double Vx;
        				double x=h*i;
        				if(x>=a && x<=b){
        		    		Vx=V0;
        				}else{
            				Vx=0;
        				}
        				d[i]=pic*4.*pow(h,2)/dt-2.-2*pow(h,2)*Vx;
    				}

    				ofstream fileg("psitot_CN_1.dat");

//LOOP SU PASSI

    				for(int n=0;n<Nsteps;n++){
        				for(int i=1;i<Nmat-1;i++){
            				double Vx;
            				double x=h*i;

            				if(x>=a && x<=b){
                				Vx=V0;
            				}else{
                				Vx=0;
            				}
            				f[i]=-psi0[i+1]+2.*psi0[i]-psi0[i-1]+(pic*4.*pow(h,2)/dt)*psi0[i]+2.*pow(h,2)*Vx*psi0[i];
            
        				}
        				
        				f[0]=-psi0[Nmat-1]-psi0[1]+2.*psi0[0]+(pic*4.*pow(h,2)/dt)*psi0[0];
        				f[Nmat-1]=-psi0[0]-psi0[Nmat-2]+2.*psi0[Nmat-1]+(pic*4.*pow(h,2)/dt)*psi0[Nmat-1];
        				f2[0]=-1.;
        				f2[Nmat-2]=-1.;
        				solve_tridiagonal(Nmat-1, d, u, l,f, psi1a );
        				solve_tridiagonal(Nmat-1, d, u, l,f2, psi1b );
        		
        				psi1[Nmat-1]=(f[Nmat-1]-psi1a[0]-psi1a[Nmat-2])/(d[Nmat-1]+psi1b[0]+psi1b[Nmat-2]);
        		
        				for(int i=0;i<Nmat-1;i++){
        					psi1[i]=psi1a[i]+psi1b[i]*psi1[Nmat-1];
        				}
        		
        
        				if(n%Nprint==0){
            				for(int i=0;i<Nmat;i++){
                				double x=h*i;
                				fileg << x << "  " << pow(abs(psi0[i]),2) << '\n';
            				}
            				fileg << '\n';
        				}
        
        
        				norma=0.;
        				for(int i=0;i<Nmat;i++){
            				psi0[i]=psi1[i];
            				norma+=pow(abs(psi0[i]),2);
        				}
        				norma=norma*L/Nmat;
        				cout << "Passo " << n << "Norma " << norma << '\n';

    				}
    
    				delete [] psi0;
    				delete [] psi1;
    				delete [] psi1a;
    				delete [] psi1b;
    				delete [] f2;
    				delete [] d;
    				delete [] u;
    				delete [] l;
    				delete [] f;
    				}break;
    			case 2:{
    				
//Set di psi0 e relativa normalizzazione

    				norm=0.;
    				for(int i=0;i<Nmat;i++){
        				double x=h*i;        
        				psi0[i]=exp(pic*q*x)*exp(-pow(x-x0,2)/(2*pow(sigma,2)));
        				norm+=pow(abs(psi0[i]),2);        
    				}
    
    				norm=norm*h;
    				norm=sqrt(norm);
    				cout << "Norma" <<norm <<'\n';
    				for(int i=0;i<Nmat;i++){
        				psi0[i]=psi0[i]/norm;
    				} 
    				
//Set della matrice M
    
   					for(int i=0;i<Nmat;i++){
        				u[i]=1.;
        				l[i]=1.;
    				}
         
    				for(int i=0;i<Nmat;i++){
        				double Vx;
        				double x=h*i;
        				if(x>=a && x<=b){
        		    		Vx=V0;
        				}else{
            				Vx=0;
        				}
        				d[i]=pic*4.*pow(h,2)/dt-2.-2*pow(h,2)*Vx;
    				}

    				ofstream fileg1("psitot_CN_2.dat");
    
//LOOP SU PASSI

    				for(int n=0;n<Nsteps;n++){
        				for(int i=1;i<Nmat-1;i++){
            				double Vx;
            				double x=h*i;

            				if(x>=a && x<=b){
                				Vx=V0;
            				}else{
                				Vx=0;
            				}
            				f[i]=-psi0[i+1]+2.*psi0[i]-psi0[i-1]+(pic*4.*pow(h,2)/dt)*psi0[i]+2.*pow(h,2)*Vx*psi0[i];
            
        				}
        				
        				f[0]=-psi0[Nmat-1]-psi0[1]+2.*psi0[0]+(pic*4.*pow(h,2)/dt)*psi0[0];
        				f[Nmat-1]=-psi0[0]-psi0[Nmat-2]+2.*psi0[Nmat-1]+(pic*4.*pow(h,2)/dt)*psi0[Nmat-1];
        				f2[0]=-1.;
        				f2[Nmat-2]=-1.;
        				solve_tridiagonal(Nmat-1, d, u, l,f, psi1a );
        				solve_tridiagonal(Nmat-1, d, u, l,f2, psi1b );
        		
        				psi1[Nmat-1]=(f[Nmat-1]-psi1a[0]-psi1a[Nmat-2])/(d[Nmat-1]+psi1b[0]+psi1b[Nmat-2]);
        		
        				for(int i=0;i<Nmat-1;i++){
        					psi1[i]=psi1a[i]+psi1b[i]*psi1[Nmat-1];
        				}
        		
        
        				if(n%Nprint==0){
            				for(int i=0;i<Nmat;i++){
                				double x=h*i;
                				fileg1 << x << "  " << pow(abs(psi0[i]),2) << '\n';              
            				}
            				fileg1 << '\n';
        				}
        
        				norma=0.;
        				for(int i=0;i<Nmat;i++){
            				psi0[i]=psi1[i];
            				norma+=pow(abs(psi0[i]),2);
        				}
        				norma=norma*L/Nmat;
        				norma=sqrt(norma);
        				for(int i=0;i<Nmat;i++){
        					psi0[i]=psi0[i]/norma;
    					}
        				cout << "Passo " << n << "Norma " << pow(norma,2) << '\n';

    				}
    
    				delete [] psi0;
    				delete [] psi1;
    				delete [] psi1a;
    				delete [] psi1b;
    				delete [] f2;
    				delete [] d;
    				delete [] u;
    				delete [] l;
    				delete [] f;
					}break;
				}
    		}break;

//METODO EULERO ESPLICITO
    
    	case 2:{
    		int p=0;
    		cout << "Senza (1) o con (2) normalizzazione?\n";
    		cin >> p;
    		switch (p) {
				case 1:{
					
//Set di psi0 e relativa normalizzazione

					ofstream fileg2("psitot_ES_1.dat");
    				norm=0.;
    				for(int i=0;i<Nmat;i++){
        				double x=h*i;       
        				psi0[i]=exp(pic*q*x)*exp(-pow(x-x0,2)/(2*pow(sigma,2)));
        				norm+=pow(abs(psi0[i]),2);
        
    				}
    
    				norm=norm*h;
    				norm=sqrt(norm);
    				cout << "Norma" <<norm <<'\n';
    				for(int i=0;i<Nmat;i++){
        				psi0[i]=psi0[i]/norm;
    				}
    				
//Loop sui passi    				
					 
    				for(int n=0;n<Nsteps;n++){
    					for (int i=1;i<Nmat-1;i++) {
    						double Vx;
        					double x=h*i;
        					if(x>=a && x<=b){
        		    			Vx=V0;
        					}else{
            					Vx=0;
        					}
    						psi1[i]=psi0[i]-(pic*dt)*(-(1./2.)*((psi0[i+1]-2.*psi0[i]+psi0[i-1])/(pow(h,2)))+Vx*psi0[i]);
						}
						
						double Vx1;
        				double x1=0.;
        				if(x1>=a && x1<=b){
        					Vx1=V0;
        				}else{
            				Vx1=0;
        				}
						psi1[0]=psi0[0]-(pic*dt)*(-(1./2.)*((psi0[1]-2.*psi0[0]+psi0[Nmat-1])/(pow(h,2)))+Vx1*psi0[0]);
						
						double Vy;
						double y=h*(Nmat-1);
						if(y>=a && y<=b){
        					Vy=V0;
        				}else{
            				Vy=0;
        				}
        				psi1[Nmat-1]=psi0[Nmat-1]-(pic*dt)*(-(1./2.)*((psi0[0]-2.*psi0[Nmat-1]+psi0[Nmat-2])/(pow(h,2)))+Vy*psi0[Nmat-1]);
        		
    					       		
        				if(n%Nprint==0){
            				for(int i=0;i<Nmat;i++){
                				double x=h*i;
                				fileg2 << x << "  " << pow(abs(psi0[i]),2) << '\n';               
            				}
            				fileg2 << '\n';
        				}
        	
        
        				norma=0.;
        				for(int i=0;i<Nmat;i++){
            				psi0[i]=psi1[i];
            				norma+=pow(abs(psi0[i]),2);
        				}
        				norma=norma*L/Nmat;
        				cout << "Passo " << n << "Norma " << norma << '\n';

    				}
    
    				delete [] psi0;
    				delete [] psi1;
    				}break;
    			case 2:{
    				
//Set di psi0 e relativa normalizzazione    				
    				
    				ofstream fileg3("psitot_ES_2.dat");
    				norm=0.;
    				for(int i=0;i<Nmat;i++){
        				double x=h*i;       
        				psi0[i]=exp(pic*q*x)*exp(-pow(x-x0,2)/(2*pow(sigma,2)));
        				norm+=pow(abs(psi0[i]),2);       
    				}
    
    				norm=norm*h;
    				norm=sqrt(norm);
    				cout << "Norma " <<norm <<'\n';
    				for(int i=0;i<Nmat;i++){
        				psi0[i]=psi0[i]/norm;
    				}
    				
    				for(int n=0;n<Nsteps;n++){
    					for (int i=1;i<Nmat-1;i++) {
    						double Vx;
        					double x=h*i;
        					if(x>=a && x<=b){
        		    			Vx=V0;
        					}else{
            					Vx=0;
        					}
    						psi1[i]=psi0[i]-(pic*dt)*(-(1./2.)*((psi0[i+1]-2.*psi0[i]+psi0[i-1])/(pow(h,2)))+Vx*psi0[i]);
						}
						
						double Vx1;
        				double x1=0.;
        				if(x1>=a && x1<=b){
        					Vx1=V0;
        				}else{
            				Vx1=0;
        				}
						psi1[0]=psi0[0]-(pic*dt)*(-(1./2.)*((psi0[1]-2.*psi0[0]+psi0[Nmat-1])/(pow(h,2)))+Vx1*psi0[0]);
						
						double Vz;
						double z=h*(Nmat-1);
						if(z>=a && z<=b){
        					Vz=V0;
        				}else{
            				Vz=0;
        				}
        				psi1[Nmat-1]=psi0[Nmat-1]-(pic*dt)*(-(1./2.)*((psi0[0]-2.*psi0[Nmat-1]+psi0[Nmat-2])/(pow(h,2)))+Vz*psi0[Nmat-1]);
        		
      		
        				if(n%Nprint==0){
            				for(int i=0;i<Nmat;i++){
                				double x=h*i;
                				fileg3 << x << "  " << pow(abs(psi0[i]),2) << '\n';            
            				}
            				fileg3 << '\n';
        				}
        	
        				double norm3=0.;
        				for(int i=0;i<Nmat;i++){
            				psi0[i]=psi1[i];
            				norm3+=pow(abs(psi0[i]),2);
        				}
        				norm3=norm3*L/Nmat;
        				norm3=sqrt(norm3);
        				
        				for(int i=0;i<Nmat;i++){
        					psi0[i]=psi0[i]/norm3;
    					}
    					double norm4=0.;
    					for(int i=0;i<Nmat;i++){
            				norm4+=pow(abs(psi0[i]),2);
        				}
        				
						norm4=norm4*h;
        				cout << "Passo " << n << "Norma " << norm4 << '\n';

    				}
    
    				delete [] psi0;
    				delete [] psi1;
    				}break;
    			}
    			}break;
    			
//METODO VERLET 
   		
		case 3:{
			int u=0;
    		cout << "Senza (1) o con (2) normalizzazione?\n";
    		cin >> u;
    		switch (u) {
				case 1:{
					
//Set di psi0 e relativa normalizzazione

					ofstream fileg4("psitot_V_1.dat");
    				norm=0.;
    				for(int i=0;i<Nmat;i++){
        				double x=h*i;        
        				psi0[i]=exp(pic*q*x)*exp(-pow(x-x0,2)/(2*pow(sigma,2)));
        				norm+=pow(abs(psi0[i]),2);
        
    				}
    
    				norm=norm*h;
    				norm=sqrt(norm);
    				cout << "Norma" <<norm <<'\n';
    				for(int i=0;i<Nmat;i++){
        				psi0[i]=psi0[i]/norm;
    				}
    				
//Il primo punto lo calcoliamo con Eulero (motivo dovuto a espressione Metodo Stormer)

					for (int i=1;i<Nmat-1;i++) {
    					double Vx;
        				double x=h*i;
        				if(x>=a && x<=b){
        		    		Vx=V0;
        				}else{
            				Vx=0;
        				}
    					psi1[i]=psi0[i]-(pic*dt)*(-(1./2.)*((psi0[i+1]-2.*psi0[i]+psi0[i-1])/(pow(h,2)))+Vx*psi0[i]);
					}
						
					double Vx1;
        			double x1=0.;
        			if(x1>=a && x1<=b){
        				Vx1=V0;
        			}else{
            			Vx1=0;
        			}
					psi1[0]=psi0[0]-(pic*dt)*(-(1./2.)*((psi0[1]-2.*psi0[0]+psi0[Nmat-1])/(pow(h,2)))+Vx1*psi0[0]);
					
					double Vy;
					double y=h*(Nmat-1);
					if(y>=a && y<=b){
        				Vy=V0;
        			}else{
            			Vy=0;
        			}
        			psi1[Nmat-1]=psi0[Nmat-1]-(pic*dt)*(-(1./2.)*((psi0[0]-2.*psi0[Nmat-1]+psi0[Nmat-2])/(pow(h,2)))+Vy*psi0[Nmat-1]);
        			
        			norma=0.;
        			for(int i=0;i<Nmat;i++){
            			norma+=pow(abs(psi1[i]),2);
        			}
        			norma=norma*L/Nmat;
        			cout << "Passo " << 0 << " Norma " << norm << '\n';
        			
        		    				
//proseguiamo ora con algoritmo Stormer 
   					 
    				for(int n=1;n<Nsteps;n++){
    					for (int i=1;i<Nmat-1;i++) {
    						double Vx;
        					double x=h*i;
        					if(x>=a && x<=b){
        		    			Vx=V0;
        					}else{
            					Vx=0;
        					}
    						psi2[i]=psi0[i]+(2.*pic*dt)*((1./2.)*((psi1[i+1]-2.*psi1[i]+psi1[i-1])/(pow(h,2)))-Vx*psi1[i]);
						}
						
						double Vx2;
        				double x2=0.;
        				if(x2>=a && x2<=b){
        					Vx2=V0;
        				}else{
            				Vx2=0;
        				}
						psi2[0]=psi0[0]+(2.*pic*dt)*((1./2.)*((psi1[1]-2.*psi1[0]+psi1[Nmat-1])/(pow(h,2)))-Vx2*psi1[0]);
			
						double Vy2;
						double y2=h*(Nmat-1);
						if(y2>=a && y2<=b){
        					Vy2=V0;
        				}else{
            				Vy2=0;
        				}
        				psi2[Nmat-1]=psi0[Nmat-1]+(2.*pic*dt)*((1./2.)*((psi1[0]-2.*psi1[Nmat-1]+psi1[Nmat-2])/(pow(h,2)))-Vy2*psi1[Nmat-1]);
        		
      		
        				if(n%Nprint==0){
            				for(int i=0;i<Nmat;i++){
                				double x=h*i;
                				fileg4 << x << "  " << pow(abs(psi1[i]),2) << '\n';
            				}
            				fileg4 << '\n';
        				}
        	
        
        				norma2=0.;
        				for(int i=0;i<Nmat;i++){
            				psi0[i]=psi1[i];
            				psi1[i]=psi2[i];
            				norma2+=pow(abs(psi1[i]),2);
        				}
        				norma2=norma2*L/Nmat;
        				cout << "Passo " << n << "Norma " << norma2 << '\n';

    				}
    
    				delete [] psi0;
    				delete [] psi1;
    				delete [] psi2;
    				
				}break;
				case 2:{
					
//Set di psi0 e relativa normalizzazione

					ofstream fileg5("psitot_V_2.dat");
    				norm=0.;
    				for(int i=0;i<Nmat;i++){
        				double x=h*i;      
        				psi0[i]=exp(pic*q*x)*exp(-pow(x-x0,2)/(2*pow(sigma,2)));
        				norm+=pow(abs(psi0[i]),2);        
    				}
     
    				norm=norm*h;
    				norm=sqrt(norm);
    				cout << "Norma" <<norm <<'\n';
    				for(int i=0;i<Nmat;i++){
        				psi0[i]=psi0[i]/norm;
    				}
    				
//Il primo punto lo calcoliamo con Eulero (motivo dovuto a espressione Metodo Verlet)

					for (int i=1;i<Nmat-1;i++) {
    					double Vx;
        				double x=h*i;
        				if(x>=a && x<=b){
        		    		Vx=V0;
        				}else{
            				Vx=0;
        				}
    					psi1[i]=psi0[i]-(pic*dt)*(-(1./2.)*((psi0[i+1]-2.*psi0[i]+psi0[i-1])/(pow(h,2)))+Vx*psi0[i]);
					}
						
					double Vx1;
        			double x1=0.;
        			if(x1>=a && x1<=b){
        				Vx1=V0;
        			}else{
            			Vx1=0;
        			}
					psi1[0]=psi0[0]-(pic*dt)*(-(1./2.)*((psi0[1]-2.*psi0[0]+psi0[Nmat-1])/(pow(h,2)))+Vx1*psi0[0]);
			
					double Vy;
					double y=h*(Nmat-1);
					if(y>=a && y<=b){
        				Vy=V0;
        			}else{
            			Vy=0;
        			}
        			psi1[Nmat-1]=psi0[Nmat-1]-(pic*dt)*(-(1./2.)*((psi0[0]-2.*psi0[Nmat-1]+psi0[Nmat-2])/(pow(h,2)))+Vy*psi0[Nmat-1]);
        			
        			norma=0.;
        			for(int i=0;i<Nmat;i++){
            			norma+=pow(abs(psi1[i]),2);
        			}
        			
					norma=norma*L/Nmat;
        			norma=sqrt(norma);
        			for(int i=0;i<Nmat;i++){
        				psi1[i]=psi1[i]/norma;
    				} 
       			        			
        			cout << "Passo " << 0 << " Norma " << pow(norma,2) << '\n';
        			       		    				
//proseguiamo ora con algoritmo Stormer   
 					 
    				for(int n=1;n<Nsteps;n++){
    					for (int i=1;i<Nmat-1;i++) {
    						double Vx;
        					double x=h*i;
        					if(x>=a && x<=b){
        		    			Vx=V0;
        					}else{
            					Vx=0;
        					}
    						psi2[i]=psi0[i]+(2.*pic*dt)*((1./2.)*((psi1[i+1]-2.*psi1[i]+psi1[i-1])/(pow(h,2)))-Vx*psi1[i]);
						}
						
						double Vx2;
        				double x2=0.;
        				if(x2>=a && x2<=b){
        					Vx2=V0;
        				}else{
            				Vx2=0;
        				}
						psi2[0]=psi0[0]+(2.*pic*dt)*((1./2.)*((psi1[1]-2.*psi1[0]+psi1[Nmat-1])/(pow(h,2)))-Vx2*psi1[0]);
			
						double Vy2;
						double y2=h*(Nmat-1);
						if(y2>=a && y2<=b){
        					Vy2=V0;
        				}else{
            				Vy2=0;
        				}
        				psi2[Nmat-1]=psi0[Nmat-1]+(2.*pic*dt)*((1./2.)*((psi1[0]-2.*psi1[Nmat-1]+psi1[Nmat-2])/(pow(h,2)))-Vy2*psi1[Nmat-1]);
        		
       		
        				if(n%Nprint==0){
            				for(int i=0;i<Nmat;i++){
                				double x=h*i;
                				fileg5 << x << "  " << pow(abs(psi1[i]),2) << '\n';               
            				}
            				fileg5 << '\n';
        				}
        				
						normaz2=0.;
						for(int i=0;i<Nmat;i++){
            				normaz2+=pow(abs(psi2[i]),2);
        				}
        			
						normaz2=normaz2*L/Nmat;
        				normaz2=sqrt(normaz2);
        				for(int i=0;i<Nmat;i++){
        					psi2[i]=psi2[i]/normaz2;
    					} 

          				normaz=0.;
        				for(int i=0;i<Nmat;i++){
            				psi0[i]=psi1[i];
            				psi1[i]=psi2[i];
            				normaz+=pow(abs(psi1[i]),2);
        				}
        				normaz=normaz*L/Nmat;

        				cout << "Passo: " << n << " Norma: " << normaz << '\n';

    				}
    
    				delete [] psi0;
    				delete [] psi1;
    				delete [] psi2;					
					
				}break;
			}

			}break;
		}
	
    
    return 0;
}













