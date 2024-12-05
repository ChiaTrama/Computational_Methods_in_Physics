
#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>

using namespace std;

//notazione come su dispense: Ma=b
void solve_tridiagonal(int n, complex<double> *d, complex<double> *u, complex<double> *l,complex<double> *b, complex<double> *a ){
  
complex<double>* alfa = new complex<double>[n];
complex<double>* beta = new complex<double>[n];

//alfa e beta primo elemento
 alfa[0]=-d[0]/u[0];
 beta[0]=b[0]/u[0];

//calcolati i due alfa e beta iniziali, calcolo con la formula nelle dispense alfa e beta per i punti centrali
 for(int i=1;i < n; i++){
 alfa[i]=(-l[i-1]/(u[i]*alfa[i-1])-d[i]/u[i]);
beta[i]=b[i]/u[i]+l[i-1]*beta[i-1]/(u[i]*alfa[i-1]);              
  }

//alfa e beta ultimo elemento
 a[n-1]=(b[n-1]/l[n-2]+beta[n-2]/alfa[n-2])/(1./alfa[n-2]+d[n-1]/l[n-2]);
    
//a questo punto avendo alfa e beta posso calcolarmi a per ogni i escluso primo ed ultimo 
for(int i=n-1;i>0;i--){
a[i-1]=a[i]/alfa[i-1]-beta[i-1]/alfa[i-1];
   }
    
    delete [] alfa; 
    delete [] beta;
     
}


int main(int argc, const char * argv[])
{
    double L=500;
    double x0=200;
    double q=2;
    double sigma=20;
    double a=250;
    double b=260;
    double V0=1.7;
    double dt=0.1;
    int Nsteps=2000;
    int Nx,Nprint;
    complex<double> pic (0.,1.0);// coppia ordinata di oggetti, entrambi di tipo double, di cui il primo rappresenta la parte reale di un numero complesso e il secondo rappresenta la parte immaginaria 

    cout << "Lunghezza L:\n";
    cout << L <<endl;
    cout << "Posizione x0:\n";
    cout << x0 <<endl;
    cout << "Momento q:\n";
    cout << q<<endl;
    cout << "Sigma: \n";
    cout << sigma<<endl;
    cout << "Limite potenziale a:\n";
    cout << a<<endl;
    cout << "Limite potenziale b:\n";
    cout << b<<endl;
    cout << "Valore potenziale V0:\n";
    cout << V0<<endl;
    cout << "Time step dt: \n";
    cout << dt<<endl;
    cout << "Numero time steps:\n";
    cout << Nsteps<<endl;
    cout << "Scrivi il numero punti asse x Nx(int), attenta che deve essere maggiore o uguale a 1000:\n";
    cin >> Nx;
    cout << "Scrivi ogni numero passi(int):\n";
    cin >> Nprint;
//Nx E' IL NUMERO TOTALE DI GRID STEPS COMPRESI GLI ESTREMI

//se ho tre gridstep per esempio ovviamente avro solo due intervalli     
    double h=L/(Nx-1);

//NMAT E' IL NUMERO DI GRID STEPS INTERNI DOVE LA FUNZIONE D'ONDA PUO'
//ESSERE DIVERSA DA ZERO
    int Nmat=Nx-2;
    double norm;
    
    
    complex<double> *psi0 = new complex<double>[Nmat];
    complex<double> *psi1 = new complex<double>[Nmat];
    complex<double> *d = new complex<double>[Nmat];
    complex<double> *u = new complex<double>[Nmat];
    complex<double> *l = new complex<double>[Nmat];
    complex<double> *f = new complex<double>[Nmat];


//SETTA PSI0 e NORMALIZZA
    norm=0.;
//setto psi, quindi mi basta definire la griglia dei punti in cui non è nulla, escludo xo e xn 
    for(int i=1;i<Nx-1;i++){
        double x=h*i;
//psi0 ha dimensione Nx-2, quindi da i=0 a i=Nx-3, a i=nx-2 ho calcolato tutte le mie psi
        psi0[i-1]=exp(pic*q*x)*exp(-pow(x-x0,2)/(2*pow(sigma,2)));
        norm=norm+pow(abs(psi0[i-1]),2);    
    }
//la norma la trovo come l'integrale del modulo quadro della funzione d'onda,poi devo farci la radice  e a quel punto ho la costante di normalizzazione per cui devo dividere la funzione d'onda in modo che sia normalizzata (uso metodo rettangoli naif)
    norm=norm*L/Nx;
    norm=sqrt(norm);
    cout << "Norma" <<norm <<'\n';
// normalizzo la funzione d'onda 
//SE VOGLIO PLOTTARE IL TUTTO NON NORMALIZZATO BASTA TOGLIERE METTERE IN COMMENTO QUESTE DUE RIGHE DEL CICLO FOR 
  for(int i=0;i<Nmat;i++){
       psi0[i]=psi0[i]/norm;
   }
    
//SETTA MATRICE M  
    for(int i=0;i<Nmat;i++){
        u[i]=1.;
        l[i]=1.;
    }
   
/*per ogni i scrivo scrivo l'elemento diagonale principale, attenta che qui stai già escludendo gli estremi perche uso Nmat=Nx-2, ovviamente per scrivere d mi serve conoscere il potenziale V in x_i */
 
for(int i=0;i<Nmat;i++){
        double Vx;
        double x=h*(i+1);// attenta che infatti non parto da x=0 ma da x=h, in x=0 ho psi nulla 

       if(x>=a && x<=b){
            Vx=V0;
       }else{
          Vx=0; }

      d[i]=pic*4.*pow(h,2)/dt-2.-2*pow(h,2)*Vx;
    }

 
// QUI CODARE LA PARTE DIAGONALE DELLA MATRICE
// STARE ATTENTI AGLI ESTREMI

    ofstream fileg;
    string nomeg;
    nomeg = (string) "psi_tutta.dat";
    fileg.open(nomeg,ios::out);
    fileg.precision(10);


  //LOOP SU PASSI
    for(int n=0;n<Nsteps;n++){ //faccio ciclo su Nsteps che sono quelli del tempo
//ovviamente in f[0] Vx è nullo
      f[0]=-psi0[1]+2.*psi0[0]+(pic*4.*pow(h,2)/dt)*psi0[0];

     for(int i=1;i<Nmat-1;i++){
           double Vx;
           double x=h*(i+1);
            if(x>=a && x<=b){
                Vx=V0;
            }else{
                Vx=0; }
            f[i]=-psi0[i+1]+2.*psi0[i]-psi0[i-1]+(pic*4.*pow(h,2)/dt)*psi0[i]+2.*pow(h,2)*Vx*psi0[i];
            
        }// qui chiudo il for sulle i ma non temporale 

  //calcolo f estremo finale, anche qui Vx è zero, sono ben lontana da b
     f[Nmat-1]=-psi0[Nmat-2]+2.*psi0[Nmat-1]+(pic*4.*pow(h,2)/dt)*psi0[Nmat-1];

       solve_tridiagonal(Nmat, d, u, l,f, psi1 );
// quindi ora ho il vettore psi riempito della funzione d'onda al tempo t=n+1

 //n che è un intero tra 0 e Nstep che sono i timestep, quindi se n, lo prendo e lo divido per N print che sarà tipo 10  ho il risultato della divisione che sarà un tot piu un resto, se il resto è uguale a zero fai le estruzioni se no break, quindi se ho che ho un resto nullo, cioe una divisione esatta allora posso stampare
        if(n%Nprint==0){

            for(int i=0;i<Nmat;i++){
                double x=h*(i+1);
                fileg << x << "  " << pow(abs(psi0[i]),2) << '\n';  
//allora stampo il modulo quadro in funzione di x     
            }

            fileg << '\n';

        }
        
 // a questo punto rimetto la norma a zero 
        norm=0.;
//e poi per i che va da 0 a Nmat-1 metto psi0 uguale a psi1, per poter andare avanti con il tempo
        for(int i=0;i<Nmat;i++){

            psi0[i]=psi1[i];
            norm+=pow(abs(psi0[i]),2);
        }
        norm=norm*L/Nx;
        norm=sqrt(norm);
// ricalcolo la norma a questo nuovo istante di tempo, la norma è data per la somma su tutti gli x a quell'istante, quindi scorro su i e ho t fermo 
        cout << "Passo " << n <<"  "<< "Norma " << norm << '\n';

    }// qui sto chiudendo il for con l'n quindi in tutto cio n è definito, al di fuori no obv 
    

    fileg.close();
    delete [] psi0;
    delete [] psi1;
    delete [] d;
    delete [] u;
    delete [] l;
    delete [] f;
 
    
    return 0;
}


