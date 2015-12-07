#include <cmath>
#include <iostream>

using namespace std;

void function(double* v, double* f);


int main () {
  
  //damit keine Verwechslungen : Vektoren a,b,c aus allgemeiner Beschreibung von RK Methoden behalten Namen
  //Konstanten a, b, c aus Aufgabe werden umbenannt zu d, e, f 
 const int N = 4; //Größe der Matrix (quadratisch)
 
  //Fülle Vektor b 
  double b[N];
  b[0] = 1.0 / 6 ;
  b[1] = 1.0 / 3 ;
  b[2] = 1.0 / 3 ;
  b[3] = 1.0 / 6 ;
  
  //Fülle Vektor c
  double c[N];
  c[0] = 0 ;
  c[1] = 0.5 ;
  c[2] = 0.5 ;
  c[3] = 1;
   

    
    
 //Zeitschritte 
    double t = 0.0; //aus Aufgabe
    double tmax = 100; //aus Aufgabe
    double dt = 0.1; //willkürlich gewählt
    int Nt = tmax/dt + 1; // Anzahl Zeitschritte 
 
 // Vektor v = (x,y,z)
    const int dim = 3;
    double v [dim];
    double v1[dim]; // für später um v = v +const zu rechnen
    double v2[dim];
    double v3[dim];
    v[0] = 1; //initial x
    v[1] = 1; // initial y
    v[2] = 1; // initial z
 
 // Vektor K1 bis K4
    double K1[dim];
    double K2[dim];
    double K3[dim];
    double K4[dim];
    
// f =(xpunkt,ypunkt,zpunkt)
    double f[dim];
    
//Ausgabe der initial values    
cout << t << "\t" << v[0] << "\t" << v[1] << "\t" << v[2] << endl; //gibt aus : t x y z

//Ausführen der Iterationen 
for (int i = 1; i < Nt; i++ ){
  t +=dt;
  
 
    function(v,f); // f hat jetzt werte von f(v)  
    
    for ( int j = 0 ; j<dim; j++){
      K1[j] = f[j];
      v1[j] = v[j] + dt*0.5*K1[j];
    }
    
    function(v1,f) ;// f hat  jetzt werte von f(v+dt/2+K1[j]) ( ich weiß magic number 0,5, aber mit Mtrix wäre es mir zu umständlich 
    
    for ( int j = 0 ; j<dim; j++){
    
      K2[j] = f[j];
      v2[j] = v[j] + dt*0.5*K2[j];
    }
       
      
    function(v2,f) ;
   
    for ( int j = 0 ; j<dim; j++){
      K3[j] = f[j];
      v3[j] = v[j] + dt*K3[j];
    
    }
       
   function(v3,f);
   
   for ( int j = 0 ; j<dim; j++){
     K4[j] = f[j];
       
   }
   
   for (int j = 0; j < dim ; j++){
      v[j] = v[j] + dt * ( b[0]*K1[j] +b[1]*K2[j] +b[2]*K3[j] +b[3]*K4[j] );
    
   } 
  cout << t << "\t" << v[0] << "\t" << v[1] << "\t" << v[2] << endl; //gibt aus : t x y z
   
  
}


    
    
    return 0;
}


void function(double* v,double* f){
    //Aus Aufgabenstellung
    //Konstanten aus Aufgabe ( a --> d, b-->e , c-->g
    const int d = 10;
    const int e = 28;
    const double g = 8.0/3;
    
    f[0] = d*(v[1] - v[0]); //xpunkt  
    f[1] = v[0]*(e-v[2]) - v[1];//ypunkt
    f[2] = v[0]*v[1]-g*v[2]; //zpunkt
  

}