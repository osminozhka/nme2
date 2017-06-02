/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: valeria
 *
 * Created on May 31, 2017, 6:22 PM
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>

using namespace std;

double p(double x){
    return  -pow(x,2.0/3.0);
}

double q(double x){
    return pow(x,-2.0/3.0)/9.0;
}




void elimi(double**A, double*B, double*X, int m);
void thomas(double**A, double*b, double*u,int m);
void Save(double*result,double*X,int m);
const double alfa( 1.0 );
const double beta( 0.0 );
const double a( 1.0 );
const double b( 2.0 );
const int m( 1000 );


int main(int argc, char** argv) {

    //matrix F m+1 x m+1
    double **F= new double *[m+1];
    for (int i= 0;i<m+1;i++){
        F[i]= new double[m+1];
    }
    
    //vectors u,ps,x
    double *U = new double[m+1];
    double *X = new double[m+1];
    double *PS = new double[m+1];
    
  //grid
    const double h= (b-a)/m;
    for (int i=0;i<=m;i++){
        X[i]=a + i*h;        
    }
    

   
    PS[0]=alfa;
    PS[m]=beta;
    
   // fill the F matrix
    for (int i = 1;i<=m-1;i++){
        F[i][i-1] = -p(X[i])/(h*h);
        F[i][i] = (p(X[i+1])+ p(X[i]))/(h*h) + q(X[i]);
       F[i][i+1] = -p(X[i+1])/(h*h); 
    }
    
    F[0][0] = 1;
    F[m][m] = 1;
   
/*
//print F
        for (int i =0; i<=m; i++){
              for(int j=0; j<=m; j++){
                    cout <<F[i][j] <<" ";
          }
                cout << endl;

       }
  */
    //zero right hand side
    for(int i=1;i<m;i++){
        PS[i]=0;
    }
   
    
    //solve AU=PS
   
  elimi(F, PS, U, m);  

   
 //thomas(F, PS, U,m);
   Save(U,X,m);
    return 0;
}



void thomas(double**A, double*b, double*u, int n){

   
   for( int k = 1; k < n+1; k++ )
   {
      b[ k ] = ( b[ k ] - A[k][k-1] * b[ k - 1 ] ) / ( A[k][k]- A[k][k-1] * A[k-1][k] );
      if( k < n  )
      {
         A[k][k+1]= A[k][k+1]/ ( A[k][k] - A[ k ][k - 1 ] * A[k - 1] [k ] );
         A[k][k] = 1.0;
         A[k][k-1]=0.0;
      }
   }
    for(int i=0;i<n+1;i++){
                for(int j=0;j<n+1;j++){
                        cout << A[i][j] << " ";
                }
                cout << endl;
        }
   
   
   
   for( int k = n - 1; k >= 0; k-- )
      b[ k ] = b[ k ] - A[k][k+1] * b[ k + 1 ];
   
   for(int i=0;i<=n;i++){
       cout<<"U["<<i<<"]="<<b[i]<<endl;
   }
   for(int i=0;i<=n;i++){
       b[i]=u[i];
   }
   
}


void Save(double*result,double*X,int m)
{
    
    fstream file;
    file.open( "vek.txt", ios::out );
      for(int i = 0; i <= m; i++)
    {
       file  << X[i] <<"   "<< result[i] << endl;
    }
    cout << "Data was successfully saved." << endl;
}

void elimi(double**A, double*B, double*X, int m){
     double **C = new double*[m+1];
        for(int i=0; i<=m; i++){
                C[i]=new double[m+2];
        }
        //===
        //split A and B to C
        for(int i=0;i<=m;i++){
                for(int j=0;j<=m;j++){
                        C[i][j]=A[i][j];
                }
                C[i][m+1]=B[i];
        }
    
     
     
     for(int k=1;k<m;k++){
         //A
           C[k][k+1]=C[k][k+1]/(C[k][k] - C[k][k-1]*C[k-1][k]);
         //ps
           C[k][m+1]= C[k][m+1]-C[k][k-1]*C[k-1][m+1]/(C[k][k] - C[k][k-1]*C[k-1][k]);
            C[k][k-1]=0.0;
            C[k][k]=1.0;
         
     }
     
     
     /*
     for(int i=0;i<m+1;i++){
                for(int j=0;j<=m+1;j++){
                        cout << C[i][j] << " ";
                }
                cout << endl;
        }*/
     
     X[10]=beta;
      for( int k = m - 1; k >= 0; k-- )
      X[ k ] = C[ k ][m+1] - C[k][k+1] * X[ k + 1 ];
     
    // for(int i=0;i<=m;i++){
      // cout<<"U["<<i<<"]="<<X[i]<<endl;
   //}
}
