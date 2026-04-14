#include<bits/stdc++.h>
using namespace std;
const int N=3;

//Forward substitution kori-> L*y=b;
void forwardSub(const double L[N][N],const double b[N],double y[N]){
    for(int i=0;i<N;i++){
        y[i]=b[i];
        for(int j=0;j<i;j++)y[i]-=L[i][j]*y[j];
        //Dhore nei L[i][i] er value non-zero 
        y[i]/=L[i][i];
    }
}

//Backward substitution kori-> U*x=y;
void backwardSub(const double U[N][N],const double y[N],double x[N]){
    for(int i=N-1;i>=0;i--){
        x[i]=y[i];
        for(int j=i+1;j<N;j++)x[i]-=U[i][j]*x[j];
        x[i]/=U[i][i];
    }
}

//Doolittle Lu decomposition
void luDecompose(const double A[N][N],double L[N][N],double U[N][N]){
    for(int i=0;i<N;i++){
        //upper triangular u
        for(int j=i;j<N;j++){
            U[i][j]=A[i][j];
            for(int k=0;k<i;k++)U[i][j]-=L[i][k]*U[k][j];
        }
        //lower triangular L
        for(int j=i;j<N;j++){
            if(i==j)L[i][i]=1;//doolittle e diagonal 1
            else{
                L[j][i]=A[j][i];
                for(int k=0;k<i;k++)L[j][i]-=L[j][k]*U[k][i];
                L[j][i]/=U[i][i];
            }
        }
    }
}

int main(){
    double A[N][N]={
        {1,1,-1},
        {1,-2,3},
        {2,3,1}
    };

    double b1[N]={4,-6,7};
    double b2[N]={2,0,5};//one more example

    double L[N][N],U[N][N];
    luDecompose(A,L,U);

    double y[N],x[N];

    //solve for b1
    forwardSub(L,b1,y);
    backwardSub(U,y,x);

    cout<<"Solution for b1:";
    for(int i=0;i<N;i++)cout<<x[i]<<" ";
    cout<<endl;

    //solve for b2
    forwardSub(L,b2,y);
    backwardSub(U,y,x);

    cout<<"Solution for b2:";
    for(int i=0;i<N;i++)cout<<x[i]<<" ";
    cout<<endl;

    return 0;
}
