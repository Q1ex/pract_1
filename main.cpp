#include <iostream>
#include <vector>
#include "plant.h"

using namespace std;

void input(int *M1, int *M2, int *H1,int *H2,int *N){
    cout<<"M1="; cin>>*M1;
    cout<<"M2="; cin>>*M2;
    cout<<"H1="; cin>>*H1;
    cout<<"H2="; cin>>*H2;
    cout<<"N="; cin>>*N;
}

void find_minmax(vector<double> X,double& min,double& max){
    min = X[0];
    max = X[0];
    for(double x: X)
    {
        if (x < min) min = x;
        else if (x > max) max = x;
    }
}


void
gran(vector<double>X, int H){
    double min, max;
    find_minmax(X, min, max);
    double bin_size = (max - min)/ H;
    for(size_t i=0; i< H;i++){
        auto lo = min + i * bin_size;
        auto hi = min + (i+1) * bin_size;
        cout<<'\t'<<i+1<<":("<<lo<<","<<hi<<")"<<endl;
    }
    cout<<endl<<endl;
}


vector<vector<double>>sorting_1(vector<double>X1, vector<double>X2, int H1, int H2, int N,vector<vector<double>> X_out){
    double min1, max1,min2,max2;
    find_minmax(X1, min1, max1);
    find_minmax(X2,min2,max2);
    double bin_size1 = (max1 - min1)/ H1;
    double bin_size2 = (max2 - min2)/ H2;
    double hi1= min1+bin_size1;
    double hi2= min2+bin_size2;

    vector<vector<double>> XXL(N, vector<double> (2)); ///созд массива
    for (size_t i = 0; i < N; i++){
        XXL[i][0]=X1[i];
        XXL[i][1]=X2[i];

    }

    for(size_t i=0;i<N;i++)
    {
        size_t k=0; size_t p=0;
         double hi1= min1+bin_size1;
        double hi2= min2+bin_size2;
        while(XXL[i][0]>hi1){
            k=k+1;
            hi1=hi1+bin_size1*k;

        }
        while(XXL[i][1]>hi2){
            p=p+1;
            hi2=hi2+bin_size2*p;
        }
        X_out[p][k]++;
    }

    return X_out;

}



///fF////////////////////////////
/*vector <size_t>
sorting(vector<double>X, int H ,int N){
    double min, max;
    find_minmax(X, min, max);
    double bin_size = (max - min)/ H;


    vector<size_t> result(H);
    size_t Xn = X.size();
    for(size_t j=0; j < H; j++){
        auto lo = min + j * bin_size;
        auto hi = min + (j+1) * bin_size;
        size_t k=0;
        for(size_t i=0; i < Xn; i++){
            if((lo<=X[i]) && (X[i]<=hi)){
                result[j]++;
            }
        }
    }

    return result;

}



void output_tab(int H1,int H2,vector<size_t>sort1, vector<size_t>sort2){
    cout<<"X2\\X1";
    for(size_t k=0; k<H1; k++)
        cout<<'\t'<<k+1<<"("<<sort1[k]<<")";
    cout<<endl<<endl;
    for(size_t i=0; i<H2;i++)
    {
        cout<<i+1<<"("<<sort2[i]<<")"<<'\t';
        for(size_t j=0;j<H1;j++)
            cout<<sort1[j]+sort2[i]<<'\t';
        cout<<endl<<endl;
    }

}
*/
void output_tab(vector<vector<double>> X_out, int H1, int H2, int N){
    cout<<"X2\\X1";
    for(size_t k=0; k<H1;k++)
        cout<<'\t'<<k+1;
    cout<<endl<<endl;
    for(size_t i=0; i<H2;i++){
        cout<<i+1;
        for(size_t j=0; j<H1;j++){
            cout<<'\t'<<X_out[i][j]<<" ";
        }
        cout<<endl<<endl;
    }
    cout<<endl<<endl;
    cout<<"%";
    for(size_t k=0; k<H1;k++)
        cout<<'\t'<<k+1;
    cout<<endl<<endl;
    for(size_t i=0; i<H2;i++){
        cout<<i+1;
        for(size_t j=0; j<H1;j++){
            cout<<'\t'<<(X_out[i][j]/N)*100<<"%"<<" ";
        }
        cout<<endl<<endl;
    }

return;
}

int main()
{
    Plant plant;
    plant_init(plant);

    int M1,M2,H1,H2,N;
    double X1_min, X1_max, X2_min, X2_max;
    input(&M1,&M2,&H1,&H2,&N);
    vector<double> X1(N);
    vector<double> X2(N);
    cout<<'\t'<<"X1"<<'\t'<<'\t'<<"X2"<<endl;
    for(int i=0;i<N;i++)
    {
        int in_channel=M1;
        X1[i]=plant_measure(in_channel,plant);

        in_channel=M2;
        X2[i]=plant_measure(in_channel,plant);

        cout<<i+1<<'\t'<<X1[i]<<'\t'<<'\t'<<X2[i]<<endl;
    }
    cout<<endl<<endl;


    vector<vector<double>> XXL(N, vector<double> (2));
    for (size_t i = 0; i < N; i++){
        XXL[i][0]=X1[i];
        XXL[i][1]=X2[i];
    }

   /* for(size_t i=0;i<N;i++) ///¬џвод массива
    {
        for(size_t j=0; j<2;j++)
            cout<<XXL[i][j]<<' ';
        cout<<endl;
    }
    cout<<endl<<endl;*/
///h//////////////////////////////////////////////////////////


    find_minmax(X1, X1_min, X1_max);
    cout<<"M1:"<<'\t'<<"min="<<X1_min<<'\t'<<"max="<<X1_max<<endl;
    gran(X1,H1);
    find_minmax(X2, X2_min, X2_max);
    cout<<"M2:"<<'\t'<<"min="<<X2_min<<'\t'<<"max="<<X2_max<<endl;
    gran(X2,H2); //вывод границ нтервалов

    ///—ортировка и вывод
    vector<vector<double>> X_out(H2, vector<double> (H1));
    for(size_t i=0;i<H2;i++)
        for(size_t j=0; j<H1;j++)
            X_out[i][j]=0;

    X_out=sorting_1(X1,X2,H1,H2,N,X_out);
    output_tab(X_out,H1,H2, N);

   /* vector<size_t> sort_1=sorting(X1, H1, N);
    vector<size_t> sort_2=sorting(X2, H2, N);
    output_tab(H1,H2,sort_1,sort_2);*/






}
