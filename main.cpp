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

vector <size_t>
sorting(vector<double>X, int H, int N){
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


    find_minmax(X1, X1_min, X1_max);
    cout<<"M1:"<<'\t'<<"min="<<X1_min<<'\t'<<"max="<<X1_max<<endl;
    find_minmax(X2, X2_min, X2_max);
    cout<<"M2:"<<'\t'<<"min="<<X2_min<<'\t'<<"max="<<X2_max<<endl<<endl;

    vector<size_t> sort_1=sorting(X1, H1, N);
    vector<size_t> sort_2=sorting(X2, H2, N);
    output_tab(H1,H2,sort_1,sort_2);






}
