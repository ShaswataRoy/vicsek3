// C++ implementation to read
// file word by word
#include <bits/stdc++.h>
#include <iostream>
#include <cstdlib>
using namespace std;


const long long int N=16382;
int steps=2000, near_count;
double x[N], y[N], vx[N], vy[N], theta[N], mean_ang[N], list_ang[N];
double L=128, vc=0.5, dt=1, eta=0.15, r=1.;
// driver code
int main()
{
    // filestream variable file
    fstream file;
    string word, t, q, filename;

    // filename of the file
    filename = "data4/1999.txt";

    // opening file
    file.open(filename.c_str());

    // extracting words form the file
    int i=0;
    while (true)
    {
        if(i%4==0)
            file>>x[i/4];
        if(i%4==1)
            file>>y[i/4];
        if(i%4==2)
            file>>vx[i/4];
        if(i%4==3)
            file>>vy[i/4];

        i++;

        if(i==16382*4)
            break;
        // displaying content
    }

    cout<<x[16382]<<"\n"<<y[16382]<<"\n"<<vx[16382]<<"\n"<<vy[16382];

    return 0;
}
