#include<iostream>
#include<map>
using namespace std;
int main()
{
    int H[21] = {1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1 , 0,0 ,0,0,0,0};
    int M = 3;
    int N = 7;
    int* rowsta = new int[M+1];
    for(int i = 0 ; i < M+1; i++)
        rowsta[i]  = 0;
    int* colsta = new int[N+1];
    for(int i = 0 ; i < N+1; i++)
        colsta[i]  = 0;
    
    map< size_t,  size_t> tree;
    size_t L = 0;
    for (size_t n=0; n<N; n++) {
        for( size_t m=0; m<M; m++) {
            size_t values = m*N + n;
            if (H[values]) {
                tree[values] = L;
                L++;
                rowsta[m+1]++;
                colsta[n+1]++;
            }
        }
    }
    int rowMax = 0;
    int colMax = 0;
    for (size_t m=0; m<M; m++) {
        rowMax = (rowsta[m+1] > rowMax) ? rowsta[m+1] : rowMax;
        rowsta[m+1] += rowsta[m];
    }
    
    for (size_t n=0; n<N; n++) {
        colMax = (colsta[n+1] > colMax) ? colsta[n+1] : colMax;
        colsta[n+1] += colsta[n];
    }
    
    
    size_t* itlver = new size_t[L];
    size_t* chkind = new size_t[L];
    size_t index = 0;
    for (map< size_t,  size_t>::iterator ptr=tree.begin(); ptr!=tree.end(); ptr++) {
        chkind[index] = ptr->first % N;
        itlver[index] = ptr->second;
        index++;
        
    }
    
    for(int i = 0; i < 21; i++){
        if (i%7 == 0) cout<<endl;
        cout<<H[i]<<" ";
       
    }
    cout<<endl;
    
    for(int i = 0; i < M+1; i++) cout<<"rowsta ["<<i<<"] ="<<rowsta[i]<<endl;
    for(int i = 0; i < N+1; i++) cout<<"colsta ["<<i<<"] ="<<colsta[i]<<endl;
    for(int i = 0; i < L; i++){
        cout<<"chkind["<<i<<"] = "<<chkind[i]<<", ";
        cout<<"itvler["<<i<<"] = "<<itlver[i]<<endl;
    }
    
    
    return 0;
}
