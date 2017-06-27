//
//  main.cpp
//  LDPC_OpenCL
//
//  Created by 賀正翔 on 2017/4/21.
//  Copyright © 2017年 tomho. All rights reserved.
//
#include <iostream>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <OpenCL/opencl.h>
#include <time.h>
#include <string>
#include <random>
#include "stdio.h"
using namespace std;


cl_program load_program(cl_context , const char* );


int main (int argc, const char * argv[])
{
    
    
    
    //----------------------------------------------------------
    //read input file (H, encoded codeword, llr)
    //----------------------------------------------------------
    
    size_t M;
    size_t N;
    int tableZ;
    int Z;
    string codeRate;
    int ParityPunc;
    int codeword;
    size_t subword;
    int block_M;
    int block_N;
    int block_col;
    float snr, snrMin, snrMax, snrSpace;
    float puncValue;
    int maxIter, maxBLER, maxwordIter;
    int CBS;         // the length of inf_bit without puncturing
    //clock_t start, stop;
    
    //----------------------------------------------------------
    // look up table
    //----------------------------------------------------------
    /*
        CSB        Z    table            CR      M   N
        40~47      12   384             0.33      34  50
        48~63      10   320             0.40      26  42
        64~95      16   512             0.50      18  34
        96~400     28   384             0.67      10  26
        1000       64   512             0.75       8  24
        2000       128  512             0.83       6  22
        4000       256  512             0.89       4  20
        6000       384  384
        8000       512  512
    */
    //----------------------------------------------------------
    // variables setting
    //----------------------------------------------------------
    codeword    = 100;
    subword     = 10;
    CBS         = 1088;
    Z           = 128;    // real block size
    tableZ      = 512;
    

    codeRate    = "089";
    ParityPunc  = 122;     //CR = CBS/(bolck_M*Z + CBS - 2*Z - ParityPunc)
    block_M     = 4;
    snrMin      = 7.5;
    snrMax      = 7.5;
    snrSpace    = 0.1;
    maxBLER     = 20;
    maxwordIter = 10000;
    
    block_N     = block_M + 16;
    M = block_M*Z;
    N = block_N*Z;
    puncValue   = 1e-20;
    block_col   = 60;    // size of table txt file 34*60
    maxIter = 50;

    
    
    default_random_engine generator;
    normal_distribution<double> distribution(0,1);  // standard normal (mean, divation)
    
    
    //----------------------------------------------------------
    // read the shift table
    //----------------------------------------------------------
    
    string tableName = "table/Codebook_A/MTK_eMBB_LDPC shift_coefficient_v5_A_CZ" + to_string(tableZ)+".txt";
    ifstream tableFile(tableName, ios::in);

    int* table = new int[block_M*block_N];
    int  tmp;
    for(int i = 0; i < block_M; i++){
        for(int j = 0; j < block_col; j++){
            tableFile>>tmp;
            if(j < block_N){
                if(tmp >= Z)
                    tmp = tmp % Z;
                table[i*block_N + j] = tmp;
            }
        }
    }
    tableFile.close();
    
    fstream outtable("table/read_revision_test", ios::out);
    for(int i  = 0; i < block_M; i++){
        for(int j = 0; j < block_N; j++){
            outtable<<table[i*block_N + j]<<" ";
        }
        outtable<<endl;
    }
    outtable.close();
    
    
    //----------------------------------------------------------
    //  create H data struct
    //----------------------------------------------------------
    int* rowsta = new int[M+1];
    for(int i = 0 ; i < M+1; i++)
        rowsta[i]  = 0;
    int* colsta = new int[N+1];
    for(int i = 0 ; i < N+1; i++)
        colsta[i]  = 0;
    
    map< size_t,  size_t> tree;
    size_t L = 0;
    int subcol;
    int shift;
    size_t m;
    
    for (size_t n=0; n<N; n++) {
        subcol = n % Z;
        for(int block = 0; block < block_M; block++){
            shift = table[block*block_N + n/Z];
            if(shift != -1 ){
                m = Z*block + (subcol + Z - shift);
                m = (m > Z*(block+1))? m-Z:m;
                size_t values = m*N + n;
                tree[values] = L;
                //debug<<"tree(m*N+n) = ("<<m<<"*"<<N<<"+"<<n<<"("<<values<<") = "<<L<<endl;
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
    
    
    cout<<"initial LDPC parameter success"<<endl;
    
    
    
    
    size_t tot_col = N*subword;
    size_t tot_row = M*subword;
    size_t tot_L = L*subword;

    
    for (snr = snrMin; snr <= snrMax; snr += snrSpace)
    {
        int lastError = 0;
        int secondError = 0;
        cout<<"snr = "<<snr<<endl;
        string outfile;
        outfile = "resultA/CBS" + to_string(CBS) + "/CR_"+codeRate +"/" + to_string(int(ceil(snr*10))) \
            + "_precise.txt";
        cout<<"output to "<<outfile<<endl;
        fstream output(outfile, ios::out);
        

        output<<"----------------------------------------------"<<endl;
        output<<" parameter"<<endl;
        output<<"----------------------------------------------"<<endl;
        output<<"SNR      : "<<snr<<endl;
        output<<"maxIter  : "<<maxIter<<endl;
        output<<"codeword : "<<codeword<<endl;
        output<<"subword  : "<<subword<<endl;
        output<<"codeRate : "<<codeRate<<endl;
        output<<"H size   : "<<M<<"*"<<N<<endl;
        output<<endl<<endl;

        
        cl_context context;
        cl_device_id *devices;
        char *devname;
        size_t cb;

        // create a GPU context
        context = clCreateContextFromType(NULL, CL_DEVICE_TYPE_GPU, NULL, NULL, NULL);
        if(context == 0) {
            printf("Can't create GPU context\n");
            return 0;
        }
        
        // get a list of devices
        clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, NULL, &cb);
        devices = (cl_device_id*) malloc(cb);
        clGetContextInfo(context, CL_CONTEXT_DEVICES, cb, devices, 0);
        
        // show the name of the first device
        clGetDeviceInfo(devices[0], CL_DEVICE_NAME, 0, NULL, &cb);
        devname = (char*) malloc(cb);
        clGetDeviceInfo(devices[0], CL_DEVICE_NAME, cb, devname, 0);
        printf("Device: %s\n", devname);
        
        
        cl_command_queue queue;
        cl_program program;
        
        cl_mem clm_rowsta;
        cl_mem clm_colsta;
        cl_mem clm_itlver;
        cl_mem clm_chkind;
        cl_mem clm_xram;
        cl_mem clm_sout;
        cl_mem clm_hout;
        cl_mem clm_llr;
        cl_mem clm_check;
        cl_mem clm_zeroCRC;
        
        cl_kernel vnd;
        cl_kernel cnd;
        cl_kernel CRC;
        cl_kernel memset;
        cl_kernel s2h; //soft output to hard output
        cl_kernel zeroCRC;
        cl_int error;
        
        // Create a command queue
        queue = clCreateCommandQueue(context, devices[0], 0, &error);
        // Create the compute program from the source buffer
        program = load_program(context, "Kernels_sp.cl");
        // Create the compute kernel in the program we wish to run
        vnd = clCreateKernel(program, "vnd", &error);
        cnd = clCreateKernel(program, "cnd", &error);
        s2h = clCreateKernel(program, "s2h", &error);
        CRC = clCreateKernel(program, "CRC", &error);
        zeroCRC = clCreateKernel(program, "zeroCRC", &error);
        memset = clCreateKernel(program, "memset", &error);
        
        clm_rowsta = clCreateBuffer( context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * (M+1), rowsta, NULL);
        clm_colsta = clCreateBuffer( context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * (N+1), colsta, NULL);
        clm_itlver = clCreateBuffer( context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(size_t) * L, itlver, NULL);
        clm_chkind = clCreateBuffer( context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(size_t) * L, chkind, NULL);
        clm_xram = clCreateBuffer( context, CL_MEM_READ_WRITE, sizeof(float) * L * subword, NULL, NULL);
        clm_sout = clCreateBuffer( context, CL_MEM_READ_WRITE, sizeof(float) * N * subword, NULL, NULL);
        clm_hout = clCreateBuffer( context, CL_MEM_READ_WRITE, sizeof(float) * N * subword, NULL, NULL);
        clm_llr = clCreateBuffer( context, CL_MEM_READ_ONLY, sizeof(float) * N * subword, NULL, NULL);
        clm_check = clCreateBuffer( context, CL_MEM_READ_WRITE, sizeof(int) * M * subword, NULL, NULL);
        clm_zeroCRC = clCreateBuffer( context, CL_MEM_READ_WRITE, sizeof(int) * subword, NULL, NULL);
        
  
        error  = clSetKernelArg( vnd, 0, sizeof(int),     &N);
        error |= clSetKernelArg( vnd, 1, sizeof(size_t),  &L);
        error |= clSetKernelArg( vnd, 2, sizeof(cl_mem),  &clm_colsta);
        error |= clSetKernelArg( vnd, 3, sizeof(cl_mem),  &clm_xram);
        error |= clSetKernelArg( vnd, 4, sizeof(cl_mem),  &clm_sout);
        error |= clSetKernelArg( vnd, 5, sizeof(cl_mem),  &clm_llr);
        error |= clSetKernelArg( cnd, 0, sizeof(int),     &M);
        error |= clSetKernelArg( cnd, 1, sizeof(size_t),  &L);
        error |= clSetKernelArg( cnd, 2, sizeof(cl_mem),  &clm_rowsta);
        error |= clSetKernelArg( cnd, 3, sizeof(cl_mem),  &clm_itlver);
        error |= clSetKernelArg( cnd, 4, sizeof(cl_mem),  &clm_chkind);
        error |= clSetKernelArg( cnd, 5, sizeof(cl_mem),  &clm_xram);
        error |= clSetKernelArg( s2h, 0, sizeof(cl_mem),  &clm_sout);
        error |= clSetKernelArg( s2h, 1, sizeof(cl_mem),  &clm_hout);
        error |= clSetKernelArg( CRC, 0, sizeof(int),     &N);
        error |= clSetKernelArg( CRC, 1, sizeof(int),     &M);
        error |= clSetKernelArg( CRC, 2, sizeof(int),     &subword);
        error |= clSetKernelArg( CRC, 3, sizeof(cl_mem),  &clm_rowsta);
        error |= clSetKernelArg( CRC, 4, sizeof(cl_mem),  &clm_chkind);
        error |= clSetKernelArg( CRC, 5, sizeof(cl_mem),  &clm_hout);
        error |= clSetKernelArg( CRC, 6, sizeof(cl_mem),  &clm_check);
        error |= clSetKernelArg( zeroCRC, 0, sizeof(int), &N);
        error |= clSetKernelArg( zeroCRC, 1, sizeof(int), &CBS);
        error |= clSetKernelArg( zeroCRC, 2, sizeof(cl_mem), &clm_hout);
        error |= clSetKernelArg( zeroCRC, 3, sizeof(cl_mem), &clm_zeroCRC);
        error |= clSetKernelArg( memset, 0, sizeof(cl_mem), &clm_xram);
        
        cout<<"initial opencl surrounding success"<<endl;
        
        
        
        int CRCBLER = 0;
        int worditer;
        for (worditer = 0; CRCBLER <= maxBLER && worditer <= maxwordIter; worditer ++)
        //for (worditer = 0; worditer < codeword/subword ; worditer ++)
        {
            cout<<"(worditer, snr, blockszie) = ("<<worditer<<", "<<snr<<", "<<Z<<")"<<endl<<endl;
 
            clEnqueueNDRangeKernel(queue, memset, 1, NULL, &tot_L, NULL, 0, NULL, NULL);
          
            
            //----------------------------------------------------------
            //  create zero input and add LLR
            //----------------------------------------------------------
            float* llr = new float[tot_col];
            int s = 1;                              // mapping: 1 to -1.0 and 0 to +1.0
            float sigma = pow(10,-(float(snr)/10));
            
            float noise;
            for(int  i = 0; i < tot_col; i++){
                noise = distribution(generator);
                llr[i] = 2*(s + noise*sqrt(sigma))/sigma;
            }
            
            // puncturing first two block in col
            for(int i = 0 ; i < tot_col; i+= N){
                for(int k = 0; k < 2*Z; k++){
                    llr[i+k] = puncValue;
                }
            }
            
            // zero padding from CBS to 16*Z
            for(int i = 0 ; i < tot_col; i+= N){
                // each codeword
                for(int k = CBS; k < 16*Z; k++){
                    llr[i+k] = 50;
                }
            }
            
            
            // puncturing parity check bit for code rate adjustment
            for(int i = int(N-1); i < tot_col; i+=N){
                for(int k = 0; k <ParityPunc; k++){
                    llr[i-k] = puncValue;
                }
            }
            

            clEnqueueWriteBuffer(queue, clm_llr, CL_FALSE, 0, sizeof(float) * N * subword, llr, 0, NULL, NULL);

            
            //----------------------------------------------------------
            //start iteration to decode
            //----------------------------------------------------------
            for (int iter = 0; iter < maxIter; iter++) {
                
                clEnqueueNDRangeKernel(queue, cnd, 1, NULL, &tot_row, NULL, 0, NULL, NULL);
                clEnqueueNDRangeKernel(queue, vnd, 1, NULL, &tot_col, NULL, 0, NULL, NULL);
                
            }

            //----------------------------------------------------------
            // check decoding result
            //----------------------------------------------------------
            
            clEnqueueNDRangeKernel(queue, s2h, 1, NULL, &tot_col, NULL, 0, NULL, NULL);
            int* hout = new int[N*subword];
            clEnqueueReadBuffer(queue, clm_hout, CL_TRUE, 0, sizeof(int) * N * subword, hout, 0, NULL, NULL);
            

            


            
            // zero input parallel checking
            clEnqueueNDRangeKernel(queue, zeroCRC, 1, NULL, &subword, NULL, 0, NULL, NULL);
            int* zeroCheck = new int [subword];
            clEnqueueReadBuffer(queue, clm_zeroCRC, CL_TRUE, 0, sizeof(int) * subword, zeroCheck, 0, NULL, NULL);
            
            for(int i = 0; i < subword; i++){
                if(zeroCheck[i] != 0){
                    CRCBLER++;
                    secondError = lastError;
                    lastError = int(worditer*subword + i);
                }
            }
            
            
            delete[] llr;
        }
        
        
        float preciseBLER;
        preciseBLER = (CRCBLER-1)/(0.5*(lastError + secondError));
        output<<endl<<endl;
        output<<CRCBLER<<", "<<lastError<<", "<<secondError<<endl;
        output<<"preciseBLER for snr "<<snr<<" = "<<preciseBLER<<endl;
        
        output.close();

    }
    return 0;
}


cl_program load_program(cl_context context, const char* filename)
{
    std::ifstream in(filename, std::ios_base::binary);
    if(!in.good()) {
        return 0;
    }
    
    // get file length
    in.seekg(0, std::ios_base::end);
    size_t length = in.tellg();
    in.seekg(0, std::ios_base::beg);
    
    // read program source
    std::vector<char> data(length + 1);
    in.read(&data[0], length);
    data[length] = 0;
    
    // create and build program
    const char* source = &data[0];
    cl_program program = clCreateProgramWithSource(context, 1, &source, 0, 0);
    if(program == 0) {
        return 0;
    }
    
    if(clBuildProgram(program, 0, 0, 0, 0, 0) != CL_SUCCESS) {
        return 0;
    }
    
    return program;
}





// unused functions
//----------------------------------------------------------
// read the H
//----------------------------------------------------------
/*
 ifstream Hin("table/revisionH", ios::in);
 int* H = new int[M*N];
 for(int i = 0; i < M; i++){
    for(int j = 0; j < N; j++)
        Hin>>H[i*N + j];
 }
 Hin.close();
 for (size_t n=0; n<N; n++) {
    for(int m = 0; m<M; m++){
        size_t values = m*N + n;
        if(H[values] == 1){
            tree[values] = L;
            L++;
            rowsta[m+1]++;
            colsta[n+1]++;
        }
    }
 }
 */

//----------------------------------------------------------
//  read LLR and input
//----------------------------------------------------------
// remeber to close LLR and in file if use this method

/*
float* llr = new float[tot_col];
for(int  i = 0; i < tot_col; i++){
    LLR>>llr[i];
}

for(int k = 0; k < subword; k++){
    int init_BER = 0;
    for(int i = 0; i < CBS; i++){
        if(llr[k*N+i] <= 0){
            init_BER++;
            tot_BER++;
        }
    }
    output<<"init_BER = "<<init_BER<<endl;
    }
    
    
    fstream llrtest("llrtest", ios::out);
    for(int i = 0; i < tot_col; i++){
        llrtest<<llr[i]<<" ";
    }
    llrtest.close();
    
    
    clEnqueueWriteBuffer(queue, clm_llr, CL_FALSE, 0, sizeof(float) * N * subword, llr, 0, NULL, NULL);
    
*/


// cnd vnd debug
/*
cout<<"xram"<<endl;
for(int k = 0; k < subword; k++){
    for(int i = 3000; i < 4000; i++){
        cout<<xram[i]<<" ";
    }
    cout<<endl;
}


cout<<"sout"<<endl;
//  soft to hard by sequential
float* sout = new float[N*subword];
clEnqueueReadBuffer(queue, clm_sout, CL_TRUE, 0, sizeof(float) * N * subword, sout, 0, NULL, NULL);
for(int k = 0; k < subword; k++){
    for(int i = 0; i < 5*Z; i++){
        cout<<sout[i]<<" ";
    }
    cout<<endl;
}

cout<<"===================================="<<endl;
cout<<" cnd process"<<endl;
cout<<"===================================="<<endl;
clEnqueueNDRangeKernel(queue, cnd, 1, NULL, &tot_row, NULL, 0, NULL, NULL);
cout<<endl<<endl<<"iter = "<<iter<<endl;
cout<<"xram"<<endl;
clEnqueueReadBuffer(queue, clm_xram, CL_TRUE, 0, sizeof(float) * tot_L, xram, 0, NULL, NULL);
for(int k = 0; k < subword; k++){
    for(int i = 3000; i < 4000; i++){
        cout<<xram[i]<<" ";
    }
    cout<<endl;
}


cout<<"sout"<<endl;
//  soft to hard by sequential
//float* sout = new float[N*subword];
clEnqueueReadBuffer(queue, clm_sout, CL_TRUE, 0, sizeof(float) * N * subword, sout, 0, NULL, NULL);
for(int k = 0; k < subword; k++){
    for(int i = 0; i < 5*Z; i++){
        cout<<sout[i]<<" ";
    }
    cout<<endl;
}
*/


//  check by H matrix itself
/*
int subwordBLER = 0;
start = clock();
clEnqueueNDRangeKernel(queue, CRC, 1, NULL, &M, NULL, 0, NULL, NULL);
int* check = new int [M*subword];
clEnqueueReadBuffer(queue, clm_check, CL_TRUE, 0, sizeof(int) * M * subword, check, 0, NULL, NULL);
for(int i = 0; i < subword; i++){
    int tmp = 0;
    for(int j = 0; j < M; j++){
        tmp += check[M*i + j];
    }
    if(tmp != 0){
        //output<<i<<" ";
        CRCBLER++;
        subwordBLER++;
    }
}
stop = clock();
cout<<"CRC cost time = "<<(stop-start)<<" cycle"<<endl;
output<<"BLER subword by CRC = "<<subwordBLER<<endl;

*/


// check by input and sequential
/*
for(int k = 0; k < subword; k++){
    int BER = 0;
    for(int i = 0; i < CBS; i++){
        if(hout[k*N+i] != 0){
            BER++;
        }
    }
    //output<<"BER by input checking = "<<BER<<endl;
    if(BER != 0) BLER ++;
}
*/





