#ifndef ROW_MAX
#define ROW_MAX (32)
#endif

#ifndef COL_MAX
#define COL_MAX (32)
#endif

#ifndef EPS
#define EPS (1.175494351e-38f)
#endif

#ifndef LIM
#define LIM (18.03f)
#endif

#ifndef ALPHA
#define ALPHA (1.0f)
#endif

__kernel void vnd(
          const int N,
          const size_t L,
		__global const int* _colsta,
		__global float* _x,
		__global float* _softout,
		__global const float* _llr) {
	float x_local[COL_MAX];
	int k = get_global_id(0);
    int codeword_ind = k/N;
    int iBegin = _colsta[k % N];
	int iLength = _colsta[k % N + 1] - iBegin;
	float sum = _llr[k];
	
	for (int i = 0; i < iLength; i++) {
        x_local[i] = _x[codeword_ind*L + i+iBegin] * ALPHA;
		sum += x_local[i];
	}
	for (int i = 0; i < iLength; i++) {
		_x[codeword_ind*L+i+iBegin] = sum - x_local[i];
	}
	_softout[k] = sum;
}

__kernel void cnd(
          const int M,
          const size_t L,
		__global const int* _rowsta,
		__global const size_t* _itlv,
		__global const size_t* _chkind,
		__global float* _x) {
	size_t p_local[ROW_MAX];
	float x_local[ROW_MAX];
	
	int k = get_global_id(0);
    int codeword_ind = k/M;
    float prod = 1.0f;
	float temp;
	int iBegin = _rowsta[k%M];
	int iLength = _rowsta[k%M + 1] - iBegin;


	// horizontal process (CND)
	for (int i = 0; i < iLength; i++) {
		p_local[i] = _itlv[i+iBegin];
		x_local[i] = tanh(_x[codeword_ind*L+p_local[i]]/2 + EPS);
		prod *= x_local[i];
	}
	for (int i = 0; i < iLength; i++) {
		temp = clamp(prod/x_local[i], -1.0f, 1.0f);
		_x[codeword_ind*L+p_local[i]] = clamp(2*atanh(temp), -LIM, LIM);
	}
}


__kernel void memset(__global float* mem) {
	mem[get_global_id(0)] = 0;
}


__kernel void s2h(

                  __global const float* _sout,
                  __global int* _hout) {

    int k = get_global_id(0);
    _hout[k] = (_sout[k] < 0)? 1:0;
}

__kernel void CRC(
                  const int N,
                  const int M,
                  const int subword,
                  __global const int* _rowsta,
                  __global const size_t* _chkind,
                  __global const int* _hout,
                  __global int* _check) {
    
    int k = get_global_id(0);
    int iBegin = _rowsta[k];
    int iLength = _rowsta[k + 1] - iBegin;
    
    
    //first reset for the last codeword set
    for (int j = 0; j < subword; j++){
        _check[M*j + k] = 0;
    }
    for (int i = 0; i < iLength; i++) {
        for (int j = 0; j < subword; j++){
            _check[M*j + k] += _hout[N*j + _chkind[i+iBegin]];
        }
    }
    for (int j = 0; j < subword; j++){
        _check[M*j + k] = _check[M*j + k] % 2;
    }
}

__kernel void zeroCRC(
                  const int N,
                  const int inf_bit,
                  __global const int* _hout,
                  __global int* _zeroCheck) {
    
    int k = get_global_id(0);
    
    //first reset for the last codeword set
    _zeroCheck[k] = 0;
    
    for (int i = 0; i < inf_bit; i++) {
        if(_hout[k*N + i] != 0){
            _zeroCheck[k] = 1;
            i = N;
        }
    }
}


