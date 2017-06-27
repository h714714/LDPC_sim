% =========================================================================
% Title       : Excution for Quasi-Cyclic LDPC codes
% File        : main.m
% -------------------------------------------------------------------------
% Description :
%   To execute functions in the source codes of sim_LDPC provided by Studer
% -------------------------------------------------------------------------                    
%   Author: Cheng-Hsiang Ho (e-mail: b02901130@ntu.edu.tw)     
% =========================================================================

%add own path to functions
addpath('/Users/tomho/Desktop/openCL_LDPC');
addpath('/Users/tomho/Desktop/openCL_LDPC/codes');
addpath('/Users/tomho/Desktop/openCL_LDPC/param');
addpath('/Users/tomho/Desktop/openCL_LDPC/LDPC_struct');
addpath('/Users/tomho/Desktop/openCL_LDPC/mediaTek');

%setting simulation parameter
codeword = 20000;
%snr = 1;
maxiter = 50;
matlab_iter = 5;

mtk = 1;    %option for MediaTek LDPC



if mtk==1

  matName = 'block_8000_2_15.mat';
  blockSize = 40;
  LF = 12;
  table_name = 'Z=384';
  LDPC = mediaTek(matName, blockSize, LF, table_name);    %for the first execution only
  return



else
%create LDPC 
  LDPC = LDPC_11nD2_648b_R23;            %can be adjusted
  RunID = 0;                      %different random seed
  TxRx = ERR_LDPC_648b_R12_LAYERED_MPA_I5(RunID, matlab_iter); %adjust param in this function

  %start simulation;
  sim_LDPC(RunID,TxRx,LDPC, codeword, snr, maxiter);

end