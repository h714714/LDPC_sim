
% =========================================================================
% Title       : Simulator for Quasi-Cyclic LDPC codes
% File        : sim_LDPC.m
% -------------------------------------------------------------------------
% Description :
%   This file performs the main Monte-Carlo simulation procedure.
%   Encodes LDPC codes described by the codes found in the codes/ folder
%   transmits bits over an AWGN channel and calls the decoding algorithm.
% ------------------------------------------------------------------------- 
% Revisions   :
%   Date       Version  Author  Description
%   20-may-11  1.3      studer  cleanup for reproducible research
%   04-jul-07  1.2      studer  multiple bug fixes
%   02-jul-07  1.1      studer  modularized & improved version
%   05-oct-06  1.0      studer  initial version 
% -------------------------------------------------------------------------
%   (C) 2006-2011 Communication Theory Group                      
%   ETH Zurich, 8092 Zurich, Switzerland                               
%   Author: Dr. Christoph Studer (e-mail: studer@rice.edu)     
% =========================================================================

function sim_LDPC_tek(RunID,LDPC, codeword, snr, maxiter) 

  randn('state',RunID)       %different RunID leads to different rand seqs
  rand('state',RunID) 
  
  tic;
  fprintf('snr %f \n', snr);
  path = '/Users/tom/Desktop/openCL_LDPC/LDPC_OpenCL/';
  str_codeword = int2str(codeword);
  str_snr = int2str(snr*10);
  str_b = int2str(LDPC.Z);
  llr_file   = strcat(path,'LLR/LLR_',str_b,'_',str_snr,'_',str_codeword,'.txt');
  input_file = strcat(path,'input/input_',str_b,'_',str_snr,'_',str_codeword,'.txt');
  llr=fopen(llr_file,'w');
  input=fopen(input_file,'w');
  inf=fopen(strcat(path,'MATLABinf.txt'),'w');
  

  sigma2 = 10^(-snr/10);   % sigma "2" for square
  [N,M] = size(LDPC.myP);
  fprintf(inf, '%d \n',maxiter);
  fprintf(inf, '%d \n',codeword);  
  fprintf(inf,'%2.5f %d %d\n',LDPC.rate, N*LDPC.Z, M*LDPC.Z); 
  time = 0;
  for trial=1:codeword 
    if mod(trial,100) == 0
      fprintf('block, snr, codeword = %d, %d, %d \n',LDPC.Z, snr, trial);
    end
    % -- draw random bits and map to symbol
    c = gf(round(rand(1,LDPC.inf_bits)),1);
    x = c*LDPC.G;
    
    %fprintf('output input signal\n')    
    c_ind = (c~=0);
    for n = 1:LDPC.inf_bits
        fprintf(input, '%d ',c_ind(n));
    end
    fprintf(input,'\n');
    

     
    noise = randn(1,length(x)); 
    s = sign((x==0)-0.5); % mapping: 1 to -1.0 and 0 to +1.0


    
      % -- AWGN channel
      y = s + noise*sqrt(sigma2);       %amp => sqrt(10^-(1/10)) = (sqrt(10))^(-1/10)      
      % -- compute LLRs & decode    
      LLR_A2 =  2*y/sigma2;             %by min sum algo paper P.12
      % -- error rate before decoder
      ber_ref = (c~=0);
      ber_ind = (LLR_A2(1:LDPC.inf_bits)<0);
      ber_bef = sum(abs(ber_ref - ber_ind));
      %fprintf('BER before decoding = %f\n', ber_bef);
      
      
      %fprintf('output llr\n');
      for n = 1:LDPC.tot_bits
          fprintf(llr, '%3.5f ',LLR_A2(n));
      end
      fprintf(llr,'\n');
      %fprintf('finish output llr\n');
      
%{      
      tic
      switch (TxRx.Decoder.LDPC.Scheduling)
        case 'Layered', % layered schedules
          [bit_output,LLR_D2,NumC,NumV] = decLDPC_layered(TxRx,LDPC,LLR_A2);
        case 'Flooding', % flooding schedule
          [bit_output,LLR_D2,NumC,NumV] = decLDPC_flooding(TxRx,LDPC,LLR_A2);
        otherwise,
          error('Unknown TxRx.Decoder.LDPC.Scheduling method.')  
      end 
      runTime = toc;
      time = runTime+time;
      
      % -- calculate BER
      ref_output = (c==1);   
      ber_after = sum(abs(ref_output-bit_output));
%}      

    %fprintf('BER after decoding is %1.10f\n.',ber_after);

  end
  fprintf(inf, '%d \n',snr);
        
  %fprintf('output running time\n');
  fprintf(inf,'%4.10f \n', time);
  
  fclose(input);
  fclose(llr);
  fclose(inf);
  
return


