
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

function sim_LDPC(RunID,TxRx,LDPC, codeword, snr, maxiter) 

  randn('state',RunID)       %different RunID leads to different rand seqs
  rand('state',RunID) 

  % -- initialize
  BER = zeros(1,length(TxRx.Sim.SNR_dB_list)); 
  FER = zeros(1,length(TxRx.Sim.SNR_dB_list));
  
  tic;
  str_snr = int2str(snr*10);
  llr_file   = strcat('/Users/tomho/Desktop/openCL_LDPC/LDPC_OpenCL/LLR/LLR_origin_',str_snr,'.txt');
  infile = strcat('/Users/tomho/Desktop/openCL_LDPC/LDPC_OpenCL/input/input_origin_',str_snr, '.txt');

  llr=fopen(llr_file,'w');
  input=fopen(infile,'w');
  inf=fopen('/Users/tomho/Desktop/openCL_LDPC/LDPC_OpenCL/MATLABinf.txt','w');
  
  
  %variance = 10^(-snr/10);
  %noise = sqrt(variance)*randn(size(x));
  sigma2 = 10^(-snr/10); 
  [N,M] = size(LDPC.H_prot);
  fprintf(inf, '%d \n',maxiter);
  fprintf(inf, '%d \n',codeword);  
  fprintf(inf,'%2.5f %d %d\n',LDPC.rate, N*LDPC.Z, M*LDPC.Z); 
  time = 0;

  for trial=1:codeword
    
    fprintf('codeword %d \n', trial);
    
    % -- draw random bits and map to symbol
    %c = gf(round(rand(1,LDPC.inf_bits)),1);
    
    c = round(rand(1,LDPC.inf_bits));
    %fprintf('output input signal\n')    
    for n = 1:LDPC.inf_bits
        fprintf(input, '%d ',c(n));
    end
    fprintf(input,'\n');
    
    %x = mod(c*LDPC.G,2); % generate codeword 
    x = c*LDPC.G;
    
    noise = randn(1,length(x)); 
    s = sign((x==0)-0.5); % mapping: 1 to -1.0 and 0 to +1.0


    for k=1:length(TxRx.Sim.SNR_dB_list)
      % -- AWGN channel
      y = s + noise*sqrt(sigma2);       %amp => sqrt(10^-(1/10)) = (sqrt(10))^(-1/10)
      %ps = sum(sum((s-mean(mean(s))).^2));
      %pn = sum(sum(noise*sqrt(sigma2)).^2);
      %snr = snr + 10*log(10*(ps/pn));
      
      % -- compute LLRs & decode    
      LLR_A2 =  2*y/sigma2;             %by min sum algo paper P.12
      %LLR_A2 = s;
      
      
      % -- error rate before decoder
      ber_ref = (LLR_A2(1:LDPC.inf_bits)<0);
      ber_bef = sum(abs(ber_ref - c));

      %fprintf('BER before decoding = %f\n', ber_bef);
      %fprintf('output llr\n');
      for n = 1:LDPC.tot_bits
          fprintf(llr, '%3.5f ',LLR_A2(n));
      end
      fprintf(llr,'\n');
            
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
      tmp = sum(abs(ref_output-bit_output));
      BER(k) = BER(k) + tmp; 
      FER(k) = FER(k) + (tmp>0);
      %}
    end
    %fprintf('BER after decoding is %1.10f\n.',1234);
  end
  fprintf(inf, '%d \n',snr);
        
  fprintf('output running time\n');
  fprintf(inf,'%4.10f \n', time);
  
  fclose(input);
  fclose(llr);
  fclose(inf);
  
return


