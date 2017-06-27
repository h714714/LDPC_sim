function TxRx = ERR_LDPC_648b_R12_LAYERED_MPA_I5(RunID, maxiter) 

  % == LDPC SETTINGS ====================================

  TxRx.Sim.name = 'ERR_LDPC_648b_R12_LAYERED_MPA_I5';
  TxRx.Sim.nr_of_channels = 1; % 1k for good results, 
                                  % 10k for accurate results
  TxRx.Sim.SNR_dB_list = [0.2];
  TxRx.Decoder.LDPC.Scheduling = 'Flooding'; % 'Layered' and 'Flooding'
  TxRx.Decoder.LDPC.Type = 'SPA'; % 'MPA' and 'SPA' (optimal)
  TxRx.Decoder.LDPC.Iterations = maxiter;  
  %load('/Users/tomho/Desktop/LDPC_OpenCL/codes/LDPC_11nD2_648b_R12.mat'); % load code
    
  % == EXECUTE SIMULATION ===============================  
  
  %sim_LDPC(RunID,TxRx,LDPC) 
  
return
  
