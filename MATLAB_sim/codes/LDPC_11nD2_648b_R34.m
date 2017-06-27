% =============================================================================
% Title       : 802.11n D2.0 LDPC Code Definition
% -----------------------------------------------------------------------------
% Revisions   :
%   Date       Version  Author  Description
%   05-jul-07  1.0      studer  file created
% =============================================================================

function LDPC = LDPC_11nD2_648b_R34

  LDPC.name = 'LDPC_11nD2_648b_R34.mat';
   
  LDPC.Z = 27; % subblock size    
            
  LDPC.H_prot = {'16' '17' '22' '24' '9' '3' '14' '-' '4' '2' '7' '-' '26' '-' '2' '-' '21' '-' '1' '0' '-' '-' '-' '-';
                 '25' '12' '12' '3' '3' '26' '6' '21' '-' '15' '22' '-' '15' '-' '4' '-' '-' '16' '-' '0' '0' '-' '-' '-';
                 '25' '18' '26' '16' '22' '23' '9' '-' '0' '-' '4' '-' '4' '-' '8' '23' '11' '-' '-' '-' '0' '0' '-' '-';
                 '9' '7' '0' '1' '17' '-' '-' '7' '3' '-' '3' '23' '-' '16' '-' '-' '21' '-' '0' '-' '-' '0' '0' '-';
                 '24' '5' '26' '7' '1' '-' '-' '15' '24' '15' '-' '8' '-' '13' '-' '13' '-' '11' '-' '-' '-' '-' '0' '0';
                 '2' '2' '19' '14' '24' '1' '15' '19' '-' '21' '-' '2' '-' '24' '-' '3' '-' '2' '1' '-' '-' '-' '-' '0'};

  LDPC = myGenmat(LDPC); % compute check and generator matrices
  
return
            