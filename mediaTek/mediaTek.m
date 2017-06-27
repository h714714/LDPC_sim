% ========================================================================
% Title       : LDPC H, G generator for mediaTek code
% ========================================================================

function LDPC = mediaTek(matName, blockSize, LF, table_name)

% parameer
LDPC.name = matName;
LDPC.Z = blockSize;

% read sheet from table
table = xlsread('MTK_eMBB_LDPC Codebook_v4_wi_qRO.xls',table_name);
table = table([2:35], [2:51]);
table(table>LF) = mod(table(table>LF), LF);


fprintf('output table\n');
f = strcat('/Users/tomho/Desktop/openCL_LDPC/LDPC_OpenCL/table/',table_name,'.txt');
fid=fopen(f,'w+');
    [P, Q] = size(table);
      for i = 1:P
        for j = 1:Q
          fprintf(fid, '%d ',table(i,j));
        end
        fprintf(fid,'\n');
     end
  fclose(fid);

LDPC.myP = table;



  disp('Create check matrix H ...');
  [N,M] = size(table);
  LDPC.H = [];
  for n=1:N
    row = [];
    for m=1:M
      fprintf('n, m = %d, %d\n', n, m);
      content = table(n,m);
      if content==-1
        row = [ row zeros(blockSize) ];
      else
        row = [ row eyeperm(blockSize,content) ];
      end
      
    end
    LDPC.H = [ LDPC.H ; row ];

  end

  
  [LDPC.par_bits,LDPC.tot_bits] = size(LDPC.H);
  LDPC.inf_bits = LDPC.tot_bits - LDPC.par_bits;
  LDPC.rate = LDPC.inf_bits/LDPC.tot_bits;       
    
  % == output H information == 
fprintf('output H\n');
str = int2str(LDPC.Z);
f = strcat('/Users/tomho/Desktop/openCL_LDPC/LDPC_OpenCL/block_H/',str,'.txt');
fid=fopen(f,'w+');
    idx = (LDPC.H~=0);
      for i = 1:N*LDPC.Z
        %fprintf('i = %d\n', i);
        for j = 1:M*LDPC.Z
          fprintf(fid, '%d ',idx(i,j));
        end
        fprintf(fid,'\n');
     end
  fclose(fid);


%{
  % == compute parity check matrix 
  disp('Compute generator matrix G ...');
  %if gf
  LDPC.H = gf(LDPC.H,1);
  H1 = LDPC.H(:,1:LDPC.inf_bits);
  H2 = LDPC.H(:,LDPC.inf_bits+1:end);    
  LDPC.P = ((H2)\(-H1))';
   
  % -- construct parity check matrix    
  LDPC.G = [ gf(eye(LDPC.inf_bits),1) LDPC.P ];
 
  disp('check H*G^t ...');
  sum(sum(LDPC.H*(LDPC.G)' ~= 0))
  % -- store LDPC construct
  disp('Storing LDPC data to disc ...');
  save(LDPC.name,'LDPC');
  
%}  
return
end
% -- generate permutation matrix
function PZ = eyeperm(Z,p)
  %p = mod(p,Z);
  EZ = eye(Z);
  PZ = [ EZ(:,Z-p+1:Z) EZ(:,1:Z-p) ];
return
end

