% ========================================================================
% Title       : LDPC H, G generator for mediaTek code
% ========================================================================

function LDPC = revisionTableGen(blockSize, LF, M, N, table_name)

fileID = fopen(table_name,'r');

formatSpec = '%d';
tmp = fscanf(fileID,formatSpec);
tmp(tmp>LF) = mod(tmp(tmp>LF), LF);



f = strcat('/Users/tomho/Desktop/openCL_LDPC/LDPC_OpenCL/table/revisionTable');
fid=fopen(f,'w+');
table = zeros(M,N);
for i = 1:M
    for j = 1:N
        table(i,j) = tmp((i-1)*50 + j);
        fprintf(fid, '%d ',table(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);



  disp('Create check matrix H ...');
  H = [];
  for m=1:M
    row = [];
    for n=1:N
      content = table(m,n);
      if content==-1
        row = [ row zeros(blockSize) ];
      else
        row = [ row eyeperm(blockSize,content) ];
      end
      
    end
    H = [ H ; row ];

  end

  % == output H information == 
    fprintf('output H\n');
    f = strcat('/Users/tomho/Desktop/openCL_LDPC/LDPC_OpenCL/table/revisionH');
    fid=fopen(f,'w+');
    for i = 1:M*blockSize
        i
        for j = 1:N*blockSize
            fprintf(fid, '%d ',H(i,j));
        end
        fprintf(fid,'\n');
    end
    fclose(fid);


end
% -- generate permutation matrix
function PZ = eyeperm(Z,p)
  %p = mod(p,Z);
  EZ = eye(Z);
  PZ = [ EZ(:,Z-p+1:Z) EZ(:,1:Z-p) ];
return
end

