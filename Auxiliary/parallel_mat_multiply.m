function AB=parallel_mat_multiply(A,B)
% Vectorized matrix-matrix and matrix-vector multiplication.
%
% INPUT
%   - A     : M-by-N-by-K array of K matrices stacked along 3rd dimension
%   - B     : 2nd input argument can be specified in one of the following 
%             formats:
%             (1) N-by-Q-by-K array of K matrices stacked along 3rd 
%                 dimension
%             (2) N-by-K array of N-dimensional K column vectors
%             (3) If B is omitted or specified as an empty array,
%                 then B will be set to A if M=N, and to 
%                 permute(A,[2 1 3]) otherwise.
%
% OUTPUT
%   - AB    : output will vary according to format of B, as described above: 
%             (1) AB is M-by-Q-by-K array so that AB(:,:,k)=A(:,:,k)*B(:,:,k)
%             (2) AB is M-by-K array so that AB(:,k)=A(:,:,k)*B(:,k)
%             (3) AB is M-by-M-by-K array so that AB(:,:,k)=A(:,:,k)^2 if 
%                 M=N and AB(:,:,k)=A(:,:,k)*A(:,:,k)' otherwise      
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%

if nargin<1 || isempty(A)
    error('Insufficient number of input arguments')
elseif ~isnumeric(A) || ndims(A)>3
    error('Invalid entry for 1st input argument (A)')
end

chk=true;
if nargin<2 || isempty(B)
    if size(A,1)==size(A,2)
        B=A;
    else
        B=permute(A,[2 1 3]);
    end
elseif ~isnumeric(B)
    error('Invalid entry for 2nd input argument (B)')
elseif ~ismatrix(B) || (ismatrix(B) && ismatrix(A))
    if size(A,2)~=size(B,1)
        error('size(B,1) must equal size(A,2) when performing matrix-matrix multiplication')
    elseif size(A,3)~=size(B,3) && size(A,3)>1
        error('size(B,3) must equal size(A,3) when performing matrix-matrix multiplication')
    end
elseif ismatrix(B)
    chk=false;
    if size(A,2)~=size(B,1)
        error('size(B,1) must equal size(A,2) when performing matrix-vector multiplication')
    elseif size(A,3)~=size(B,2)
        error('size(B,2) must equal size(A,3) when performing matrix-vector multiplication')        
    end
end

if chk % matrix-matrix multiplication
    A=permute(A,[1 2 4 3]);
    B=permute(B,[4 1 2 3]);
    AB=sum(bsxfun(@times,A,B),2);
    if size(AB,3)==1
        AB=permute(AB,[1 2 4 3]);
    else 
        AB=permute(AB,[1 3 4 2]);
    end
else % matrix-vector multiplication
    B=permute(B,[3 1 2]);
    AB=sum(bsxfun(@times,A,B),2);
    AB=permute(AB,[1 3 2]);
end


