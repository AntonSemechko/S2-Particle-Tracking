function [FV,VA]=FaceVertexMat(TR)
% Generate a sparse face-vertex attachment matrix. 
%
% INPUT:
%   - TR    : surface mesh represented as an object of 'TriRep' class,
%             'triangulation' class, or a cell such that TR={Tri,V}, where
%             Tri is an M-by-3 array of faces and V is an N-by-3 array of 
%             vertex coordinates.
%
% OUTPUT:
%   - FV    : Nv-by-Nf sparse logical matrix where Nv and Nf are the number
%             of vertices and faces, respectively. FV(i,j)=1 if face j is 
%             attached to vertex i, and FV(i,j)=0 otherwise.
%   - VA    : N-by-1 cell containing list of faces attached to the mesh 
%             vertices.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Face and vertex lists
[Tri,V,fmt]=GetMeshData(TR);
if fmt>1, TR=triangulation(Tri,V); end

if size(Tri,2)~=3
    error('This function is intended ONLY for triangular surface meshes')
end

Nv=size(V,1);
Nf=size(Tri,1);

VA=vertexAttachments(TR)';
i=cell(Nv,1);
for n=1:Nv
    i{n}=n*ones(numel(VA{n}),1);
end

if exist('Cell2Vec','file')==3
    i=Cell2Vec(i);
    j=Cell2Vec(VA);    
else
    i=cell2mat(i);
    j=cell2mat(VA);
end

FV=sparse(i,j,true(size(i)),Nv,Nf);
