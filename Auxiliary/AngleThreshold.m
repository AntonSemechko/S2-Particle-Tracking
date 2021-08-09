function [a_thr,D]=AngleThreshold(SM)
% Compute angle threshold parameter used to speed-up particle tracking on a
% unit sphere.
%
% INPUT:
%   - SM    : spherical mesh represented as an object of 'TriRep' class,
%             'triangulation' class, or a cell such that TR={Tri,V}, where
%             Tri is an M-by-3 array of faces and V is an N-by-3 array of 
%             vertex coordinates.
%
% OUTPUT:
%   - a_thr : angle threshold (in degrees). This parameter corresponds to 
%             the maximum angular separation between two vertices connected
%             by an edge.  
%   - D     : M-by-1 array of maximum geodesic distances between vertices 
%             of mesh triangles. 
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Face and vertex lists
[Tri,X,fmt]=GetMeshData(SM);

if nargout<2    
    
    % Edge lengths
    if fmt>1, SM=triangulation(Tri,X); end
    E=edges(SM);
    L=EvalGeoDist(X(E(:,1),:),X(E(:,2),:),false);
    
    % Angle threshold
    a_thr=max(L)*(180/pi);    
    
else
    
    % Maximum angular separation (in degrees) between vertices of every face
    Nf=size(Tri,1);
    D=zeros(Nf,3);
    for i=1:3
        Tri=circshift(Tri,[0 1]);
        D(:,i)=EvalGeoDist(X(Tri(:,1),:),X(Tri(:,2),:),false);
    end
    D=max(D,[],2);
    D=D*(180/pi);   
    
    % Angle threshold
    a_thr=max(D);
    
end

