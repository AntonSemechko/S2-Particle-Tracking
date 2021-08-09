function BinStruct=BinSphericalTriangles(SM,n,vis)
% Bin triangles of a high-resolution spherical mesh in order to speed-up 
% particle tracking on a unit sphere. Bins are determined by partitioning 
% the sphere using latitude-longitude grid. Triangle is assigned to a bin
% if any of the bin vertices or its centroid fall inside the triangle, and
% vise versa.
%
% INPUT:
%   - SM        : spherical mesh represented as an object of 'TriRep' 
%                 class, 'triangulation' class, or a cell such that 
%                 SM={Tri,V}, where Tri is an M-by-3 array of faces and V 
%                 is an N-by-3 array of vertex coordinates.
%   - n         : grid size such that n=[n_lat n_lon], where 3<=n_lat<=200
%                 and 6<=n_lon<=200 is the number of 'vertical' and 
%                 'horizontal' bins, respectively. n=[10 20] is the default
%                 setting. To account for non-uniformly distributed data
%                 n can be specified as [n_lat n_lon 1]. When the latter 
%                 setting is used, the spacing between bin edges is 
%                 automatically adjusted to match marginal density
%                 functions (along lat and lon) of point masses 
%                 corresponding to the vertices of SM.
%   - vis       : set vis=true to visualize spherical grid superimposed on
%                 top of SM. vis=false is the default setting.                 

% OUTPUT:
%   - BinStruct : structure comprised of the following fields:
%                   - mesh              <-- input spherical mesh specified as a 'triangulation' object 
%                   - normals           <-- normals of SM's faces
%                   - grid_size         <-- n=[n_lat n_lon]
%                   - BinFaceList       <-- n(1)*n(2)-by-1 cell, where each cell entry provides a list of SM faces contained in a given spherical bin
%                   - BinFaceNormalList <-- n(1)*n(2)-by-1 cell, where each cell entry contains a list of face normals corresponding to faces in BinFaceList
%                   - a_thr             <-- n(1)*n(2)-by-1 array of maximum triangle edge lengths (in degrees) contained in the bins
%                   - VA                <-- 1-by-N cell, where each cell entry contains a list of faces attached to mesh vertices
%                   - grid              <-- structure containing faces and vertices of the spherical binning grid
%                   
%                 * Note that bin indices for BinFaceList are determined by the 'AssignSphericalBins' function
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


if nargin<2 || isempty(n), n=[10 20]; end
if nargin<3 || isempty(vis), vis=false; end

n=round(n);
if ~isnumeric(n) || numel(n)<2 || numel(n)>3 || sum(isnan(n) | isinf(n))>0 || n(1)<3 || n(2)<6 || sum(n(1:2)>200)>0 
    error('Invalid entry for 2nd input argument (n)')
end

if numel(n)==2 
    n=[n 0];
elseif ~(n(3)==0 || n(3)==1)
    error('Invalid entry for 2nd input argument (n). n(3) must be 0 or 1.')
end

% Face and vertex lists
[Tri,X,fmt]=GetMeshData(SM);
if fmt>1, SM=triangulation(Tri,X); end

% Vertex attachments and face normals 
VA=vertexAttachments(SM)';
%FN=TriangleNormals(SM); % inaccurate; use triangle centroids instead 

% Construct spherical grid
grid_tol=1E-6; % minimum bin width
if n(3)==0 % uniform grid; in spherical coordinate system
    
    dt_lat=pi/n(1);
    dt_lon=(2*pi)/n(2);
    
    lat=0:dt_lat:pi;
    lon=0:dt_lon:(2*pi);
    
else % non-uniform grid adapted to data

    SM2=SM;
    Nx=size(X,1);
    while Nx<4*prod(n(1:2))
        SM2=SubdivideSphericalMesh(SM2,1);        
        Nx=size(SM2.Points,1);
    end
    X2=SM2.Points;
    
    % Cumulative mass distribution functions along latitude and longitude
    lat_data=EvalGeoDist(X2,[0 0 1]);
    lon_data=atan2(X2(:,2),X2(:,1));
    idx=lon_data<0;
    lon_data(idx)=2*pi+lon_data(idx); % unwrap from [-pi,pi] to [0,2*pi] interval
    
    % Edges of latitude bins
    lat_data=sort(lat_data,'ascend');
    W=cumsum(ones(1,Nx)/Nx);
    W_lat=n(1)*W;
    
    lat=zeros(1,n(1));
    for i=1:n(1)
        id=find(W_lat<=i,1,'last');
        lat(i)=lat_data(id);
    end
    lat(end)=pi;
    lat=[0 lat];
    
    % Make sure edges are separated by at least grid_tol radians
    de_lat=lat(2:end)-lat(1:(end-1));
    if sum(de_lat<grid_tol)>0 
        for i=1:(n(1)-1)
            if de_lat(i)<grid_tol                
                dw=grid_tol-de_lat(i);
                de_lat(i)=grid_tol;                                
                
                [y,w]=deal(de_lat((i+1):end));                
                id=y>2*grid_tol;
                w(~id)=0;                                
                if sum(w)>0 % forward error diffusion
                    de_lat((i+1):end)=y-dw*w/sum(w);
                else % backward error diffusion
                    [y,w]=deal(de_lat(1:(i-1)));
                    id=y>2*grid_tol;
                    w(~id)=0;
                    de_lat(1:(i-1))=y-dw*w/sum(w);
                end
            elseif de_lat(i)>(4/9*pi) % make sure no bin is wider than 80 degrees
                dw=de_lat(i)-(4/9*pi);
                de_lat(i)=4/9*pi;
                id=de_lat<(4/9*pi) & de_lat>2*grid_tol;
                id(i)=false;
                de_lat(id)=de_lat(id)+dw/sum(id);
            end
        end
        lat=[0 cumsum(de_lat/sum(de_lat))]*pi;        
    end    
    % ---------------------------------------------------------------------
    
            
    % Edges of longitude bins
    lon_data=sort(lon_data,'ascend');
    W_lon=n(2)*W;
    
    lon=zeros(1,n(2));
    for i=1:n(2)
        id=find(W_lon<=i,1,'last');
        lon(i)=lon_data(id);
    end
    lon(end)=2*pi;
    lon=[0 lon];
    
    % Make sure edges are separated by at least grid_tol radians
    de_lon=lon(2:end)-lon(1:(end-1));
    if sum(de_lon<grid_tol)>0 
        for i=1:(n(2)-1)
            if de_lon(i)<grid_tol                
                dw=grid_tol-de_lon(i);
                de_lon(i)=grid_tol;                                
                
                [y,w]=deal(de_lon((i+1):end));                
                id=y>2*grid_tol;
                w(~id)=0;                                
                if sum(w)>0 % forward error diffusion
                    de_lon((i+1):end)=y-dw*w/sum(w);
                else % backward error diffusion
                    [y,w]=deal(de_lon(1:(i-1)));
                    id=y>2*grid_tol;
                    w(~id)=0;
                    de_lon(1:(i-1))=y-dw*w/sum(w);
                end
            elseif de_lon(i)>(4/9*pi) % make sure no bin is wider than 80 degrees
                dw=de_lon(i)-(4/9*pi);
                de_lon(i)=4/9*pi;
                id=de_lon<(4/9*pi) & de_lon>2*grid_tol;
                id(i)=false;
                de_lon(id)=de_lon(id)+dw/sum(id);
            end
        end
        lon=[0 cumsum(de_lon/sum(de_lon))]*(2*pi);        
    end
    
    clear SM2 X2 W_lat W_lon
end

[lon_g,lat_g]=meshgrid(lon,lat);
lon_g(:,end)=0;

s11=sin(lat_g); 
x=s11.*cos(lon_g); x(end,:)=0;
y=s11.*sin(lon_g); y(end,:)=0;
z=cos(lat_g);
clear lat_g log_g

v_id=reshape(1:numel(x),n(1)+1,[]);
v_id(:,end)=v_id(:,1);

v1=v_id(1:end-1,1:end-1);
v2=v_id(2:end,1:end-1);
v3=v_id(2:end,2:end);
v4=v_id(1:end-1,2:end);

Fg=[v1(:) v2(:) v3(:) v4(:)];
Xg=[x(:) y(:) z(:)];
G=struct('faces',Fg,'vertices',Xg);

% Assign grid vertices to mesh triangles
[~,D]=AngleThreshold(SM);
D=1.001*D;
a_thr_max=max(D);
Fv_idx=SphericalTriangleIntersection(SM,Xg,a_thr_max);
Fv_idx=Fv_idx(:,1);
Fv_idx=reshape(Fv_idx(Fg(:)),[],4);

% Assign grid centroids to mesh triangles
Cg=Xg(Fg(:,1),:)+Xg(Fg(:,2),:)+Xg(Fg(:,3),:)+Xg(Fg(:,4),:);
Cg=ProjectOnSn(Cg);
Ff_idx=SphericalTriangleIntersection(SM,Cg,a_thr_max);
Ff_idx=Ff_idx(:,1);
F_idx=[Fv_idx,Ff_idx];

% Assign mesh vertices to bins 
Bv_idx=AssignSphericalBins({lat lon n(3)},X);

% Assign mesh centroids to bins
C=X(Tri(:,1),:)+X(Tri(:,2),:)+X(Tri(:,3),:);
C=ProjectOnSn(C);
FN=C;
Bf_idx=AssignSphericalBins({lat lon n(3)},C);

% Assign mesh triangles to spherical bins
[BinFaceList,BinFaceNormalList]=deal(cell(prod(n(1:2)),1));

idx_unq=unique(Bv_idx);
idx_skip=true(prod(n(1:2)),1);
idx_skip(idx_unq)=false; 

idx_F=(1:size(Tri,1))';
a_thr=zeros(prod(n(1:2)),1);

%chk_cell2vec=exist('Cell2Vec','file')==3;
for i=1:prod(n(1:2))

    idx_f=F_idx(i,:)';    
    if ~idx_skip(i)        
        idx_v=Bv_idx==i;              % mesh vertices contained in the i-th bin
        idx_f_v=cell2mat(VA(idx_v))'; % faces attached to idx_v vertices
        idx_f=cat(1,idx_f,idx_f_v);        
    end
    idx_f=cat(1,idx_f,idx_F(Bf_idx==i));
    idx_f=unique(idx_f);
    
    BinFaceList{i}=idx_f;
    BinFaceNormalList{i}=FN(idx_f,:);
    a_thr(i)=max(D(idx_f));
    
end

% Bin structure
BinStruct.mesh=SM;
BinStruct.normals=FN;
BinStruct.grid_size={lat lon n(3)};
BinStruct.BinFaceList=BinFaceList;
BinStruct.BinFaceNormalList=BinFaceNormalList;
BinStruct.a_thr=a_thr;
BinStruct.VA=VA;
BinStruct.grid=G;

% Triangle determinants 
% ------------------------

% Cross-products
X1=X(Tri(:,1),:); 
X2=X(Tri(:,2),:); 
X3=X(Tri(:,3),:);

A=X2-X1;
B=X3-X1;
C=-X1;

TD.cAB=cross(A,B,2);
TD.cCB=cross(C,B,2);
TD.cAC=cross(A,C,2);

% Determinant of [X1;X2;X3]'
A=permute(cat(3,X1,X2,X3),[2 3 1]);
detA=sum(cross(A(:,1,:),A(:,2,:),1).*A(:,3,:),1);

TD.A=A;
TD.detA=detA(:);
TD.c12=cross(X1,X2,2);

BinStruct.triangle_data=TD;



% Visualize spherical grid and mesh
% -------------------------------------------------------------------------
if ~vis(1), return; end

figure('color','w')
h=trimesh(SM); 
set(h,'EdgeColor','k','FaceColor',0.75*[1 1 1],'FaceAlpha',0.9,'EdgeAlpha',0.25,'LineWidth',0.5);
axis equal, hold on

h=patch('faces',Fg,'vertices',Xg);
set(h,'FaceColor','none','EdgeColor','b','LineWidth',2)

