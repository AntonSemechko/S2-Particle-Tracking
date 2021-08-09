function [Tuv,MDS,uvw]=SphericalTriangleIntersection_UsingStereoProj(TR,D,MDS)
% Given a set of rays, D, emanating from the origin and a triangular 
% surface mesh of a unit sphere, find mesh triangles intersected by D, 
% as well as planar and spherical barycentric coordinates of the points
% of intersection. 
%
% INPUT:
%   - TR    : spherical mesh represented as an object of 'TriRep' class,
%             'triangulation' class, or a cell such that TR={Tri,V}, where
%             Tri is an M-by-3 array of faces and V is an N-by-3 array of 
%             vertex coordinates. 
%             -------------------------- IMPORTANT------------------------- 
%             Counter-clockwise orientation of mesh triangles is assumed
%             (i.e., normals pointing towards exterior of the surface).
%             -------------------------------------------------------------
%   - D     : M-by-3 array of ray direction vectors emanating from the 
%             origin, where M is the total number of rays. Alternatively, D
%             can be thought of as points lying on the surface of a unit 
%             sphere.  
%   - MDS   : optional input argument. MDS is a data structure containing
%             various info about TR used during ray-triangle intersection 
%             tests. The purpose of MDS is to reduce the run-time 
%             'SphericalTriangleIntersection_UsingSterepProj' when calling
%             the function repeatedly with the same TR, but different D. 
%
% OUTPUT:
%   - Tuv   : N-My-3 array containing a list of triangles intersected by D
%             (1st column) and planar barycentric coordinates of the points
%             of intersection (2nd and 3rd columns). For example, 
%             ray D(i,:) intersects triangle Tuv(i,1) and the barycentric
%             coordinates of the point of intersection are u=Tuv(i,2) and
%             v=Tuv(i,3) (note: w=1-u-v). To compute the point of 
%             intersection (P) explicitly use the following formula: 
%             P=(1-u-v)*P1 + u*P2 + v*P3, where P1, P2 and P3 are the 
%             coordinates of the triangle vertices.
%   - MDS   : data structure of TR; useful when making repeated calls to
%             'SphericalTriangleIntersection_UsingSterepProj'.
%   - uvw   : M-by-3 array of spherical barycentric coordinates of the 
%             queried rays/points in D. For example, i-th ray D(i,:) 
%             intersects spherical triangle Tuv(i,1), and spherical 
%             barycentric coordinates of this point are [u,v,w]=uvw(i,:),
%             so D(i,:)=u*P1 + v*P2 + w*P3, where P1, P2 and P3 are the 
%             vertices of triangle Tuv(i,1). Note that unlike planar 
%             barycentric coordinates, u+v+w>=1.
%
% REFERENCES:
%  [1] Möller, T., Trumbore, B. (1997) 'Fast, minimum storage ray-triangle 
%      intersection', Journal of Graphics Tools, Vol. 2, pp. 21-28.
%  [2] Langer, T., Belyaev, A., Seidel, H-P. (2006) 'Spherical barycentric 
%      coordinates', In Proceedings of the 4th Eurographics Symposium on
%      Geometry Processing (SGP 2006), pp.81–88. 
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Inputs
if nargin<1 || isempty(TR)
    error('Insufficient number of input arguments')
end
[Tri,X]=GetMeshData(TR);

if nargin<3 || isempty(MDS)    
    MDS.chart_data=BreakMeshIntoCharts(Tri,X);
    MDS.triangle_data=PrecomputeTriangleData(Tri,X);    
elseif ~isstruct(MDS) || ~isfield(MDS,'chart_data') || ~isfield(MDS,'triangle_data')
    error('Invalid entry for 3rd input argument (MDS)')    
end

if nargin<2 || isempty(D)
    Tuv=[]; uvw=[];
    return
elseif ~ismatrix(D) || ~isnumeric(D) || size(D,2)~=3
    error('Invalid entry for 2nd input argument (D)')
end    

% Assign D to "chart generators"; six in total 
D=ProjectOnSn(D);
D2G=bsxfun(@minus,permute(D,[1 3 2]),permute(MDS.chart_data.generators,[3 1 2]));
D2G=sum(D2G.^2,3);
[~,G_id]=min(D2G,[],2);

% Find triangles intersected by P and (planar) barycentric coordinates of 
% the points of intersection
LUT=MDS.chart_data.LUT; 
Np=size(D,1);
Tuv=nan(Np,3);
if nargout>2, uvw=nan(Np,3); end
[tol1,tol2]=deal(1E-12,1E-12);
for g=1:6    
    
    idx=G_id==g;
    if sum(idx)==0 || isempty(MDS.chart_data.charts{g}), continue; end

    % Chart faces intersected by P(idx,:)
    Pg=D(idx,:);
    p=StereoProj(Pg,g);    
    T=pointLocation(MDS.chart_data.charts{g},p); 
    T_nan=isnan(T);
    if sum(T_nan)>0
        T(T_nan)=[];
        if isempty(T), continue; end
        Pg(T_nan,:)=[];
        idx=find(idx);
        idx(T_nan)=[];
    end    

    % Corresponding mesh faces intersected by P(idx,:)
    T=LUT{g}(T);
    T_nan=isnan(T);
    if sum(T_nan)>0
        T(T_nan)=[];
        if isempty(T), continue; end
        Pg(T_nan,:)=[];
        if islogical(idx)
            idx=find(idx);
        end
        idx(T_nan)=[];
    end
    T=MDS.chart_data.F_sub{g}(T);
    
    % Planar barycentric coordinates of the points of intersection
    uv_g=GetPlanarBCs(MDS.triangle_data,T,Pg);
    w=1-(uv_g(:,1)+uv_g(:,2));
    
    % Make sure all faces in T are indeed intersected by Pg
    idx_u=uv_g(:,1)<-tol1 | uv_g(:,1)>(1+tol1);    
    idx_v=uv_g(:,2)<-tol1 | uv_g(:,2)>(1+tol1);
    idx_w=w<-tol1;
    idx_redo=idx_u | idx_v | idx_w;
    
    Tuv(idx,:)=[T uv_g];    
    if sum(idx_redo)==0, continue; end

    % Some of the faces were assigned incorrectly so need to look at their
    % neighbours; this happens because spherical triangles are approximated
    % with planar ones during stereographic projection        
    T=T(idx_redo);
    Pg=Pg(idx_redo,:);
    if islogical(idx)
        idx=find(idx);
        idx=idx(idx_redo);
    end

    for i=1:numel(T)
        
        % Triangles in the neighbourhood of T(i)
        Ti_ngb=MDS.chart_data.VA(Tri(T(i),:));
        Ti_ngb=cell2mat(Ti_ngb);
        Ti_ngb=unique(Ti_ngb(:));
        
        % Get barycentric co-ords 
        n=numel(Ti_ngb);
        uv_i=GetPlanarBCs(MDS.triangle_data,Ti_ngb,ones(n,1)*Pg(i,:));
        idx_val=uv_i(:,1)>-tol2 & uv_i(:,1)<(1+tol2) & ...
                uv_i(:,2)>-tol2 & uv_i(:,2)<(1+tol2) & ... 
                (uv_i(:,1)+uv_i(:,2))<(1+tol2);
            
        % Select one
        idx_val=find(idx_val,1); 

        if isempty(idx_val)
          
            % No triangles in the neighbourhood contain Pg(i,:), look at
            % the entire chart
            v_id=MDS.chart_data.V_sub{g};
            Xg=X(v_id,:);

            % Triangles in the neighbourhood of Pg(i,:)
            d=sum(bsxfun(@minus,Xg,Pg(i,:)).^2,2);
            [~,id_min]=min(d);
            id_min=v_id(id_min);
            Ti_ngb=MDS.chart_data.VA{id_min}';

            % Get barycentric co-ords
            n=numel(Ti_ngb);
            uv_i=GetPlanarBCs(MDS.triangle_data,Ti_ngb,ones(n,1)*Pg(i,:));
            idx_val=uv_i(:,1)>-tol2 & uv_i(:,1)<(1+tol2) & ...
                    uv_i(:,2)>-tol2 & uv_i(:,2)<(1+tol2) & ...
                   (uv_i(:,1)+uv_i(:,2))<(1+tol2);
            idx_val=find(idx_val,1); 
            
            if isempty(idx_val)
                
                % No triangles in the chart contain Pg(i,:), look at
                % the entire mesh
                d=sum(bsxfun(@minus,X,Pg(i,:)).^2,2);
                [~,id_srt]=sort(d);
                id_srt=id_srt(1:3);
                id_srt(id_srt==id_min)=[];
                
                % Triangles in the neighbourhood of Pg(i,:)
                Ti_ngb=MDS.chart_data.VA(id_srt);
                Ti_ngb=cell2mat(Ti_ngb);
                Ti_ngb=unique(Ti_ngb(:));
                
                % Get barycentric co-ords
                n=numel(Ti_ngb);
                uv_i=GetPlanarBCs(MDS.triangle_data,Ti_ngb,ones(n,1)*Pg(i,:));
                idx_val=uv_i(:,1)>-tol2 & uv_i(:,1)<(1+tol2) & ...
                    uv_i(:,2)>-tol2 & uv_i(:,2)<(1+tol2) & ...
                    (uv_i(:,1)+uv_i(:,2))<(1+tol2);
                idx_val=find(idx_val,1);
                
            end
            
        end
        
        % Update Tuv
        if ~isempty(idx_val)
            tuv=[Ti_ngb(idx_val) uv_i(idx_val,:)];
        else
            tuv=[NaN NaN NaN];
        end            
        Tuv(idx(i),:)=tuv;
        
        
    end
    
end

id_nan=isnan(Tuv(:,1));
if any(id_nan)
    [~,tuv]=SphericalTriangleIntersection(TR,D(id_nan,:),MDS.triangle_data.a_thr);
    Tuv(id_nan,:)=tuv;
end

uv=Tuv(:,2:3);
uv(uv<0)=0;
uv(uv>1)=1;
uv_net=sum(uv,2);
idx=uv_net>1;
if any(idx)
    uv(idx,:)=bsxfun(@rdivide,uv(idx,:),uv_net(idx));
end
Tuv(:,2:3)=uv;

% Spherical barycentric coordinates
if nargout>2
    idx_val=~isnan(Tuv(:,1));
    if any(idx_val)
        uvw(idx_val,:)=GetSphericalBCs(MDS.triangle_data,Tuv(idx_val,1),D(idx_val,:));
    end
end


function uvw=GetSphericalBCs(TD,T,P)
% Spherical barycentric coordinates of the points of intersection

A=TD.A(:,:,T);
detA=TD.detA(T);
c12=TD.c12(:,:,T);

P=permute(P,[2 3 1]);

detAu=sum(cross(P,A(:,2,:),1).*A(:,3,:),1);
detAv=sum(cross(A(:,1,:),P,1).*A(:,3,:),1);
detAw=sum(c12.*P,1);

uvw=bsxfun(@rdivide,[detAu(:) detAv(:) detAw(:)],detA(:));


function uv=GetPlanarBCs(TD,T,P)
% Planar barycentric coordinates of the points of intersection

det_M=-sum(TD.cAB(T,:).*P,2);

u=-sum(TD.cCB(T,:).*P,2);
u=u./det_M;

v=-sum(TD.cAC(T,:).*P,2);
v=v./det_M;
uv=[u v];


function TD=PrecomputeTriangleData(Tri,X)

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
TD.c12=cross(A(:,1,:),A(:,2,:),1);

% Largest separation between two vertices, used as a back-up in case 
% 'SphericalTriangleIntersection' has to be called
TD.a_thr=1.01*AngleThreshold({Tri X});


function CD=BreakMeshIntoCharts(Tri,X)

flag=CheckFoldovers(Tri,X);
if flag
    error('Input mesh contains foldovers. Unable to proceed.')
end

% Chart generators
x=[1 0 0;-1 0 0];
y=circshift(x,[0 1]);
z=circshift(y,[0 1]);
G=cat(1,x,y,z);

% Assign mesh vertices to generators 
Dx=bsxfun(@minus,permute(X,[1 3 2]),permute(G,[3 1 2]));
Dx=sum(Dx.^2,3); 
[~,G_id_x]=min(Dx,[],2);
G_id_x=G_id_x(Tri);

% Assign triangle centroids to generators 
C=X(Tri(:,1),:,1)+X(Tri(:,2),:,1)+X(Tri(:,3),:,1);
C=ProjectOnSn(C);
Dc=bsxfun(@minus,permute(C,[1 3 2]),permute(G,[3 1 2]));
Dc=sum(Dc.^2,3);
[~,G_id_c]=min(Dc,[],2);

% Face-vertex adjacency matrix  
[FV,VA]=FaceVertexMat({Tri X});

%figure('color','w')
%axis equal off, hold on
%h=trimesh(triangulation(Tri,X));
%set(h,'EdgeColor','k','FaceColor',[0.9 0.9 0.9],'EdgeAlpha',0.5);

% Partition mesh into overlapping charts
[F_sub,V_sub,TR_sub,LUT]=deal(cell(6,1));
sc=[1 -1 -1 1 1 -1]; % sign correction
warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
warning('off','MATLAB:delaunayTriangulation:ConsConsSplitWarnId')
gen_val=true(6,1);
for g=1:6
    
    % Assign mesh triangles to generators    
    id_f=sum(G_id_x==g,2)>0 | G_id_c==g;
    if ~any(id_f), continue; end 
    
    % Expand neighbourhood of the chart 
    id_v=sum(FV(:,id_f),2)>0;
    id_f=full(sum(FV(id_v,:),1)>0)';
    
    % List of faces
    F_sub{g}=find(id_f);  
    
    % Remove non-referenced vertices
    [tr,V_sub{g}]=RemoveNonRefVerts({Tri(id_f,:) X});
    
    % Stereographic projection
    tri=tr{1};
    x=StereoProj(tr{2},g);

    % Remove inverted triangles
    x21=x(tri(:,2),:)-x(tri(:,1),:);
    x31=x(tri(:,3),:)-x(tri(:,1),:);
    s=sign(x21(:,1).*x31(:,2) - x21(:,2).*x31(:,1));
    s=sc(g)*s; 
    s=s<0;
    if sum(s)>0
        %Xg=tri(s,:);
        %Xg=unique(Xg(:));
        %Xg=tr{2}(Xg,:);
        tri(s,:)=[];
    end
    tr=triangulation(tri,x);
    
           
    % Convert tr, a 'triangulation' object, to a 'delaunayTriangulation' 
    % object. This is done for two reasons: 
    % 1) Older versions of Matlab (<R2015a) do not support 'pointLocation' 
    %    with 'triangulation' objects, and
    % 2) 'pointLocation' works much faster when called with 
    %    'delaunayTriangulation' objects than 'triangulation' objects
    DT=delaunayTriangulation(tr.Points,edges(tr));
    TR_sub{g}=DT;
    
    % Mapping between DT and tr faces; some DT faces will have no counterparts in tr
    [LUT{g},flag]=dt2tr_map(DT,tr);        
    if flag
        % Unable to generate reliable charts due to self-intersections caused by stereographic projection.
        % g-th chart will be ignored.
        gen_val(g)=false;
    end
    
    %disp(g)
    if g==1 && false
        plot3(G(g,1),G(g,2),G(g,3),'.b','MarkerSize',25)
        [f,x]=GetMeshData(tr);
        x=cat(2,ones(size(x,1),1),x);        
        h=trimesh(triangulation(f,x));
        set(h,'EdgeColor','r','FaceColor','w','EdgeAlpha',0.5,'FaceAlpha',0.5);
        %v=X(V_sub{g},:);
        %plot3(v(:,1),v(:,2),v(:,3),'.r','MarkerSize',20)
        
        if ~isempty(Xg)
            plot3(Xg(:,1),Xg(:,2),Xg(:,3),'.g','MarkerSize',20)
            view([90 0])
        end
    end
    %disp('------------------')
    
end
warning('on','MATLAB:triangulation:PtsNotInTriWarnId')
warning('on','MATLAB:delaunayTriangulation:ConsConsSplitWarnId')

if sum(gen_val)<6
    idx=~gen_val;
    G(idx,:)=10;
end

CD.generators=G;
CD.charts=TR_sub;
CD.F_sub=F_sub;
CD.V_sub=V_sub;
CD.LUT=LUT;
CD.VA=VA;


function [LUT,flag,X_new]=dt2tr_map(DT,TR)
% Map between corresponding triangles

F_dt=DT.ConnectivityList;
F_tr=TR.ConnectivityList;

F_dt=sort(F_dt,2);
F_tr=sort(F_tr,2);

LUT=nan(size(F_dt,1),1);
[~,ia,ib]=intersect(F_dt,F_tr,'rows','stable');
LUT(ia)=ib;

if numel(ia)~=size(F_tr,1)
    flag=true; % one or more triangles have been split due to self-intersections caused by stereogrphic projection
else
    flag=false;
end

if nargout<3, return; end

X_dt=DT.Points;
X_tr=TR.Points;
if ~isequal(X_dt,X_tr)
    Nx_tr=size(X_tr,1);
    X_new=X_dt;
    X_new(1:Nx_tr,:)=[];
else
   X_new=[];
end


function x=StereoProj(X,G_id)
% Stereographic projection from unit sphere onto a plane tangent to 
% generator G(G_id,:)

switch G_id
    case 1 % x=1
        Y=X(:,[2 3]);
        D=1+X(:,1);
    case 2 % x=-1
        Y=X(:,[2 3]);
        D=1-X(:,1);
    case 3 % y=1
        Y=X(:,[1 3]);
        D=1+X(:,2);
    case 4 % y=-1
        Y=X(:,[1 3]);
        D=1-X(:,2);
    case 5 % z=1
        Y=X(:,[1 2]);
        D=1+X(:,3);
    case 6 % z=-1
        Y=X(:,[1 2]);
        D=1-X(:,3);
end

x=bsxfun(@rdivide,Y,D);


function [flag,idx]=CheckFoldovers(Tri,X)
% Determine if spherical mesh contains inverted faces

FN=TriangleNormals({Tri X});

C=X(Tri(:,1),:)+X(Tri(:,2),:)+X(Tri(:,3),:);
C=ProjectOnSn(C); % face centroids

idx=sum(FN.*C,2)<=cos(80*pi/180); % find inverted faces
flag=sum(idx)>0;

