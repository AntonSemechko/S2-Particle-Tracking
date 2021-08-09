function s2_particle_tracking_demo(alg,spl)
% Demo illusrating one potential use of S^2 paraticle tracking functions
%   (a) SphericalTriangleIntersection_UsingStereoProj
%   (b) SphericalTriangleIntersection 
%
% INPUT:
%   - alg    : (optional) set alg=1 to use (a) and set ald=2 to use (b) for
%              particle tracking. alg=1 is the default setting.
%   - spl    : (optional) integer in the range [1,100] specifying index of 
%              the sample velocity field. Default setting is to select
%              sample field randomly.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Which particle tracking algo to use?
if nargin<1 || isempty(alg)
    alg=1;
elseif numel(alg)~=1 || ~isnumeric(alg) || ~ismember(alg,[1 2])
    error("First input argument ('alg') must be either 1 or 2.")
end

% Index of sample velocity field
if nargin<2 || isempty(spl)
    spl=randi(100);
elseif numel(spl)~=1 || ~isnumeric(spl) || spl<1 || spl>100 || spl~=round(spl)
    error("Second input argument ('spl') must be an integer in the range [1,100].")
end
fprintf("Running demo with sample velocity field #%u\n",spl)


% Load sample mesh and tangential velocity field
persistent TR VF

if isempty(TR)
    try
        TR=load('sample_velocity_fields.mat');
        [TR,VF]=deal(TR.TR,TR.VF);
    catch
        error("File 'sample_velocity_fields.mat' is not on path. Unable to run demo wihtout it.")
    end
end

[Tri,X]=GetMeshData(TR); % mesh faces and vertices
Vt=VF(:,:,spl); % tangential velocity field defined at the mesh vertices 


% Initialize data structure that will contains a bunch of pre-computed
% variables to speed-up particle tracking
if alg==1 
    MDS=[]; % for (a)
else
    BS=BinSphericalTriangles(TR); % for (b)
end


% Convert tangential velocities to angular velocities for the purpose of
% intergration
W=cross(X,Vt,2);

% Visualize the velocity field
hf=figure('color','w');
if isprop(hf,'WindowState'), hf.WindowState='maximized'; end

hp=trimesh(SubdivideSphericalMesh(IcosahedronMesh,5));
set(hp,'FaceColor',0.8*[1 1 1],'EdgeColor','none','FaceAlpha',0.8);
hold on
axis equal off

quiver3(X(:,1),X(:,2),X(:,3),Vt(:,1,1),Vt(:,2,1),Vt(:,3,1),'-r','LineWidth',1);
drawnow
pause(0.5)


% Initialize particles 
Np=642;
P=X(1:Np,:);

ro=1+1E-3;
[xD,yD,zD]=deal(P(:,1)',P(:,2)',P(:,3)');

nan_pad=nan(1,Np);
h0=plot3(ro*P(:,1),ro*P(:,2),ro*P(:,3),'.b','MarkerSize',12);
h1=plot3(reshape(cat(1,ro*xD,nan_pad),[],1),reshape(cat(1,ro*yD,nan_pad),[],1),reshape(cat(1,ro*zD,nan_pad),[],1),'-k','LineWidth',1);

% Camera view and lighting 
ha=gca;
ha.CameraPosition= [-9.1422 -11.9 8.7179];
ha.CameraTarget= [0.0071631 0.023905 0.04066];
ha.CameraUpVector= [0 0 1];
ha.CameraViewAngle= 5.838;

set(hp,'SpecularExponent',35,'SpecularStrength',0.15)
hl1=camlight('headlight');
set(hl1,'style','infinite','position',10*get(hl1,'position'))
hl2=light('position',-get(hl1,'position'));
set(hl2,'style','infinite')
lighting phong


%pause(1)
%im=export_fig('-nocrop','-r150','-silent');
%imwrite(im(:,701:2300,:),'img_00.jpg');
%[A,map] = rgb2ind(im(:,701:2300,:),256);
%imwrite(A,map,'s2_particle_tracking_demo.gif','gif','LoopCount',Inf,'DelayTime',0.1);


% Integrate particle positions forward in time  
[Tuv,tri,r]=deal(zeros(Np,3));
[u,v]=deal(zeros(Np,1));
R=zeros(3,3,Np);
ds=1*pi/180;
cnt=0;

m=200; % # of most recent particle positions to show
k=0;
while cnt<1E5

    cnt=cnt+1;
    
    % Localize the partilces
    if alg==1
        [Tuv(:),MDS]=SphericalTriangleIntersection_UsingStereoProj(TR,P,MDS);
    else
        [~,Tuv(:)]=SphericalTriangleIntersection(TR,P,[],BS);
    end        
    
    % Interpolate particle velocities
    [u(:),v(:)]=deal(Tuv(:,2),Tuv(:,3));
    tri(:)=Tri(Tuv(:,1),:);    
    r(:)=(1-u-v).*W(tri(:,1),:) + u.*W(tri(:,2),:) + v.*W(tri(:,3),:);
    r=ProjectOnSn(r);
    
    % Update particle positions
    R(:)=r2R(ds*r);
    P(:)=parallel_mat_multiply(R,P')';
    if mod(cnt,20)==0
        P=ProjectOnSn(P);
    end

    % Visualize
    if isvalid(hf)
        
        if cnt<m
            xD=cat(1,xD,P(:,1)');
            yD=cat(1,yD,P(:,2)');
            zD=cat(1,zD,P(:,3)');
        else
            xD=cat(1,xD(2:m,:),P(:,1)');
            yD=cat(1,yD(2:m,:),P(:,2)');
            zD=cat(1,zD(2:m,:),P(:,3)');
        end                
        
        set(h1,'XData',reshape(cat(1,ro*xD,nan_pad),[],1),...
               'YData',reshape(cat(1,ro*yD,nan_pad),[],1),...
               'ZData',reshape(cat(1,ro*zD,nan_pad),[],1))

        set(h0,'XData',ro*P(:,1),...
               'YData',ro*P(:,2),...
               'ZData',ro*P(:,3))
           
        if mod(cnt,max(50,m))==0
            drawnow
        end
        pause(1E-3)

        
        % grab a gif frame
        if mod(cnt,5)==0 && k<30 && false
            k=k+1;
            
            pause(1)
            im=export_fig('-nocrop','-r150','-silent');
            if k<10
                imwrite(im(:,701:2300,:),sprintf('img_0%u.jpg',k));
            elseif k<100
                imwrite(im(:,701:2300,:),sprintf('img_%u.jpg',k));
            end
            
            [A,map] = rgb2ind(im(:,701:2300,:),256);
            imwrite(A,map,'s2_particle_tracking_demo.gif','gif','WriteMode','append','DelayTime',0.1);
        end
        
        
    else
        break
    end    

end



function R=r2R(r,A)
% Exponential map from so(3) to SO(3) 

N=size(r,1);
a2=zeros(1,1,N);
a2(:)=sum(r.^2,2); % rotation angle (radians^2)
a2(:)=max(a2,1E-12);
a=sqrt(a2);

S=sinc(a/pi);
C=(1-cos(a))./a2;

% Angular velocity tensor
A(1,2,:)=-r(:,3);
A(1,3,:)= r(:,2);
A(2,1,:)= r(:,3);
A(2,3,:)=-r(:,1);
A(3,1,:)=-r(:,2);
A(3,2,:)= r(:,1);

% Rotation matrix
R=bsxfun(@times,S,A) + bsxfun(@times,C,parallel_mat_multiply(A));
R(1,1,:)=1 + R(1,1,:);
R(2,2,:)=1 + R(2,2,:);
R(3,3,:)=1 + R(3,3,:);
