function [Ta,Tb,Ta_oh,Tb_oh,Nt,Nv]=s2_ray_triangle_intersection_benchmark_test(Nt,Nr,n_rep,vis)
% Run a Monte Carlo simulation to evaluate performance of two functions
%   (A) 'SphericalTriangleIntersection_UsingSttereoProj.m', and
%   (B) 'SphericalTriangleIntersection.m'
% used for performing ray/triangle intersection tests on surface meshes 
% of a unit sphere. Here, both the rays and mesh vertices are drawn 
% from uniform distributions.
%
% INPUT:
%   - Nt    : 1-by-M vector specifying complexity of the test meshes in 
%             terms of number of triangles, where 1<=M<=20. Every entry 
%             in Nt must be a positive integer >= 100. Default setting is
%             Nt=[1 1 2.5 5 10 25 50 100 250 500 1E3]*1E3;
%   - Nr    : 1-by-N vector specifying number of simultaneous ray/triangle
%             intersection tests per mesh. Every entry in Nr must a 
%             positive integer >=1. Default setting is 
%             Nr=[0.1 0.5 1 2.5 5 10]*1E3;       
%   - n_rep : positive integer specifying the number of tests to perform 
%             for each unique combination of two entries from Nt and Nr. 
%             Total number of tests will be M*N*n_rep. n_rep=3 is the
%             defautl setting.
%   - vis   : set vis=false to disable plotting. vis=true is the default
%             setting.
%
% OUTPUT:
%   - Ta    : M-by-N-by-n_rep array containing total run-times (in seconds)
%             for function (A), so that Ta(i,j,k) is the run-time 
%             corresponding to Nt(i) and Nr(j) at repetition 1<=k<=n_rep
%   - Tb    : M-by-N-by-n_rep array containing total run-times (in seconds)
%             for function (B).
%   - Ta_oh : 1-by-M vector of "overhead" times used to initialize data 
%             structure used by function (A) to speed-up ray/triangles 
%             intersection tests.   
%   - Tb_oh : 1-by-M vector of "overhead" times used to initialize data 
%             structure used by function (B) to speed-up ray/triangles 
%             intersection tests.   
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Experiment settings
% -------------------------------------------------------------------------
if nargin<1 || isempty(Nt)
    Nt=[1 1 2.5 5 10 25 50 100 250 500 1E3]*1E3; % number of mesh triangles
elseif ~isvector(Nt) || any(Nt<100 | ~isfinite(Nt) | round(Nt)~=Nt)
    error('Invalid entry for 1st input argument (Nt)')
end
Nt=unique(Nt);
    
if nargin<2 || isempty(Nr)
    Nr=[0.1 0.5 1 2.5 5 10 25 50 100]*1E3;       % number of rays
elseif ~isvector(Nr) || any(Nr<1 | ~isfinite(Nr) | round(Nr)~=Nr)
    error('Invalid entry for 2nd input argument (Nr)')
end

if nargin<3 || isempty(n_rep)
    n_rep=3;                                     % number of repetitions for each 2-combination of Nv and Nr
elseif numel(n_rep)~=1 || ~isnumeric(n_rep) || ~isfinite(n_rep) 
    error('Invalid entry for 3rd input argument (n_rep)')
end

if nargin<4 || isempty(vis)
    vis=true;                                    % visualization on 
elseif numel(vis)~=1 || ~islogical(vis) 
    error('Invalid entry for 4th input argument (vis)')
end


% Begin visualization
% -------------------------------------------------------------------------
hf=NaN;
if vis
    try
        c_map=cbrewer_accent_colormap;
        
        hf=figure('color','w');
        ha=gca;
        
        title('A - stereo projection (○) | B - spherical binning (☐)','FontSize',25,'FontWeight','bold')
        xlabel('mesh complexity (# of triangles)','FontSize',25)
        ylabel('search time per ray (msec)','FontSize',25)
        
        set(ha,'FontSize',20,'XLim',[0.9*Nt(1) 1.1*Nt(end)],'XScale','log','XLimMode','manual','YScale','log','Box','on')        
        if isprop(ha,'Toolbar')
            disableDefaultInteractivity(ha);
            ha.Toolbar.Visible=false;
            hf.WindowState='maximized';
        end
        
        hold(ha,'on');
        
        an_1a = animatedline(ha);
        set(an_1a, 'Linestyle','-','LineWidth',2,'Color','k');
        
        an_1b = animatedline(ha);
        set(an_1b, 'Linestyle','--','LineWidth',2,'Color','k');
                
        
    catch
        vis=false;
    end
end




Nr_str=cell(length(Nr),1);
for i=1:length(Nr)
    Nr_str{i}=sprintf('%.1E',Nr(i));
end


% Run the experiment
% -------------------------------------------------------------------------
Nv=ceil((Nt+4)/2);                                % corresponding number of vertices
N=n_rep*length(Nv)*length(Nr);                    % total number of tests
[Ta,Tb]=deal(zeros(length(Nv),length(Nr),n_rep)); % total run-time per test
[Ta_oh,Tb_oh]=deal(zeros(1,length(Nv)));
for nv=1:length(Nv) % number of mesh vertices
    
    % Generate mesh of a unit sphere with Nv(nv) randomly distributed vertices
    X=RandSampleSphere(Nv(nv),'stratified');
    Tri=convhull(X);
    if ClosedMeshVolume({Tri X})<0
        Tri=fliplr(Tri); % make sure faces have counter-clockwise orientation
    end 
    TR=triangulation(Tri,X);
    
    % Pre-compute data structure used by 'SphericalTriangleIntersection_UsingStereoProj'; this structure can be omitted if function is being called only once
    tic 
    [~,MDS]=SphericalTriangleIntersection_UsingStereoProj(TR);
    Ta_oh(nv)=toc;
    
    % Pre-compute data structure used by 'SphericalTriangleIntersection'; ; this structure can be omitted if function is being called only once
    n_lat=20; % number of latitude bins
    n_lon=35; % number of longitude bins
    adp=1;    % set adp=1 to automatically adjust bin sizes to data 
    tic 
    BS=BinSphericalTriangles(TR,[n_lat n_lon adp]);
    Tb_oh(nv)=toc;
    
    % Do ray/trianle intersection tests using varying number of rays
    if ~vis
        
        if nv==1
            fprintf(2,'\nLegend:\n')
            fprintf("Function A = SphericalTriangleIntersection_UsingStereoProj.m (uses stereographic projection to speed-up search for the ray-triangle intersection pairs)\n");
            fprintf("Function B = SphericalTriangleIntersection.m (uses binning of longitude & latitude to speed-up search for the ray-triangle intersection pairs)\n");
        end
        
        fprintf(2,'\nExperiment %2d : # of mesh faces = %d\n',nv,Nt(nv));
        fprintf('--------------+--------------+--------------\n')
        fprintf('# rays        |  Function A  |  Function B\n');
        fprintf('--------------+--------------+--------------\n')
        
    end
        
    for nr=1:length(Nr) % number of rays
                    
        for n=1:n_rep

            if vis, fprintf('.'); end
            
            % Generate random set of rays emanating from the origin (these can be thought as points on unit sphere)
            P=RandSampleSphere(Nr(nr),'uniform');

            % Find spherical triangles intersected by the rays
            tic
            [~,~]=SphericalTriangleIntersection_UsingStereoProj(TR,P,MDS);
            Ta(nv,nr,n)=toc;
            
            tic
            [~,~]=SphericalTriangleIntersection(TR,P,[],BS);
            Tb(nv,nr,n)=toc;
                        
            N=N-1;
        end
        if vis, fprintf(' '); end
        
        % Average search times per ray 
        vis=vis & ishandle(hf);
        if vis 
            
            if ~ishandle(hf), return; end

            t1=Ta(nv,nr,:)/Nr(nr);
            t1=mean(t1(:));
            hl_1=plot(ha,Nt(nv),1E3*t1,'o','MarkerSize',10,'MarkerFaceColor',c_map(nr,:),'MarkerEdgeColor','k');
            
            t2=Tb(nv,nr,:)/Nr(nr);
            t2=mean(t2(:));
            hl_2=plot(ha,Nt(nv),1E3*t2,'s','MarkerSize',10,'MarkerFaceColor',c_map(nr,:),'MarkerEdgeColor','k'); %#ok<*NASGU>
 
            if nr==length(Nr)
                addpoints(an_1a,Nt(nv),1E3*t1);
                addpoints(an_1b,Nt(nv),1E3*t2);
                %disp(t2/t1)
            end            
            
            if nr==1
                hl=hl_1;
            else
                hl=[hl;hl_1]; %#ok<*AGROW>
            end
            
            pause(0.1)
            drawnow
        else
            t2=Tb(nv,nr,:)/Nr(nr);
            t1=Ta(nv,nr,:)/Nr(nr);            
            fprintf('%-13.3E |  %-12.3E|  %-12.3E(msec/ray)\n',Nr(nr),1E3*mean(t1(:)),1E3*mean(t2(:)));
        end
        
    end
    
    if N>0
        if vis
            fprintf('| %u tests remaining\n',N);            
        end
    else
        fprintf('\n');
    end
    
    vis=vis & ishandle(hf);
    if nv==1 && vis
        hl=legend(ha,hl(:),Nr_str,'Location','EastOutside');
        try
            hl.AutoUpdate=false;
            hl.Title='# rays / function call';
        catch
        end
    end
    
end

if nargout<1, clear T1 T2 Ta_oh Tb_oh Nt Nv; end


function map=cbrewer_accent_colormap
% Color map used to distinguish among run-times with varying numbers of rays

map=[127,201,127;...
     144,185,180;...
     166,177,205;...
     190,174,212;...
     228,180,173;...
     253,192,134;...
     254,214,137;...
     255,241,145;...
     255,255,153;...
     203,226,162;...
     108,164,172;...
     56,108,176;...
     104,62,169;...
     192,20,151;...
     240,2,127;...
     222,45,69;...
     191,91,23;...
     163,97,27;...
     134,101,49;...
     102,102,102]/255;

