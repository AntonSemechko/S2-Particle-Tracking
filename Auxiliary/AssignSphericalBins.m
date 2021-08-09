function Bidx=AssignSphericalBins(n,X)
% Assign spherical bin indices to points on a unit sphere.
%
% INPUT:
%   - n     : 1-by-2 cell containing histogram bin edges. Latitude and 
%            logitude edges must be contained in n{1} and n{2}, 
%            respectively. n{3} is optional and specifies whether bin
%            edges are distributed uniformly (i.e., n{3}==0) or not
%            (i.e., n{3}=1);
%   - X     : N-by-3 array of spherical point coordinates.
%
% OUTPUT:
%   - Bidx  : N-by-1 array of bin indices
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Bin edges
E_lat=n{1};
E_lon=n{2};
chk_uni=false;
if numel(n)==3
    chk_uni=n{3}==0;
    if chk_uni
        dt_lat=E_lat(2)-E_lat(1);
        dt_lon=E_lon(2)-E_lon(1);
    end
end
n=[numel(E_lat) numel(E_lon)]-1; % # of bins along latitude and longitude

% Covert X from Cartesian to spherical co-ords. Convention used here:  
%   latitude  range [0,  pi]    : lat=0 is at the north pole
%   longitude range [0,2*pi]    : lon=0 is along the meridian passing through (1,0,0)

lat=EvalGeoDist(X,[0 0 1]);
lon=atan2(X(:,2),X(:,1));

idx=lon<0;
lon(idx)=2*pi+lon(idx); % unwrap from [-pi,pi] to [0,2*pi] interval
%lon(lon<0)=0;
%lon(lon>(2*pi))=2*pi;

% Bin the vertices
if chk_uni % uniform grid

    idx_lat=floor(lat/dt_lat)+1;
    idx_lon=floor(lon/dt_lon)+1;    
    idx_lat=min(idx_lat,n(1));
    idx_lon=min(idx_lon,n(2));

else % non-uniform grid    

    % First implementation; comparable to the one below that uses 'histc'
    % -------------------------
    %idx_lat=BinValues(lat,E_lat);
    %idx_lon=BinValues(lon,E_lon);
    % -------------------------
    
    
    % This seems to work reasonably well
    % ----------------------------------
    [~,idx_lat]=histc(lat,E_lat(:)); %#ok<*HISTC>
    [~,idx_lon]=histc(lon,E_lon(:));
    idx_lat=max(min(idx_lat,n(1)),1);
    idx_lon=max(min(idx_lon,n(2)),1);
    % ----------------------------------
    
    
    % Antother implementation; doesn't work for Matlab versions <= R2014a
    % ---------------------------------------------------------------------
    %flag=false;
    %try 
    %    idx_lat=discretize(lat,E_lat(:));
    %    idx_lon=discretize(lon,E_lon(:));
    %    if sum(isnan([id_lat;id_lon]))>0
    %        flag=true;
    %    end
    %catch
    %    flag=true;
    %end
    
    % 'discretize' function is either unavailable (i.e., Matlab version 
    % 2014 or older) OR returned NaNs due to limited precision
    %if flag
    %    [~,~,idx_lat]=histcounts(lat,E_lat(:));
    %    [~,~,idx_lon]=histcounts(lon,E_lon(:));
    %    idx_lat=max(min(idx_lat,n(1)),1);
    %    idx_lon=max(min(idx_lon,n(2)),1);
    %end
    % ---------------------------------------------------------------------
        
end

% Assign bin indices
Bidx=sub2ind(n(1:2),idx_lat,idx_lon);


function bin_id=BinValues(X,E)
% Assign values into non-uniform bins without using any of the built-in
% Matlab functions like 'histc', 'discretize', or 'hiscounts'

Ne=numel(E);

X=X(:)';
E=E(:);

E(1)=-Inf;
E(end)=Inf;
idx=bsxfun(@ge,X,E(1:(Ne-1))) & bsxfun(@lt,X,E(2:Ne));
[bin_id,~]=find(idx);

