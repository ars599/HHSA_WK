function [IMF] = dt_EEMD_D(data,TNM,NoiseLevel,NE,maxcp)
% Input
% data: 1D data
% fs: frequency rate
% TNM: number of IMF;
% ifmethod:first layer IMF, method of determining an instantaneous frequency (and amplitude)
% Output
% fm:   - 2-D matrix that specifies the frequency values
% am:   - 2-D matrix that specifies the amplitude values
% IMF:  - first EMD 2D, IMF <L,TNM> L: data length,  TNM : number of IMF;
%----- Define default parameters
if nargin<2
    TNM = []; % TNM : number of required imf; if it is less than zero, then automatically determine the number
end


if nargin<3
    NoiseLevel = [];
end
if nargin<4
    NE = [];
end
if nargin<5
    maxcp = [];
end
if isempty(TNM)
    TNM = -1;
end

if isempty(maxcp)
    maxcp=5;
end

if (TNM <= 0) % automatic estimating number of imf
    TNM=fix(log2(length(data)));
end

% EMD Parameter
if isempty(NoiseLevel)
    NoiseLevel = 0;  % NoiseLevel : level of added noise
end
if isempty(NE)
    NE = 1;         % NE : number of ensemble
end 

toFlip = 0;         % toFlip : 0=> Original EEMD, References[2] ; 1=> Add anti-phase noise into signal, the way is the same as CEEMD, References[3]
numIteration = 10;  % numIteration : number of sifting iteration
typeSpline = 3;     % typeSpline : 1=> clamped spline; 2=> not a knot spline;
toModify = 1;       % toModify : 0=> None ; 1=> Apply modified linear extrapolation to boundary ; 2 => Mirror Boundary
randType = 1;       % randType : 1=> uniformly distributed white noise; 2=> gaussian white noise
%        seedNo : random seed used for white noise; The value of seed must be an integer between 0 and 2^32 - 1
%        checkSignal : 1=> verify if input signal had NaN or Infinity elements; 0=> Not verify input signal
% EEMD
%IMF = rcada_eemd_scn(data,NoiseLevel,NE,TNM,toFlip,numIteration); % EMD
L =length(data);
sd=std(data,1);
warning('off', 'all');
IMF = rcada_eemd_scn(data,NoiseLevel,NE,TNM); % EMD
%[fm,am] = fa(IMF,1/fs,ifmethod);
cpoint =[];
%IMF=IMF';
for i_imf = 1:TNM
    [indmin, indmax,indzc] = extr(IMF(:,i_imf));
    if ~isempty(indmax) && ~isempty(indmin) && sum(IMF(indmax,i_imf))>(1e-10*sd)
        if length(indmax)+length(indmin)+length(indzc)>=maxcp
            cpoint=[cpoint i_imf];
        end
    else
        break
    end
end
cpoint=1:max(cpoint);
CL = length(cpoint);
% if isequal(ifmethod, 'qzc')
%    fm(:,1:CL)= FAqzc(IMF(:,1:CL),1/fs);
%    am(:,1:CL) = envelope(IMF(:,1:CL));
% elseif isequal(ifmethod, 'qzc2')
%    fm(:,1:CL)= FAqzc(IMF(:,1:CL),1/fs);
%    am(:,1:CL) = envelope(IMF(:,1:CL));
% else
%    am(:,1:CL) = envelope(IMF(:,1:CL));
%    fm(:,1:CL) = fa(IMF(:,1:CL),1/fs,ifmethod);
% end


if TNM>CL  
    IMF = [IMF(:,1:CL) sum(IMF(:,CL+1:end),2)];
%     am(:,CL+1) = abs(IMF(:,end));
%     fm(:,CL+1) = zeros(L,1);
end


%%
function [indmin, indmax, indzer] = extr(x,t)
%extracts the indices corresponding to extrema

if(nargin==1)
    t=1:length(x);
end
if size(x,1)>size(x,2)
    x=x';
end
m = length(x);

if nargout > 2
    x1=x(1:m-1);
    x2=x(2:m);
    indzer = find(x1.*x2<0);
    
    if any(x == 0)
        iz = find( x==0 );
        indz = [];
        if any(diff(iz)==1)
            zer = x == 0;
            dz = diff([0 zer 0]);
            debz = find(dz == 1);
            finz = find(dz == -1)-1;
            indz = round((debz+finz)/2);
        else
            indz = iz;
        end
        indzer = sort([indzer indz]);
    end
end

d = diff(x);

n = length(d);
d1 = d(1:n-1);
d2 = d(2:n);
indmin = find(d1.*d2<0 & d1<0)+1;
indmax = find(d1.*d2<0 & d1>0)+1;


% when two or more consecutive points have the same value we consider only one extremum in the middle of the constant area

if any(d==0)
    
    imax = [];
    imin = [];
    
    bad = (d==0);
    dd = diff([0 bad 0]);
    debs = find(dd == 1);
    fins = find(dd == -1);
    if debs(1) == 1
        if length(debs) > 1
            debs = debs(2:end);
            fins = fins(2:end);
        else
            debs = [];
            fins = [];
        end
    end
    if length(debs) > 0
        if fins(end) == m
            if length(debs) > 1
                debs = debs(1:(end-1));
                fins = fins(1:(end-1));
                
            else
                debs = [];
                fins = [];
            end
        end
    end
    lc = length(debs);
    if lc > 0
        for k = 1:lc
            if d(debs(k)-1) > 0
                if d(fins(k)) < 0
                    imax = [imax round((fins(k)+debs(k))/2)];
                end
            else
                if d(fins(k)) > 0
                    imin = [imin round((fins(k)+debs(k))/2)];
                end
            end
        end
    end
    
    if length(imax) > 0
        indmax = sort([indmax imax]);
    end
    
    if length(imin) > 0
        indmin = sort([indmin imin]);
    end
    
end

