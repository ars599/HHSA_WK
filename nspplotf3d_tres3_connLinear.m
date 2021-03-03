function [plv,plph]=nspplotf3d_tres3_connLinear(ph, phr, tres, half_winsize, S)

% Input-
%	fm	    - 2-D matrix that specifies the frequency values
%	FM      - 3-D matrix that specifies the frequency values of envelope 
%	PH      - 3-D matrix that specifies the phases of envelope 

%   tres    - the time resolution
% Output-

%	All_nt	- 4-D matrix of the HHS spectrum, where
%		    1st dimension specifies the number of frequencies,
%		    2nd dimension specifies the number of frequencies of envelope
%           3rd dimension specifies the number of time points
%           4th dimension what number IMF
%	fscale	- vector that specifies the frequency-axis values
%	Fscale	- vector that specifies the AM frequency-axis values


%----- Check the input arguments
if nargin<3
    error('nspplot: both frequency and amplitude matrices required');
end

%----- Specify default variables
% if nargin < 5
%     ntp0 = [];
%     ntp1 = [];
% end

if nargin<3
    tres=[];
end
if nargin<4
    half_winsize=[];
end
if nargin<5
    S=[];
end
[npt,nimf]=size(ph);

%----- Initialize default variables
% if isempty(ntp0)
%     ntp0=1;
%     ntp1=npt;
% end

if isempty(tres)
   tres=ceil(npt/25);
end

if isempty(half_winsize)
   half_winsize=150;
end
%winsize=2*half_winsize+1;

t=ceil((1:npt)*tres/npt); %t is the mapping position of time values into the time axis grid
[~,firstid]=unique(t);
intv=round(mean(diff(firstid)));
if mod(intv,2)==1
    dist_to_med=round(0.5*(intv-1));
else
    dist_to_med=round(0.5*intv);
end
medid=(firstid+dist_to_med)';

%All_nt=zeros(tres,1);

%for i_imf = 1:nimf
    %nt=zeros(fres,Freq);
%     f = fm(:,1);
    p = ph(:,1);
%     fr = fmr(:,1);
    pr = phr(:,1);
%     medf=median(f);
%     medF=median(F);
%     [~,amid]=min(abs(medF-medf));
    %F_am=F(:,amid);
%     P_am=P(:,amid);
    ph_diff=pr-p;
plv_vec= zeros(1,length(medid));  
plv_vec=complex(plv_vec);
for ires= 1:length(medid)
    if (medid(ires)- half_winsize)>=1
        ph_initid=medid(ires)- half_winsize;
    else
        ph_initid=1;
    end
    if (medid(ires)+ half_winsize)<=npt
        ph_endid=medid(ires)+ half_winsize;
    else
        ph_endid=npt;
    end
    ph_series=ph_diff(ph_initid:ph_endid);
    ph_vecs = exp(sqrt(-1)*ph_series);
    plv_vec(ires)=mean(ph_vecs);
end
plv=abs(plv_vec);
plph=angle(plv_vec);    

