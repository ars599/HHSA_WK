function [fm,am,FM,AM,IMF,IMF2,mx_nIMF2, ph, PH] = multi_EEMD_SV(data,fs,TNM,TNM2,S)
% This variation of multi_emd allows using different number of ensembles and noise level in the two level of eemd 
%[fm,am,FM,AM,IMF,IMF2,mx_nIMF2] = multi_AEEMD_SV(data,fs,TNM,TNM2,ifmethod,ifmethod2,NoiseLevel,NE,NoiseLevel2,NE2,ines,shiftLevel)
% Input
% data: 1D data
% fs: frequency ratem
% TNM: number of IMF;
% TNM2: number of imf of envelope of IMF;
% ifmethod:first layer IMF, method of determining an instantaneous frequency (and amplitude)
% ifmethod2: second layer IMF, method of determining an instantaneous frequency (and amplitude)
% Output
% fm:   - 2-D matrix that specifies the frequency values
% am:   - 2-D matrix that specifies the amplitude values
% FM:   - 3-D matrix that specifies the frequency values of envelope
% AM:   - 3-D matrix that specifies the amplitude values of envelope
% IMF:  - first EMD 2D, IMF <L,TNM> L: data length,  TNM : number of IMF;
% IMF2: - second EMD 3D, IMF2 <L,TNM2, TNM> TNM2 : number of imf of envelope of IMF;
%----- Define default parameters
if nargin<3
    TNM = []; % TNM : number of required imf; if it is less than zero, then automatically determine the number
end
if nargin<4
    TNM2 = []; % TNM : number of required imf of envelope; if it is less than zero, then automatically determine the number
end
if nargin<5
     ifmethod=[];
     ifmethod2=[];
     NoiseLevel = [];
     NE = [];
     NoiseLevel2 = [];
     NE2 = [];
     shiftLevel=[];
     phase=[];
     phase2=[];
else
     try ifmethod=S.ifmethod; catch ifmethod=[]; end
     try ifmethod2=S.ifmethod2; catch ifmethod2=[]; end
     try NoiseLevel = S.ENoise; catch NoiseLevel=[]; end
     try NE = S.NEnsemble; catch NE=[]; end
     try NoiseLevel2 = S.ENoise2; catch NoiseLevel2=[]; end
     try NE2 = S.NEnsemble2; catch NE2=[]; end
     try shiftLevel=S.shiftLevel; catch shiftLevel=[]; end
     try phase=S.phase; catch phase=[]; end
     try phase2=S.phase2; catch phase2=[]; end
end
ines=[];
if isempty(TNM)
    TNM = -1;
end
if isempty(TNM2)
    TNM2 = -1;
end
if isempty(ifmethod)
    ifmethod = 'zc';
end
if isempty(ifmethod2)
    ifmethod2 = 'zc';
end

if (TNM <= 0) % automatic estimating number of imf
    TNM=fix(log2(length(data)));
end
if (TNM2 <= 0) % automatic estimating number of imf
    TNM2=fix(log2(length(data)));
end
% EMD Parameter
if isempty(NoiseLevel)
    NoiseLevel = 0.2;  % NoiseLevel : level of added noise
end
if isempty(NE)
    NE = 50;         % NE : number of ensemble
end 
if isempty(NoiseLevel2)
    NoiseLevel2 = 0.2;  % NoiseLevel : level of added noise
end
if isempty(NE2)
    NE2 = 10;         % NE : number of ensemble
end
if isempty(phase)
    phase = 0;
end
if isempty(phase2)
    phase2 = 0;
end

if NE2>NE
    NE2=NE;
end
if isempty(shiftLevel)
    shiftLevel=0;
end
mx_nIMF2=zeros(TNM,1);
TNMs=TNM+shiftLevel;
% if shiftLevel>0
%     Y2 = spmmhh_resample(data,2^shiftLevel); % Y2: upsample data
% else
    Y2=data;
% end
L=length(Y2);
sL = (length(Y2)-1)*(2^shiftLevel)+1;
   ines = zeros(sL,TNMs,NE);
            for iii=1:NE  % imn ensemble loop    
                temp = randn(1,sL); % std temp is 1
                %imn = rcada_emd(temp, 1, 2, TNMs, 10);
                imn = dt_EEMD_D(temp,TNMs,0,1); % EMD,
                ni=size(imn,2);
                while ni< TNMs  % in this algorithm we require the noise should be scale-complete
                    temp = randn(1,sL); % std temp is 1
                    imn = dt_EEMD_D(temp,TNMs,0,1); % EMD
                    ni=size(imn,2);
                end
                %imn=cmask_emdn(temp,-1,10);
                cmn = cum_mx(imn);
                %cmn=imn;
                sdim=std(cmn,1);
                cmn=cmn./repmat(sdim,size(cmn,1),1);
                ines(:,1:size(cmn,2),iii)=cmn; 
            end
if shiftLevel>0
    ines2=ines(:,(shiftLevel+1):end,:);
    ines2=spmm_downsample(ines2,2^shiftLevel);
else
    ines2=ines;
end
% EMD Parameter
% 
% toFlip = 0;         % toFlip : 0=> Original EEMD, References[2] ; 1=> Add anti-phase noise into signal, the way is the same as CEEMD, References[3]
% numIteration = 10;  % numIteration : number of sifting iteration
% typeSpline = 3;     % typeSpline : 1=> clamped spline; 2=> not a knot spline;
% toModify = 1;       % toModify : 0=> None ; 1=> Apply modified linear extrapolation to boundary ; 2 => Mirror Boundary
% randType = 1;       % randType : 1=> uniformly distributed white noise; 2=> gaussian white noise
%        seedNo : random seed used for white noise; The value of seed must be an integer between 0 and 2^32 - 1
%        checkSignal : 1=> verify if input signal had NaN or Infinity elements; 0=> Not verify input signal
%% first layer EMD 
%IMF = rcada_eemd(data,NoiseLevel,NE,TNM,toFlip,numIteration); % EMD
%L =length(data);
sd=std(Y2,1);
% IMF = rcada_eemd(data,NoiseLevel,NE,TNM,toFlip,numIteration); % EEMD
%IMF=rcada_eiemd2C(Y2, NoiseLevel, NE, TNM, ines, shiftLevel); % Similar to Flandrin CEEMD
          IMF =rcada_eemd_scn(Y2,NoiseLevel,NE,TNM, ines, shiftLevel); % EEMD
%[fm,am] = fa(IMF,1/fs,ifmethod);
cpoint =[];
for i_imf = 1:TNM
    [indmin, indmax, indzc] = extr(IMF(:,i_imf));
    if ~isempty(indmax) && ~isempty(indmin) && sum(IMF(indmax,i_imf))>(1e-10*sd)
        if length(indmax)+length(indmin)+length(indzc)>=5
            cpoint=[cpoint i_imf];
        end
%     else
%         break  % aeemd allows an IMF near zero, but not zeros for the following IMFs
    end
end
cpoint=1:max(cpoint);
CL = length(cpoint);

%am(:,1:CL) = envlp(IMF(:,1:CL));
% if all(am(:,1)==IMF(:,1))
%     am(:,2:CL) = envlp(IMF(:,2:CL));
%     [temp, am(:,1)] = fa(IMF(:,1),1/fs,'hilbert');
% end
% if ~(isequal(ifmethod, 'qzc'))
if CL>0
    am(:,1:CL) = envelope(IMF(:,1:CL));
    if ~ phase
      try
        tmpfm = fa(IMF(:,1:CL),1/fs,ifmethod); %fm(:,1:CL)
        if ~isreal(tmpfm)
          tmpfm = fa(IMF(:,1:CL),1/fs,'zc');
        end
        fm(:,1:CL)=tmpfm;
      catch
        fm(:,1:CL)= FAqzc(IMF(:,1:CL),1/fs); 
        disp('unknown error');
      end
    else
      try
        [tmpfm,~,tmpph] = fa(IMF(:,1:CL),1/fs,ifmethod); %fm(:,1:CL)
        if ~isreal(tmpfm)
          [tmpfm,~,tmpph] = fa(IMF(:,1:CL),1/fs,'zc');
        end
        fm(:,1:CL)=tmpfm;
        ph(:,1:CL)=tmpph;
      catch
        [fm(:,1:CL),ph(:,1:CL)]= FAqzc(IMF(:,1:CL),1/fs); 
        disp('unknown error');
      end
    end
end

if TNM>CL  
    IMF = [IMF(:,1:CL) sum(IMF(:,CL+1:end),2)];
    am(:,CL+1) = abs(IMF(:,end));
    fm(:,CL+1) = zeros(L,1);
    if phase
       ph(:,CL+1) = zeros(L,1);
    end
    IMF2 = zeros(L,TNM2,CL+1);
    FM = zeros(L,TNM2,CL+1);
    if phase2
       PH = zeros(L,TNM2,CL+1); 
    end
    AM = zeros(L,TNM2,CL+1);
else
    IMF2 = zeros(L,TNM2,CL);
    FM = zeros(L,TNM2,CL);
    AM = zeros(L,TNM2,CL);
    if phase2
       PH = zeros(L,TNM2,CL); 
    end
end
%% Second layer EMD 

maxC = 0;
for i_imf=1:CL
              IMF2(:,:,i_imf)=rcada_eemd_scn(am(:,i_imf)', NoiseLevel2, NE2, TNM2, ines2,0); 
    cpoint2 =[];
    for j_imf = 1:TNM2
        [indmin, indmax,indzc] = extr(IMF2(:,j_imf,i_imf));
        if ~isempty(indmax) && ~isempty(indmin)  && sum(IMF2(indmax,j_imf,i_imf))>(1e-10*sd)
            if length(indmax)+length(indmin)+length(indzc)>=5
                cpoint2=[cpoint2 j_imf];
            end
%         else
%             break %% aeemd allows an IMF near zero, but not zeros for the following IMFs
        end
    end
    cpoint2=1:max(cpoint2);
    CL2 = length(cpoint2);
    if size(IMF2,2)>CL2
        IMF2(:,1:CL2+1,i_imf) = [IMF2(:,1:CL2,i_imf) sum(IMF2(:,CL2+1:end,i_imf),2)];
        IMF2(:,CL2+2:end,i_imf) = zeros(L,size(IMF2,2)-CL2-1);
        AM(:,CL2+1,i_imf) = abs(IMF2(:,CL2+1,i_imf)); 
        FM(:,CL2+1,i_imf) = zeros(L,1);
        mx_nIMF2(i_imf)=CL2+1;
    else
        mx_nIMF2(i_imf)=CL2;
    end
    if CL2>maxC
         maxC = CL2;
    end
    if CL2>0
        AM(:,cpoint2,i_imf) = envelope(IMF2(:,cpoint2,i_imf));
      if ~ phase2  
        try
             tmpFM=fa(IMF2(:,cpoint2,i_imf),1/fs,ifmethod2);
             if ~isreal(tmpFM)
                 tmpFM=fa(IMF2(:,cpoint2,i_imf),1/fs,'zc');
             end
             FM(:,cpoint2,i_imf) = tmpFM;
        catch
             FM(:,cpoint2,i_imf) = FAqzc(IMF2(:,cpoint2,i_imf),1/fs);
             disp('unknown error');
        end
      else
        try
             [tmpFM,~,tmpPH]=fa(IMF2(:,cpoint2,i_imf),1/fs,ifmethod2);
             if ~isreal(tmpFM)
                 [tmpFM,~,tmpPH]=fa(IMF2(:,cpoint2,i_imf),1/fs,'zc');
             end
             FM(:,cpoint2,i_imf) = tmpFM;
             PH(:,cpoint2,i_imf) = tmpPH;
        catch
             [FM(:,cpoint2,i_imf),PH(:,cpoint2,i_imf)] = FAqzc(IMF2(:,cpoint2,i_imf),1/fs);
             disp('unknown error');
        end

      end

    end
end
if TNM2>maxC
    IMF2 = IMF2(:,1:maxC+1,:);
    AM = AM(:,1:maxC+1,:);
    FM = FM(:,1:maxC+1,:);
    if phase2
       PH = PH(:,1:maxC+1,:); 
    end
end
%%
function [envmax] = envelope(data,INTERP)
%computes envelopes and mean with various interpolations

NBSYM = 2;
DEF_INTERP = 'spline';

if nargin < 2
    t = 1:length(data);
    INTERP = DEF_INTERP;
end

if ~ischar(INTERP)
    error('interp parameter must be ''linear'''', ''cubic'' or ''spline''')
end

if ~any(strcmpi(INTERP,{'linear','cubic','spline'}))
    error('interp parameter must be ''linear'''', ''cubic'' or ''spline''')
end


s = size(data);

if s(1) > s(2)
    data = data';
end
envmax=zeros(size(data));
for ijk = 1:size(data,1)
    x = data(ijk,:);
    lx = length(x);
    [indmin,indmax,~] = extr(x);
    
    if (length(indmin) + length(indmax) < 3)
      %  error('not enough extrema')
        envmax(ijk,:) = abs(x);
        %envmin(ijk,:) = -abs(x);
    else
        %boundary conditions for interpolation
        
        [tmin,tmax,xmin,xmax] = boundary_conditions(indmin,indmax,t,x,NBSYM);
        
        % definition of envelopes from interpolation
%         envmax(ijk,:) = interp1(tmax,xmax,t,INTERP);
%         envmin(ijk,:) = interp1(tmin,xmin,t,INTERP);
        [tminmax,tid]=sort([tmin tmax]);
        xm=[xmin xmax];
        xminmax=abs(xm(tid));
        envmax(ijk,:) = interp1(tminmax,xminmax,t,INTERP);
    end
end
if s(1) > s(2)
    envmax = envmax';
    %envmin = envmin';
end
%---------------------------------------------------------------------------------------

function [tmin,tmax,xmin,xmax] = boundary_conditions(indmin,indmax,t,x,nbsym)
% computes the boundary conditions for interpolation (mainly mirror symmetry)


lx = length(x);

if (length(indmin) + length(indmax) < 3)
    error('not enough extrema')
end

if indmax(1) < indmin(1)
    if x(1) > x(indmin(1))
        lmax = fliplr(indmax(2:min(end,nbsym+1)));
        lmin = fliplr(indmin(1:min(end,nbsym)));
        lsym = indmax(1);
    else
        lmax = fliplr(indmax(1:min(end,nbsym)));
        lmin = [fliplr(indmin(1:min(end,nbsym-1))),1];
        lsym = 1;
    end
else
    
    if x(1) < x(indmax(1))
        lmax = fliplr(indmax(1:min(end,nbsym)));
        lmin = fliplr(indmin(2:min(end,nbsym+1)));
        lsym = indmin(1);
    else
        lmax = [fliplr(indmax(1:min(end,nbsym-1))),1];
        lmin = fliplr(indmin(1:min(end,nbsym)));
        lsym = 1;
    end
end

if indmax(end) < indmin(end)
    if x(end) < x(indmax(end))
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
        rmin = fliplr(indmin(max(end-nbsym,1):end-1));
        rsym = indmin(end);
    else
        rmax = [lx,fliplr(indmax(max(end-nbsym+2,1):end))];
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
        rsym = lx;
    end
else
    if x(end) > x(indmin(end))
        rmax = fliplr(indmax(max(end-nbsym,1):end-1));
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
        rsym = indmax(end);
    else
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
        rmin = [lx,fliplr(indmin(max(end-nbsym+2,1):end))];
        rsym = lx;
    end
end

tlmin = 2*t(lsym)-t(lmin);
tlmax = 2*t(lsym)-t(lmax);
trmin = 2*t(rsym)-t(rmin);
trmax = 2*t(rsym)-t(rmax);

% in case symmetrized parts do not extend enough
if tlmin(1) > t(1) | tlmax(1) > t(1)
    if lsym == indmax(1)
        lmax = fliplr(indmax(1:min(end,nbsym)));
    else
        lmin = fliplr(indmin(1:min(end,nbsym)));
    end
    if lsym == 1
        error('bug')
    end
    lsym = 1;
    tlmin = 2*t(lsym)-t(lmin);
    tlmax = 2*t(lsym)-t(lmax);
end

if trmin(end) < t(lx) | trmax(end) < t(lx)
    if rsym == indmax(end)
        rmax = fliplr(indmax(max(end-nbsym+1,1):end));
    else
        rmin = fliplr(indmin(max(end-nbsym+1,1):end));
    end
    if rsym == lx
        error('bug')
    end
    rsym = lx;
    trmin = 2*t(rsym)-t(rmin);
    trmax = 2*t(rsym)-t(rmax);
end

xlmax =x(lmax);
xlmin =x(lmin);
xrmax =x(rmax);
xrmin =x(rmin);

tmin = [tlmin t(indmin) trmin];
tmax = [tlmax t(indmax) trmax];
xmin = [xlmin x(indmin) xrmin];
xmax = [xlmax x(indmax) xrmax];

%---------------------------------------------------------------------------------------------------

function [indmin, indmax, indzer] = extr(x)
%extracts the indices corresponding to extrema

% if(nargin==1)
%     t=1:length(x);
% end
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
    if ~isempty(debs)
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
    
    if ~isempty(imax)
        indmax = sort([indmax imax]);
    end
    
    if ~isempty(imin)
        indmin = sort([indmin imin]);
    end
    
end
function cmn = cum_mx( imn )
cmn=imn;
for i=1:size(imn,2)
    cmn(:,i)=sum(imn(:,i:end),2);
end