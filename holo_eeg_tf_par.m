function [] = holo_eeg_tf_par(S)
% Compute instantaneous power and phase in peri-stimulus time and frequency
% FORMAT [Dtf, Dtph] = eemd_eeg_tf(S)
%
% S                     - input structure
% fields of S:
%   S.D                 - MEEG object or filename of M/EEG mat-file with
%   
%   S.channels          - cell array of channel names. Can include generic
%                         wildcards: 'All', 'EEG', 'MEG' etc.
%
%   S.frequencies      - vector of frequencies of interest
%
%   S.timewin          - time window of interest in PST in ms. 
%   S.phase            - also save phase dataset (1) or not (0)

% Output:
% Dtf                   - M/EEG object with power (also written on disk)
% Dtph                  - M/EEG object with phase (also written on disk)
%__________________________________________________________________________
%S.ph_imf

if nargin == 0
    S = [];
end
if ~isfield(S, 'parallel')
    S.parallel = 0;
end
if ~isfield(S, 'collapse')
    S.collapse = 0;
end
if ~isfield(S, 'smooth')
    S.smooth = 1;
end

if ~isfield(S, 'CoreN')
    S.CoreN = feature('numcores') - 1;
end

%--------------------------------------------------------------------------
try
    Dfile = S.D;
catch
    [Dfile, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, Dfile = []; return; end
    S.D = Dfile;
end
D=spm_eeg_load(Dfile);
nd=length(size(D));
[pathstr, name, ext] = fileparts(Dfile); sep=filesep;
if nd==3   
    Dfile_fm=strcat(pathstr,sep,'timfm_',name,ext);
    Dfile_fm2=strcat(pathstr,sep,'timfm2_',name,ext);
    Dfile_am2=strcat(pathstr,sep,'timam2_',name,ext);
    dname=D.fname;
elseif nd==4
    dname=D.fname;
    dname=dname(6:end);
    Dfile_fm=strcat(pathstr,sep,'timfm_',dname);
    Dfile_fm2=strcat(pathstr,sep,'timfm2_',dname);
    Dfile_am2=strcat(pathstr,sep,'timam2_',dname);
else
   error('unknown filetype'); 
end
%D=spm_eeg_load(Dfile);
try
    D_fm=spm_eeg_load(Dfile_fm);
    D_fm2=spm_eeg_load(Dfile_fm2);
    D_am2=spm_eeg_load(Dfile_am2);
catch
    disp('2 layer EEMD files missing'); 
    return;
end
try
    if D_fm2.sdownLV>0
        sdownLV=round(D_fm2.sdownLV);
    else
        sdownLV=0;
    end
catch
    sdownLV=0;
end
if S.dyadic==1
   dyadic=1;
else
    dyadic=0;
end

%-Configure the analysis
%--------------------------------------------------------------------------
%samplerate=D.fsample;
timesol=S.settings.timesol;
freqsol=S.settings.freqsol;
FREQsol=S.settings.FREQsol;
fw0= S.settings.fw0;
fw1=S.settings.fw1;
Fw0= S.settings.Fw0;
Fw1=S.settings.Fw1;

if S.collapse==1
    if dyadic==0 && Fw0==0
       collapse=0;
    else
        collapse=1;
    end
else
    collapse=0;
end
if S.smooth==0
    smooth=0;
else
    smooth=1;
end
fscale=linspace(fw0,fw1,freqsol)';
if collapse
    Fscale=linspace(Fw0,Fw1,FREQsol)';
    if dyadic
        Fscale=[-Inf;Fscale];
    else
        Fscale=[0;Fscale];
    end
else
    Fscale=linspace(Fw0,Fw1,FREQsol)';
end
Nfrequencies=freqsol;
if collapse
    NFREQuencies=FREQsol+1;
else
    NFREQuencies=FREQsol;
end
%S.srcsImfs=0;
if ~isfield(S, 'srcImfs')
    S.srcsImfs=0;
end

if ~isfield(S, 'channels')
    S.channels = 'All';
end

chanind = D.selectchannels(S.channels);

if isempty(chanind)
    error('No channels selected.');
end
S.frequencies = fscale;
% use second layer fm to confine the time window of the HHS
if ~isfield(S, 'timewin')
    S.timewin = 1e3*[D_fm2.time(1) D_fm2.time(end)];
end
if isinf(S.timewin(1))||(S.timewin(1)<1e3*D_fm2.time(1))
    S.timewin(1) = 1e3*D_fm2.time(1);
end
if isinf(S.timewin(2))||(S.timewin(2)>1e3*D_fm2.time(end))
    S.timewin(2) = 1e3*D_fm2.time(end);
end
% t0=S.timewin(1)/1000;
% t1=S.timewin(2)/1000;
%tfsamplerate=timesol/(t1-t0);
timeind = D_fm2.indsample(1e-3*min(S.timewin)):D_fm2.indsample(1e-3*max(S.timewin));
if sdownLV>0
    tind_all=1:sdownLV:nsamples(D_fm);
    timeindf=tind_all(timeind);
else
    timeindf=timeind;
end
Nchannels = length(chanind);        
% Nsamples = length(timeind);
%dt=(t1-t0)/(Nsamples-1);
%nPT=Nsamples-1;


%%%For holo spectra

if S.fFt
      %%% prepare for averaging
    cl   = D.condlist;
    ni = zeros(1,D.nconditions);
    for i = 1:D.nconditions
        %w = pickconditions(D, deblank(cl{i}), 1)';
        w=indtrial(D,deblank(cl{i}),'GOOD');
        ni(i) = length(w);
        if ni(i) == 0
            warning('%s: No trials for trial type %d', D.fname, cl{i});
            continue;
        end
    end
    %%%%
    %Dmse = clone(D, ['mse_', 'm', int2str(S.m), '_', 'r', str_r, '_', 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)), D.fnamedat], [Nchannels Nscales D.ntrials]);
    if S.AllImfs
%         if collapse
%             DfFt = clone(D, ['mtfhh3d_', 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)), dname], [Nchannels (NFREQuencies+1)*Nfrequencies,timesol ,D.nconditions]);
%             DfFt = DfFt.frequencies(:, 1:(NFREQuencies+1)*Nfrequencies);
%         else           
            DfFt = clone(D, ['mtfhh3d_', 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)), dname], [Nchannels NFREQuencies*Nfrequencies,timesol ,D.nconditions]);
            DfFt = DfFt.frequencies(:, 1:NFREQuencies*Nfrequencies);
%         end
        %DfFt = timeonset(DfFt, D_fm.timeonset);
        %DfFt = timeonset(DfFt, 1e-3*min(S.timewin));
        DfFt = fsample(DfFt, timesol/(1e-3*max(S.timewin)-1e-3*min(S.timewin)));
        DfFt = timeonset(DfFt, 1e-3*(min(S.timewin)+0.5*(max(S.timewin)-min(S.timewin))/timesol));
        DfFt = transformtype(DfFt, 'TF');
        DfFt = chanlabels(DfFt, 1:Nchannels, D.chanlabels(chanind));
        DfFt = badchannels(DfFt, 1:Nchannels, D.badchannels(chanind));
        DfFt = chantype(DfFt, 1:Nchannels, D.chantype(chanind));
        DfFt = coor2D(DfFt, 1:Nchannels, coor2D(D,chanind));
    end
    if S.PartImfs
        holoImfs=S.holoImfs;
        hi1=int2str(holoImfs(1)); hi2=int2str(holoImfs(end));
%         if collapse
%             DfFtp = clone(D, ['mtfhh3d',hi1,hi2, 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_', dname], [Nchannels (NFREQuencies+1)*Nfrequencies,timesol ,D.nconditions]);
%             DfFtp = DfFtp.frequencies(:, 1:(NFREQuencies+1)*Nfrequencies);
%         else
            DfFtp = clone(D, ['mtfhh3d',hi1,hi2, 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_', dname], [Nchannels NFREQuencies*Nfrequencies,timesol ,D.nconditions]);
            DfFtp = DfFtp.frequencies(:, 1:NFREQuencies*Nfrequencies);
%         end
        %DfFtp = timeonset(DfFtp, 1e-3*min(S.timewin));
        DfFtp = fsample(DfFtp, timesol/(1e-3*max(S.timewin)-1e-3*min(S.timewin)));
        DfFtp = timeonset(DfFtp, 1e-3*(min(S.timewin)+0.5*(max(S.timewin)-min(S.timewin))/timesol));
        DfFtp = transformtype(DfFtp, 'TF');
        DfFtp = chanlabels(DfFtp, 1:Nchannels, D.chanlabels(chanind));
        DfFtp = badchannels(DfFtp, 1:Nchannels, D.badchannels(chanind));
        DfFtp = chantype(DfFtp, 1:Nchannels, D.chantype(chanind));
        DfFtp = coor2D(DfFtp, 1:Nchannels, coor2D(D,chanind));
    end
end
if S.fF
    if S.AllImfs
%         if collapse
%             DfF = clone(D, ['tomegaf_', 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)), dname], [Nchannels NFREQuencies+1 Nfrequencies D.ntrials]);
%             DfF = DfF.frequencies(:, Fscale);
%         else
            DfF = clone(D, ['tomegaf_', 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)), dname], [Nchannels NFREQuencies Nfrequencies D.ntrials]);
            DfF = DfF.frequencies(:, Fscale);
%         end
        DfF = timeonset(DfF, fw0);
        DfF = fsample(DfF, (freqsol-1)/(fw1-fw0));
        DfF = transformtype(DfF, 'TF');
        DfF = chanlabels(DfF, 1:Nchannels, D.chanlabels(chanind));
        DfF = badchannels(DfF, 1:Nchannels, D.badchannels(chanind));
        DfF = chantype(DfF, 1:Nchannels, D.chantype(chanind));
        DfF = coor2D(DfF, 1:Nchannels, coor2D(D,chanind));
    end
    if S.PartImfs
        holoImfs=S.holoImfs;
        hi1=int2str(holoImfs(1)); hi2=int2str(holoImfs(end));
        if collapse
            DfFp = clone(D, ['tomegaf',hi1,hi2, 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_', dname], [Nchannels NFREQuencies Nfrequencies D.ntrials]);
            DfFp = DfFp.frequencies(:, Fscale);
        else
            DfFp = clone(D, ['tomegaf',hi1,hi2, 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_', dname], [Nchannels NFREQuencies Nfrequencies D.ntrials]);
            DfFp = DfFp.frequencies(:, Fscale);
        end
        DfFp = timeonset(DfFp, fw0);
        DfFp = fsample(DfFp, (freqsol-1)/(fw1-fw0));
        DfFp = transformtype(DfFp, 'TF');
        DfFp = chanlabels(DfFp, 1:Nchannels, D.chanlabels(chanind));
        DfFp = badchannels(DfFp, 1:Nchannels, D.badchannels(chanind));
        DfFp = chantype(DfFp, 1:Nchannels, D.chantype(chanind));
        DfFp = coor2D(DfFp, 1:Nchannels, coor2D(D,chanind));
    end
end
if S.tF
    if S.AllImfs
%         if collapse
%            DtF = clone(D, ['tfhh2_', 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)), dname], [Nchannels NFREQuencies+1 timesol D.ntrials]);
%            DtF = DtF.frequencies(:, Fscale);
%         else
           DtF = clone(D, ['tfhh2_', 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)), dname], [Nchannels NFREQuencies timesol D.ntrials]);
           DtF = DtF.frequencies(:, Fscale);
%         end
       %DtF = timeonset(DtF, 1e-3*min(S.timewin));
       DtF = fsample(DtF, timesol/(1e-3*max(S.timewin)-1e-3*min(S.timewin)));
       DtF = timeonset(DtF, 1e-3*(min(S.timewin)+0.5*(max(S.timewin)-min(S.timewin))/timesol));
       DtF = transformtype(DtF, 'TF');
       DtF = chanlabels(DtF, 1:Nchannels, D.chanlabels(chanind));
       DtF = badchannels(DtF, 1:Nchannels, D.badchannels(chanind));
       DtF = chantype(DtF, 1:Nchannels, D.chantype(chanind));
       DtF = coor2D(DtF, 1:Nchannels, coor2D(D,chanind));
    end
    if S.PartImfs
       holoImfs=S.holoImfs; 
       hi1=int2str(holoImfs(1)); hi2=int2str(holoImfs(end));
%        if collapse
%            DtFp = clone(D, ['tfhh2_',hi1,hi2, 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_', dname], [Nchannels NFREQuencies+1 timesol D.ntrials]);
%            DtFp = DtFp.frequencies(:, Fscale);
%        else
           DtFp = clone(D, ['tfhh2_',hi1,hi2, 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_', dname], [Nchannels NFREQuencies timesol D.ntrials]);
           DtFp = DtFp.frequencies(:, Fscale);
%        end
       %DtFp = timeonset(DtFp, 1e-3*min(S.timewin));
       DtFp = fsample(DtFp, timesol/(1e-3*max(S.timewin)-1e-3*min(S.timewin)));
       DtFp = timeonset(DtFp, 1e-3*(min(S.timewin)+0.5*(max(S.timewin)-min(S.timewin))/timesol));
       DtFp = transformtype(DtFp, 'TF');
       DtFp = chanlabels(DtFp, 1:Nchannels, D.chanlabels(chanind));
       DtFp = badchannels(DtFp, 1:Nchannels, D.badchannels(chanind));
       DtFp = chantype(DtFp, 1:Nchannels, D.chantype(chanind));
       DtFp = coor2D(DtFp, 1:Nchannels, coor2D(D,chanind));
    end
end
if S.ft
    if S.AllImfs
       Dft = clone(D, ['tfhh1_', 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)), dname], [Nchannels Nfrequencies timesol D.ntrials]);
       Dft = Dft.frequencies(:, fscale);
       %Dft = timeonset(Dft, 1e-3*min(S.timewin));
       Dft = fsample(Dft, timesol/(1e-3*max(S.timewin)-1e-3*min(S.timewin)));
       Dft = timeonset(Dft, 1e-3*(min(S.timewin)+0.5*(max(S.timewin)-min(S.timewin))/timesol));
       Dft = transformtype(Dft, 'TF');
       Dft = chanlabels(Dft, 1:Nchannels, D.chanlabels(chanind));
       Dft = badchannels(Dft, 1:Nchannels, D.badchannels(chanind));
       Dft = chantype(Dft, 1:Nchannels, D.chantype(chanind));
       Dft = coor2D(Dft, 1:Nchannels, coor2D(D,chanind));
    end
    if S.PartImfs
       holoImfs=S.holoImfs; 
       hi1=int2str(holoImfs(1)); hi2=int2str(holoImfs(end));
       Dftp = clone(D, ['tfhh1_',hi1,hi2, 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_', dname], [Nchannels Nfrequencies timesol D.ntrials]);
       Dftp = Dftp.frequencies(:, fscale);
       %Dftp = timeonset(Dftp, 1e-3*min(S.timewin));
       Dftp = fsample(Dftp, timesol/(1e-3*max(S.timewin)-1e-3*min(S.timewin)));
       Dftp = timeonset(Dftp, 1e-3*(min(S.timewin)+0.5*(max(S.timewin)-min(S.timewin))/timesol));
       Dftp = transformtype(Dftp, 'TF');
       Dftp = chanlabels(Dftp, 1:Nchannels, D.chanlabels(chanind));
       Dftp = badchannels(Dftp, 1:Nchannels, D.badchannels(chanind));
       Dftp = chantype(Dftp, 1:Nchannels, D.chantype(chanind));
       Dftp = coor2D(Dftp, 1:Nchannels, coor2D(D,chanind));
    end
end

%%%smooth the data
%q_pre=fspecial('gaussian', 3, 0.5);
q_pre2=fspecial('gaussian', 5, 0.5);
q_pre3=fspecial('gaussian', [1 7], 0.75);
%q_pre4=fspecial('gaussian', 5, 0.75);
%q2=fspecial('average',[1,3]);
vec_nimf=D_fm2.nimf;
vec_nIMF=D_fm2.nIMF;
max_imf=D_fm2.maximf;
max_IMF=D_fm2.maxIMF;

%-Prepare parallel comp
if S.parallel
    p=gcp('nocreate');
    if ~isempty(p)
        delete(p);
    end
    parpool(S.CoreN);
end

%-Run the analysis on all trials
%--------------------------------------------------------------------------
cl   = D.condlist;
for i = 1:D.nconditions
%     w = pickconditions(D, deblank(cl{i}), 0);
%     w_nb = pickconditions(D, deblank(cl{i}), 1);
%     btrials=setdiff(w,w_nb);
    w=indtrial(D,deblank(cl{i}));
    btrials=indtrial(D,deblank(cl{i}),'BAD');
    if S.fFt
        if S.AllImfs
%             if collapse
%                 fFt_trialdata=zeros(Nchannels,Nfrequencies,NFREQuencies+1,timesol);
%             else
                fFt_trialdata=zeros(Nchannels,Nfrequencies,NFREQuencies,timesol);
%             end
        end
        if S.PartImfs
%             if collapse
%                 fFtp_trialdata=zeros(Nchannels,Nfrequencies,NFREQuencies+1,timesol);
%             else
                fFtp_trialdata=zeros(Nchannels,Nfrequencies,NFREQuencies,timesol);
%             end
        end
    end
        
    for k = w
        fm_trialdata=D_fm(:,:,timeindf,k);
        fm_trialdata=permute(fm_trialdata,[3,2,1]);
        pfm_trialdata=mat2cell(fm_trialdata,size(fm_trialdata,1),size(fm_trialdata,2),ones(1,size(fm_trialdata,3)));
        fm2_trialdata=D_fm2(:,:,timeind,k);
        fm2_trialdata=reshape(fm2_trialdata,size(fm2_trialdata,1),max_IMF,max_imf,[]);
        fm2_trialdata=permute(fm2_trialdata,[4,2,3,1]);
        pfm2_trialdata=mat2cell(fm2_trialdata,size(fm2_trialdata,1),size(fm2_trialdata,2),size(fm2_trialdata,3),ones(1,size(fm2_trialdata,4)));
        am2_trialdata=D_am2(:,:,timeind,k);
        am2_trialdata=reshape(am2_trialdata,size(am2_trialdata,1),max_IMF,max_imf,[]);
        am2_trialdata=permute(am2_trialdata,[4,2,3,1]);
        pam2_trialdata=mat2cell(am2_trialdata,size(am2_trialdata,1),size(am2_trialdata,2),size(am2_trialdata,3),ones(1,size(am2_trialdata,4)));
%         am2_trialdata=D_am2(:,:,timeind,k); 
%         am2_trialdata_r=reshape(am2_trialdata,size(am2_trialdata,1),max_IMF,max_imf,[]); 
        if ~ S.srcsImfs
            pcurr_nimf = vec_nimf(k,:)-1;
        else
            pcurr_nimf = vec_nimf(k,:);
        end
        pcurr_nIMF = vec_nIMF(k,:);
        if S.fFt
              if ~ismember(k,btrials)
                if S.AllImfs
                    parfor ch=1:Nchannels  
                           curr_nimf = pcurr_nimf(ch);
                           curr_nIMF = pcurr_nIMF(ch);
                           if curr_nimf>=1 && curr_nIMF>=1
                               [All_nt]=nspplotf3d_tres3x(pfm_trialdata{ch}(:,1:curr_nimf),pfm2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),pam2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),[],[],freqsol,FREQsol,fw0,fw1,Fw0, Fw1,timesol,S);
                                All_ntc=squeeze(sum(All_nt,4));
                                if smooth
                                nts = smooth3(All_ntc,'gaussian',3,0.65);
                                nts = smooth3(nts,'gaussian',3,0.65);
                                nts = smooth3(nts,'gaussian',[1,1,7],0.65);
                                else
                                    nts=All_ntc;
                                end
                                %nts = smooth3(nts_pre,'box',3);
                               fFt_trialdata(ch,:,:,:)=fFt_trialdata(ch,:,:,:)+shiftdim(nts,-1);  % remove all 0.5 power on 20160530
                           end
                    end

                end
                if S.PartImfs
                    parfor ch=1:Nchannels  
                           curr_nimf = pcurr_nimf(ch);
                           curr_nIMF = pcurr_nIMF(ch);
                           if curr_nimf>=1 && curr_nIMF>=1
                           [All_nt]=nspplotf3d_tres3x(pfm_trialdata{ch}(:,1:curr_nimf),pfm2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),pam2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),[],[],freqsol,FREQsol,fw0,fw1,Fw0, Fw1,timesol,S);
                            All_ntc=squeeze(sum(All_nt(:,:,:,holoImfs),4));
                            if smooth
                            nts = smooth3(All_ntc,'gaussian',3,0.65);
                            nts = smooth3(nts,'gaussian',3,0.65);
                            nts = smooth3(nts,'gaussian',[1,1,7],0.65);
                            else
                                 nts=All_ntc;
                            end
                            %nts = smooth3(nts_pre,'box',3);
                           fFtp_trialdata(ch,:,:,:)=fFtp_trialdata(ch,:,:,:)+shiftdim(nts,-1);  % remove all 0.5 power on 20160530
                           end
                    end
                end
             end
        end
        if S.fF
            if S.AllImfs
%                 if collapse
%                     fF_trialdata=zeros(Nchannels,NFREQuencies+1,Nfrequencies);
%                 else
                    fF_trialdata=zeros(Nchannels,NFREQuencies,Nfrequencies);
%                 end
               parfor ch=1:Nchannels  
                  curr_nimf = pcurr_nimf(ch);
                  curr_nIMF = pcurr_nIMF(ch);
                  if curr_nimf>=1 && curr_nIMF>=1
                  [All_nt]=nspplotf3d_tres3x(pfm_trialdata{ch}(:,1:curr_nimf),pfm2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),pam2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),[],[],freqsol,FREQsol,fw0,fw1,Fw0, Fw1,timesol,S);
                  All_ntc=squeeze(mean(All_nt,3)); % 20160420 change from sum to mean
                  All_ntc2=squeeze(sum(All_ntc,3));
%                   nts_pre=filter2(q_pre2,All_ntc2);
%                   nts=filter2(q_pre2,nts_pre);
%                   fF_trialdata(ch,:,:)=nts';
                  fF_trialdata(ch,:,:)=All_ntc2';
                  end
               end
               if collapse
                    for ch=1:Nchannels
                        am_ntc=shiftdim(fF_trialdata(ch,2:end,:));
                        if smooth
                        nts_pre=filter2(q_pre2,am_ntc);
                        nts=filter2(q_pre2,nts_pre);
                        else
                             nts=am_ntc;
                        end
                        %nts=filter2(q2, nts);
                        fF_trialdata(ch,2:end,:)=shiftdim(nts,-1);
                        am_dc=shiftdim(fF_trialdata(ch,1,:));
                        if smooth
                        ndc_pre=filter2(q_pre2,am_dc);
                        ndc=filter2(q_pre2,ndc_pre);
                        else
                            ndc=am_dc;
                        end
                        fF_trialdata(ch,1,:)=shiftdim(ndc,-2);
                    end
                else
                    for ch=1:Nchannels
                        am_ntc=shiftdim(fF_trialdata(ch,:,:));
                        if smooth
                        nts_pre=filter2(q_pre2,am_ntc);
                        nts=filter2(q_pre2,nts_pre);
                        else
                            nts=am_ntc;
                        end
                        %nts=filter2(q2, nts);
                        fF_trialdata(ch,:,:)=shiftdim(nts,-1);
                    end
                end
               DfF(:,:,:,k)=fF_trialdata;
            end
            if S.PartImfs
%                if collapse
%                     fFp_trialdata=zeros(Nchannels,NFREQuencies+1,Nfrequencies);
%                 else
                    fFp_trialdata=zeros(Nchannels,NFREQuencies,Nfrequencies);
%                 end
               parfor ch=1:Nchannels  
                  curr_nimf = pcurr_nimf(ch);
                  curr_nIMF = pcurr_nIMF(ch);
                  if curr_nimf>=1 && curr_nIMF>=1
                  [All_nt]=nspplotf3d_tres3x(pfm_trialdata{ch}(:,1:curr_nimf),pfm2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),pam2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),[],[],freqsol,FREQsol,fw0,fw1,Fw0, Fw1,timesol,S);
                  All_ntc=squeeze(mean(All_nt,3)); % 20160420 change from sum to mean
                  All_ntc2=squeeze(sum(All_ntc(:,:,holoImfs),3));
%                   nts_pre=filter2(q_pre2,All_ntc2);
%                   nts=filter2(q_pre2,nts_pre);
                  fFp_trialdata(ch,:,:)=All_ntc2';
                  end
               end
               if collapse
                    for ch=1:Nchannels
                        am_ntc=shiftdim(fFp_trialdata(ch,2:end,:));
                        if smooth
                        nts_pre=filter2(q_pre2,am_ntc);
                        nts=filter2(q_pre2,nts_pre);
                        else
                            nts=am_ntc;
                        end
                        %nts=filter2(q2, nts);
                        fFp_trialdata(ch,2:end,:)=shiftdim(nts,-1);
                        am_dc=shiftdim(fFp_trialdata(ch,1,:));
                        if smooth
                        ndc_pre=filter2(q_pre2,am_dc);
                        ndc=filter2(q_pre2,ndc_pre);
                        else
                            ndc=am_dc;
                        end
                        fFp_trialdata(ch,1,:)=shiftdim(ndc,-2);
                    end
                else
                    for ch=1:Nchannels
                        am_ntc=shiftdim(fFp_trialdata(ch,:,:));
                        if smooth
                        nts_pre=filter2(q_pre2,am_ntc);
                        nts=filter2(q_pre2,nts_pre);
                        else
                            nts=am_ntc;
                        end
                        %nts=filter2(q2, nts);
                        fFp_trialdata(ch,:,:)=shiftdim(nts,-1);
                    end
                end
               DfFp(:,:,:,k)=fFp_trialdata;
            end
        end
        if S.tF
            if S.AllImfs
%                 if collapse
%                     tF_trialdata=zeros(Nchannels,NFREQuencies+1,timesol); 
%                 else
                    tF_trialdata=zeros(Nchannels,NFREQuencies,timesol); 
%                 end
               parfor ch=1:Nchannels  
                  curr_nimf = pcurr_nimf(ch);
                  curr_nIMF = pcurr_nIMF(ch);
                  if curr_nimf>=1 && curr_nIMF>=1
                  [All_nt]=nspplotf3d_tres3x(pfm_trialdata{ch}(:,1:curr_nimf),pfm2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),pam2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),[],[],freqsol,FREQsol,fw0,fw1,Fw0, Fw1,timesol,S);
                  [Omega_time_spec]=produce_Omega_time(All_nt);
%                   nts_pre=filter2(q_pre,Omega_time_spec);
%                   nts_pre=filter2(q_pre,nts_pre);
%                   nts=filter2(q2, nts_pre);
%                   tF_trialdata(ch,:,:)=nts;
                    tF_trialdata(ch,:,:)=Omega_time_spec; 
                  end
               end
               if collapse
                    for ch=1:Nchannels
                        am_ntc=shiftdim(tF_trialdata(ch,2:end,:));
                        if smooth
                        nts_pre=filter2(q_pre3,am_ntc);
                        nts_pre=filter2(q_pre2,nts_pre);
                        nts_pre=filter2(q_pre2,nts_pre);
                        nts=nts_pre;
                        else
                            nts=am_ntc;
                        end
                        tF_trialdata(ch,2:end,:)=shiftdim(nts,-1);
                        am_dc=shiftdim(tF_trialdata(ch,1,:));
                        if smooth
                        ndc_pre=filter2(q_pre3,am_dc);
                        ndc_pre=filter2(q_pre2,ndc_pre);
                        ndc_pre=filter2(q_pre2,ndc_pre);
                        ndc=ndc_pre;
                        else
                            ndc=am_dc;
                        end
                        tF_trialdata(ch,1,:)=shiftdim(ndc,-2);
                    end
                else
                    for ch=1:Nchannels
                        am_ntc=shiftdim(tF_trialdata(ch,:,:));
                        if smooth
                        nts_pre=filter2(q_pre3,am_ntc);
                        nts_pre=filter2(q_pre2,nts_pre);
                        nts_pre=filter2(q_pre2,nts_pre);
                        nts=nts_pre;
                        else
                            nts=am_ntc;
                        end
                        tF_trialdata(ch,:,:)=shiftdim(nts,-1);
                    end
                end
               DtF(:,:,:,k)=tF_trialdata;
            end
            if S.PartImfs
%                 if collapse
%                     tFp_trialdata=zeros(Nchannels,NFREQuencies+1,timesol); 
%                 else
                    tFp_trialdata=zeros(Nchannels,NFREQuencies,timesol); 
%                 end

               parfor ch=1:Nchannels  
                  curr_nimf = pcurr_nimf(ch);
                  curr_nIMF = pcurr_nIMF(ch);
                  if curr_nimf>=1 && curr_nIMF>=1
                  [All_nt]=nspplotf3d_tres3x(pfm_trialdata{ch}(:,1:curr_nimf),pfm2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),pam2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),[],[],freqsol,FREQsol,fw0,fw1,Fw0, Fw1,timesol,S);
                  [Omega_time_spec]=produce_Omega_time(All_nt(:,:,:,holoImfs));
%                   nts_pre=filter2(q_pre,Omega_time_spec);
%                   nts_pre=filter2(q_pre,nts_pre);
%                   nts=filter2(q2, nts_pre);
                  tFp_trialdata(ch,:,:)=Omega_time_spec;
                  end
               end
               if collapse
                    for ch=1:Nchannels
                        am_ntc=shiftdim(tFp_trialdata(ch,2:end,:));
                        if smooth
                        nts_pre=filter2(q_pre3,am_ntc);
                        nts_pre=filter2(q_pre2,nts_pre);
                        nts_pre=filter2(q_pre2,nts_pre);
                        nts=nts_pre;
                        else
                            nts=am_ntc;
                        end
                        tFp_trialdata(ch,2:end,:)=shiftdim(nts,-1);
                        am_dc=shiftdim(tFp_trialdata(ch,1,:));
                        if smooth
                        ndc_pre=filter2(q_pre3,am_dc);
                        ndc_pre=filter2(q_pre2,ndc_pre);
                        ndc_pre=filter2(q_pre2,ndc_pre);
                        ndc=ndc_pre;
                        else
                            ndc=am_dc;
                        end
                        tFp_trialdata(ch,1,:)=shiftdim(ndc,-2);
                    end
                else
                    for ch=1:Nchannels
                        am_ntc=shiftdim(tFp_trialdata(ch,:,:));
                        if smooth
                        nts_pre=filter2(q_pre3,am_ntc);
                        nts_pre=filter2(q_pre2,nts_pre);
                        nts_pre=filter2(q_pre2,nts_pre);
                        nts=nts_pre;
                        else
                            nts=am_ntc;
                        end
                        tFp_trialdata(ch,:,:)=shiftdim(nts,-1);
                    end
                end
               DtFp(:,:,:,k)=tFp_trialdata;
            end
        end
        if S.ft
            if S.AllImfs
               ft_trialdata=zeros(Nchannels,Nfrequencies,timesol); 
               parfor ch=1:Nchannels  
                  curr_nimf = pcurr_nimf(ch);
                  curr_nIMF = pcurr_nIMF(ch);
                  if curr_nimf>=1 && curr_nIMF>=1
                  [All_nt]=nspplotf3d_tres3x(pfm_trialdata{ch}(:,1:curr_nimf),pfm2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),pam2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),[],[],freqsol,FREQsol,fw0,fw1,Fw0, Fw1,timesol,S);
                  [freq_time_spec]=produce_freq_time(All_nt);
                  if smooth
                  nts_pre=filter2(q_pre3,freq_time_spec);
                  nts_pre=filter2(q_pre2,nts_pre);
                  nts_pre=filter2(q_pre2,nts_pre);
                  nts=nts_pre;
                  else
                      nts=freq_time_spec;
                  end
                  ft_trialdata(ch,:,:)=nts;
                  end
               end
               Dft(:,:,:,k)=ft_trialdata;
            end
            if S.PartImfs
               ftp_trialdata=zeros(Nchannels,Nfrequencies,timesol); 
               parfor ch=1:Nchannels  
                  curr_nimf = pcurr_nimf(ch);
                  curr_nIMF = pcurr_nIMF(ch);
                  if curr_nimf>=1 && curr_nIMF>=1
                  [All_nt]=nspplotf3d_tres3x(pfm_trialdata{ch}(:,1:curr_nimf),pfm2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),pam2_trialdata{ch}(:,1:curr_nIMF,1:curr_nimf),[],[],freqsol,FREQsol,fw0,fw1,Fw0, Fw1,timesol,S);
                  [freq_time_spec]=produce_freq_time(All_nt(:,:,:,holoImfs));
                  if smooth
                  nts_pre=filter2(q_pre3,freq_time_spec);
                  nts_pre=filter2(q_pre2,nts_pre);
                  nts_pre=filter2(q_pre2,nts_pre);
                  nts=nts_pre;
                  else
                      nts=freq_time_spec;
                  end
                  ftp_trialdata(ch,:,:)=nts;
                  end
               end
               Dftp(:,:,:,k)=ftp_trialdata;
            end
        end
    end
    if S.fFt
        if (length(w)-length(btrials))>=1
            if S.AllImfs
                if collapse
                    for ch=1:Nchannels
                        am_ntc=shiftdim(fFt_trialdata(ch,:,2:end,:));
                        nts_pre = smooth3(am_ntc,'gaussian',3,0.65);
                        nts = smooth3(nts_pre,'gaussian',[3,3,1],0.65);
                        fFt_trialdata(ch,:,2:end,:)=shiftdim(nts,-1);
                        am_dc=shiftdim(fFt_trialdata(ch,:,1,:));
                        ndc_pre = smooth3(am_dc,'gaussian',3,0.65);
                        ndc = smooth3(ndc_pre,'gaussian',[3,1,1],0.65);
                        fFt_trialdata(ch,:,1,:)=shiftdim(ndc,-1);
                    end
                else
                    for ch=1:Nchannels
                        am_ntc=shiftdim(fFt_trialdata(ch,:,:,:));
                        nts_pre = smooth3(am_ntc,'gaussian',3,0.65);
                        nts = smooth3(nts_pre,'gaussian',[3,3,1],0.65);
                        fFt_trialdata(ch,:,:,:)=shiftdim(nts,-1);
                    end
                end
                    DfFt(:,:,:,i)=reshape(fFt_trialdata,Nchannels,Nfrequencies*NFREQuencies,timesol)/(length(w)-length(btrials));
                
            end
            if S.PartImfs
                if collapse
                    for ch=1:Nchannels
                        am_ntc=shiftdim(fFtp_trialdata(ch,:,2:end,:));
                        nts_pre = smooth3(am_ntc,'gaussian',3,0.65);
                        nts = smooth3(nts_pre,'gaussian',[3,3,1],0.65);
                        fFtp_trialdata(ch,:,2:end,:)=shiftdim(nts,-1);
                        am_dc=shiftdim(fFtp_trialdata(ch,:,1,:));
                        ndc_pre = smooth3(am_dc,'gaussian',3,0.65);
                        ndc = smooth3(ndc_pre,'gaussian',[3,1,1],0.65);
                        fFtp_trialdata(ch,:,1,:)=shiftdim(ndc,-1);
                    end
                else
                    for ch=1:Nchannels
                        am_ntc=shiftdim(fFtp_trialdata(ch,:,:,:));
                        nts_pre = smooth3(am_ntc,'gaussian',3,0.65);
                        nts = smooth3(nts_pre,'gaussian',[3,3,1],0.65);
                        fFtp_trialdata(ch,:,:,:)=shiftdim(nts,-1);
                    end
                end
                    DfFtp(:,:,:,i)=reshape(fFtp_trialdata,Nchannels,Nfrequencies*NFREQuencies,timesol)/(length(w)-length(btrials));
            end
        end
    end
end

%-Save new M/EEG dataset(s)
if S.fFt
    if S.AllImfs
        DfFt.Nfrequencies=Nfrequencies;
        DfFt.NFREQuencies=NFREQuencies;
        DfFt = type(DfFt, 'evoked');
        DfFt = conditions(DfFt, '', cl);
        DfFt = repl(DfFt, '', ni);
        DfFt.dyadic=dyadic;
        DfFt.collapse=collapse;
        DfFt.fscale=fscale;
        DfFt.Fscale=Fscale;
        DfFt = DfFt.history(mfilename, S);
        save(DfFt);
    end
    if S.PartImfs
        DfFtp.Nfrequencies=Nfrequencies;
        DfFtp.NFREQuencies=NFREQuencies;
        DfFtp = type(DfFtp, 'evoked');
        DfFtp = conditions(DfFtp, '', cl);
        DfFtp = repl(DfFtp, '', ni);
        DfFtp.dyadic=dyadic;
        DfFtp.collapse=collapse;
        DfFtp.fscale=fscale;
        DfFtp.Fscale=Fscale;
        DfFtp = DfFtp.history(mfilename, S);
        save(DfFtp);
    end
end
if S.fF
    if S.AllImfs
       DfF.dyadic=dyadic;
       DfF.collapse=collapse;
       DfF = DfF.history(mfilename, S);
       save(DfF);
    end
    if S.PartImfs
        DfFp.dyadic=dyadic; 
        DfFp.collapse=collapse;
        DfFp = DfFp.history(mfilename, S);
       save(DfFp);
    end
end
if S.tF
    if S.AllImfs
        DtF.dyadic=dyadic; 
        DtF.collapse=collapse;
        DtF = DtF.history(mfilename, S);
       save(DtF);
    end
    if S.PartImfs
        DtFp.dyadic=dyadic; 
        DtFp.collapse=collapse;
        DtFp = DtFp.history(mfilename, S);
       save(DtFp);
    end
end
if S.ft
    if S.AllImfs
        Dft.dyadic=dyadic; 
        %Dft.collapse=collapse;
        Dft = Dft.history(mfilename, S);
       save(Dft);
    end
    if S.PartImfs
        Dftp.dyadic=dyadic;
        %Dftp.collapse=collapse;
        Dftp = Dftp.history(mfilename, S);
       save(Dftp);
    end
end
