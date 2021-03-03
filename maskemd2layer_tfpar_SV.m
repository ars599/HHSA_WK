function [Dimf, DIMF, Dfm, DFM, DAM] = maskemd2layer_tfpar_SV(S)
% This code implements the SPM wrapper of the 2-layer EMD with the enhanced algorithm of EMD(masking EMD)
% History: This code was written by Wei-Kuang Liang
% FORMAT [Dimf, DIMF, Dfm, DFM, DAM] = maskemd2layer_tfpar_SV(S)
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
%   S.eemd
% Output:
% Dimf                   - M/EEG object with first-layer IMFs (also written on disk)
% Dfm                  - M/EEG object with instantaneous frequency (also written on disk)
% DIMF                   - M/EEG object with second-layer IMFs (also written on disk)
% DFM                  - M/EEG object with instantaneous AM frequency (also written on disk)
% DAM                   - M/EEG object with amplitudes of second-layer IMFs (also written on disk)


if nargin == 0
    S = [];
end


%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end
if ~ isa(D,'meeg')
    D = spm_eeg_load(D);
end
if isequal(D.type, 'continuous')
    error('Time-frequency analysis can only be applied to epoched data');
end

%-Configure the analysis
%--------------------------------------------------------------------------
samplerate=D.fsample;
if ~isfield(S, 'channels')
    S.channels = 'All';
end

if ~isfield(S, 'from1stEmd')
    S.from1stEmd = 0;
end
if ~isfield(S, 'sdownLV')
    S.sdownLV = 0;
end
if ~isfield(S, 'parallel')
    S.parallel = 0;
end
if ~isfield(S, 'CoreN')
    S.CoreN = feature('numcores') - 1;
end

sdownLV=round(S.sdownLV); % downsample when store data for 2nd Layer IMFs
if sdownLV==1
   sdownLV=0; 
end

chanind = D.selectchannels(S.channels);

if isempty(chanind)
    error('No channels selected.');
end

if ~isfield(S, 'timewin')
    S.timewin = 1e3*[D.time(1) D.time(end)];
end
if isinf(S.timewin(1))|| (S.timewin(1)<1e3*D.time(1))
    S.timewin(1) = 1e3*D.time(1);
end
if isinf(S.timewin(2))|| (S.timewin(2)>1e3*D.time(end))
    S.timewin(2) = 1e3*D.time(end);
end
    
t0=S.timewin(1)/1000;
t1=S.timewin(2)/1000;
timeind = D.indsample(1e-3*min(S.timewin)):D.indsample(1e-3*max(S.timewin));

if ~isfield(S, 'phase')
    S.phase = 0;
end
if ~isfield(S, 'phase2')
    S.phase2 = 0;
end
if ~isfield(S, 'NEnsemble')
    S.NEnsemble = 0;
end
if ~isfield(S, 'ENoise')
    S.ENoise = 2;
end
if ~isfield(S, 'NEnsemble2')
    S.NEnsemble2 = 0;
end
if ~isfield(S, 'ENoise2')
    S.ENoise2 = 0.25;
end
if ~isfield(S, 'shiftLevel')
    S.shiftLevel = 0;
end
if ~isfield(S, 'ifmethod')
    S.ifmethod = 'zc';
end
if ~isfield(S, 'ifmethod2')
    S.ifmethod2 = 'zc';
end
if ~isfield(S, 'nyquist')
    S.nyquist = 0;
end

maximf=fix(log2(length(timeind)));
maxIMF=maximf;
Nchannels = length(chanind);        
Nsamples = length(timeind);
if sdownLV>0
   dNsamples=length(timeind(1:sdownLV:end));
   dsamplerate=samplerate/sdownLV;
else
   dNsamples=Nsamples; 
   dsamplerate=samplerate;
end

dt=(t1-t0)/(Nsamples-1);
Dimf = clone(D, ['timf_' D.fname], [Nchannels maximf Nsamples D.ntrials]);
Dimf = Dimf.frequencies(:, 1:maximf);
Dimf = timeonset(Dimf, t0);
Dimf = fsample(Dimf, samplerate);
Dimf = transformtype(Dimf, 'TF');
Dimf = chanlabels(Dimf, 1:Nchannels, D.chanlabels(chanind));
Dimf = badchannels(Dimf, 1:Nchannels, D.badchannels(chanind));
Dimf = chantype(Dimf, 1:Nchannels, D.chantype(chanind));
Dimf = coor2D(Dimf, 1:Nchannels, coor2D(D,chanind));
% fm & am
Dfm = clone(Dimf, ['timfm_' D.fname]);
Dam = clone(Dimf, ['timam_' D.fname]);
%imf phase
if S.phase 
    Dimph = clone(Dimf, ['timph_' D.fname]);
    Dimph = transformtype(Dimph, 'TFphase');
end
%IMF 2
DIMF = clone(Dimf, ['timf2_' D.fname], [Nchannels maxIMF*maximf dNsamples D.ntrials]);
DIMF = DIMF.frequencies(:, 1:(maxIMF*maximf));
DIMF = fsample(DIMF, dsamplerate);
DIMF.sdownLV=sdownLV;
% FM & AM
DFM = clone(DIMF, ['timfm2_' D.fname]);
DAM = clone(DIMF, ['timam2_' D.fname]);
if S.phase2
    DIMPH = clone(DIMF, ['timph2_' D.fname]);
    DIMPH = transformtype(DIMPH, 'TFphase');
end
if S.parallel
    try
        p=gcp('nocreate');
        if ~isempty(p)
            delete(p);
        end
        parpool(S.CoreN);
    catch
    end
end
phase=S.phase;
phase2=S.phase2;

vec_nimf=zeros(D.ntrials,Nchannels);
vec_nIMF2=zeros(D.ntrials,Nchannels);
mx_nIMF2=zeros(D.ntrials,Nchannels,maximf);
subD=D(chanind, timeind, :);
%-Run the analysis on all trials
%--------------------------------------------------------------------------
for k = 1:D.ntrials
        trial_data=subD(:,:,k);
        trial_data=squeeze(trial_data);
        trial_imfdata=zeros(Nchannels,maximf, Nsamples);
        trial_imfmdata=zeros(Nchannels,maximf, Nsamples);
        trial_imamdata=zeros(Nchannels,maximf, Nsamples);
        trial_IMFdata=zeros(Nchannels,maxIMF,maximf, dNsamples);
        trial_IMFMdata=zeros(Nchannels,maxIMF,maximf, dNsamples);
        trial_IMAMdata=zeros(Nchannels,maxIMF,maximf, dNsamples);
        ptrial_imfdata = cell(Nchannels,1);
        ptrial_imfmdata = cell(Nchannels,1);
        ptrial_imamdata = cell(Nchannels,1);
        ptrial_IMFdata= cell(Nchannels,1);
        ptrial_IMFMdata= cell(Nchannels,1);
        ptrial_IMAMdata= cell(Nchannels,1);
        if phase
           trial_imfphdata=zeros(Nchannels,maximf, Nsamples);
           ptrial_imfphdata = cell(Nchannels,1);
        end
        if phase2
           trial_IMFPHdata=zeros(Nchannels,maxIMF,maximf, dNsamples);
           ptrial_IMFPHdata = cell(Nchannels,1);
        end
        SD=std(reshape(trial_data,1,[]),1);
        tmpmx_nIMF2=cell(Nchannels,1);
        parfor ch=1:Nchannels  
            tmpdata=trial_data(ch,:);
            sdc=std(tmpdata,1);
            if sdc< (1e-10*SD)
               vec_nimf(k,ch)=0;
               vec_nIMF2(k,ch)=0; 
               continue;
            end
            if ~phase
                if ~phase2
                   [ptrial_imfmdata{ch},ptrial_imamdata{ch},ptrial_IMFMdata{ch},ptrial_IMAMdata{ch},ptrial_imfdata{ch},ptrial_IMFdata{ch},tmpmx_nIMF2{ch}] = multi_EMD_DCM_SV(tmpdata,samplerate,[],[],S);
                else
                   [ptrial_imfmdata{ch},ptrial_imamdata{ch},ptrial_IMFMdata{ch},ptrial_IMAMdata{ch},ptrial_imfdata{ch},ptrial_IMFdata{ch},tmpmx_nIMF2{ch},~,ptrial_IMFPHdata{ch}] = multi_EMD_DCM_SV(tmpdata,samplerate,[],[],S);
                end
            else
                if ~phase2
                   [ptrial_imfmdata{ch},ptrial_imamdata{ch},ptrial_IMFMdata{ch},ptrial_IMAMdata{ch},ptrial_imfdata{ch},ptrial_IMFdata{ch},tmpmx_nIMF2{ch},ptrial_imfphdata{ch}] = multi_EMD_DCM_SV(tmpdata,samplerate,[],[],S);
                else
                   [ptrial_imfmdata{ch},ptrial_imamdata{ch},ptrial_IMFMdata{ch},ptrial_IMAMdata{ch},ptrial_imfdata{ch},ptrial_IMFdata{ch},tmpmx_nIMF2{ch},ptrial_imfphdata{ch},ptrial_IMFPHdata{ch}] = multi_EMD_DCM_SV(tmpdata,samplerate,[],[],S);
                end
            end
            curImf=size(ptrial_imfdata{ch},2);
            vec_nimf(k,ch)=curImf;
            curIMF2=size(ptrial_IMFdata{ch},2);
            vec_nIMF2(k,ch)=curIMF2;
        end
        if sdownLV==0
            for kk = 1:Nchannels
                trial_imfdata(kk,1:size(ptrial_imfdata{kk},2),1:size(ptrial_imfdata{kk},1)) = ptrial_imfdata{kk}';
                trial_imfmdata(kk,1:size(ptrial_imfmdata{kk},2),1:size(ptrial_imfmdata{kk},1)) = ptrial_imfmdata{kk}';
                trial_imamdata(kk,1:size(ptrial_imamdata{kk},2),1:size(ptrial_imamdata{kk},1)) = ptrial_imamdata{kk}';
                trial_IMFdata(kk,1:size(ptrial_IMFdata{kk},2),1:size(ptrial_IMFdata{kk},3),1:size(ptrial_IMFdata{kk},1)) = permute(ptrial_IMFdata{kk},[2,3,1]);
                trial_IMFMdata(kk,1:size(ptrial_IMFMdata{kk},2),1:size(ptrial_IMFMdata{kk},3),1:size(ptrial_IMFMdata{kk},1)) = permute(ptrial_IMFMdata{kk},[2,3,1]);
                trial_IMAMdata(kk,1:size(ptrial_IMAMdata{kk},2),1:size(ptrial_IMAMdata{kk},3),1:size(ptrial_IMAMdata{kk},1)) = permute(ptrial_IMAMdata{kk},[2,3,1]);
                tmpn2=length(tmpmx_nIMF2{kk});
                mx_nIMF2(k,kk,1:tmpn2)=tmpmx_nIMF2{kk};
                if phase
                    trial_imfphdata(kk,1:size(ptrial_imfphdata{kk},2),1:size(ptrial_imfphdata{kk},1)) = ptrial_imfphdata{kk}';
                end
                if phase2
                    trial_IMFPHdata(kk,1:size(ptrial_IMFPHdata{kk},2),1:size(ptrial_IMFPHdata{kk},3),1:size(ptrial_IMFPHdata{kk},1)) = permute(ptrial_IMFPHdata{kk},[2,3,1]);
                end
            end
        else
            for kk = 1:Nchannels
                trial_imfdata(kk,1:size(ptrial_imfdata{kk},2),1:size(ptrial_imfdata{kk},1)) = ptrial_imfdata{kk}';
                trial_imfmdata(kk,1:size(ptrial_imfmdata{kk},2),1:size(ptrial_imfmdata{kk},1)) = ptrial_imfmdata{kk}';
                trial_imamdata(kk,1:size(ptrial_imamdata{kk},2),1:size(ptrial_imamdata{kk},1)) = ptrial_imamdata{kk}';
                trial_IMFdata(kk,1:size(ptrial_IMFdata{kk},2),1:size(ptrial_IMFdata{kk},3),1:dNsamples) = permute(spmm_downsample(ptrial_IMFdata{kk},sdownLV),[2,3,1]);
                trial_IMFMdata(kk,1:size(ptrial_IMFMdata{kk},2),1:size(ptrial_IMFMdata{kk},3),1:dNsamples) = permute(spmm_downsample(ptrial_IMFMdata{kk},sdownLV),[2,3,1]);
                trial_IMAMdata(kk,1:size(ptrial_IMAMdata{kk},2),1:size(ptrial_IMAMdata{kk},3),1:dNsamples) = permute(spmm_downsample(ptrial_IMAMdata{kk},sdownLV),[2,3,1]);
                tmpn2=length(tmpmx_nIMF2{kk});
                mx_nIMF2(k,kk,1:tmpn2)=tmpmx_nIMF2{kk};
                if phase
                    trial_imfphdata(kk,1:size(ptrial_imfphdata{kk},2),1:size(ptrial_imfphdata{kk},1)) = ptrial_imfphdata{kk}';
                end
                if phase2
                    trial_IMFPHdata(kk,1:size(ptrial_IMFPHdata{kk},2),1:size(ptrial_IMFPHdata{kk},3),1:dNsamples) = permute(spmm_downsample(ptrial_IMFPHdata{kk},sdownLV),[2,3,1]);
                end
            end

        end
        Dimf(:,:,:,k)=trial_imfdata;
        Dfm(:,:,:,k)=trial_imfmdata;
        Dam(:,:,:,k)=trial_imamdata;
        DIMF(:,:,:,k)=reshape(trial_IMFdata,Nchannels,maxIMF*maximf, dNsamples);
        DFM(:,:,:,k)=reshape(trial_IMFMdata,Nchannels,maxIMF*maximf, dNsamples);
        DAM(:,:,:,k)=reshape(trial_IMAMdata,Nchannels,maxIMF*maximf, dNsamples);
        if phase
            Dimph(:,:,:,k)=trial_imfphdata;
        end
        if phase2
            DIMPH(:,:,:,k)=reshape(trial_IMFPHdata,Nchannels,maxIMF*maximf, dNsamples);
        end

end
%-Save new M/EEG dataset(s)
 Dimf.nimf=vec_nimf;
 Dimf = Dimf.history('maskemd2layer_tfpar_SV', S);
 save(Dimf);
    %%%%%%%%%%%%%%%%%%
 Dfm.nimf=vec_nimf;
 Dfm = Dfm.history('maskemd2layer_tfpar_SV', S);
 save(Dfm);
    %%%%%%%%%%%%%%%%%%
 Dam.nimf=vec_nimf;
 Dam = Dam.history('maskemd2layer_tfpar_SV', S);
 save(Dam);
%  Dtf.nimf=vec_nimf;
%  save(Dtf);
 if S.phase
        Dimph.nimf=vec_nimf;
        Dimph = Dimph.history('maskemd2layer_tfpar_SV', S);
        save(Dimph);
 end
     
    
 %%%%%%%%%%%%%%%%%%
DIMF.nimf=vec_nimf;
DIMF.maximf=maximf;
DIMF.maxIMF=maxIMF;
DIMF.nIMF=vec_nIMF2;
DIMF.mnIMF=mx_nIMF2;
DIMF = DIMF.history('maskemd2layer_tfpar_SV', S);
save(DIMF);
%%%%%%%%%%%%%%%%%%
DFM.nimf=vec_nimf;
DFM.maximf=maximf;
DFM.maxIMF=maxIMF;
DFM.nIMF=vec_nIMF2;
DFM.mnIMF=mx_nIMF2;
DFM = DFM.history('maskemd2layer_tfpar_SV', S);
save(DFM);
%%%%%%%%%%%%%%%%%%%
DAM.nimf=vec_nimf;
DAM.maximf=maximf;
DAM.maxIMF=maxIMF;
DAM.nIMF=vec_nIMF2;
DAM.mnIMF=mx_nIMF2;
DAM = DAM.history('maskemd2layer_tfpar_SV', S);
save(DAM);
if phase2
    DIMPH.nimf=vec_nimf;
    DIMPH.maximf=maximf;
    DIMPH.maxIMF=maxIMF;
    DIMPH.nIMF=vec_nIMF2;
    DIMPH.mnIMF=mx_nIMF2;
    DIMPH = DIMPH.history('maskemd2layer_tfpar_SV', S);
    save(DIMPH);
end
%--------------------------------------------------------------------------



