function [D] = spmm_eeg_tf_rescale(S)
% Rescale spectrogram, modified by W-K Liang from spm_eeg_tf_rescale.m in SPM
% The modification is to add regularization parameter to denominator
% Because of the high resolution of frequency in HHSA/HHT, there are lots
% of zeros in spectrum
% FORMAT [D] = spm_eeg_tf_rescale(S)
%
% S                    - input structure (optional)
% (optional) fields of S:
%   S.D                - MEEG object or filename of M/EEG mat-file
%   S.tf               - structure with (optional) fields:
%     S.tf.method      - 'LogR', 'Diff', 'Rel', 'Log', 'Sqrt'
%     S.tf.Sbaseline   - 2-element vector: start and stop of baseline 
%                        (need to specify this for LogR and Diff)
%     S.tf.Db          - MEEG object or filename of M/EEG mat-file to use
%                        for the baseline (if different from the input dataset).
% 
% D                    - MEEG object with rescaled power data (also
%                        written to disk with prefix r)
%
% For 'Log' and 'Sqrt', these functions are applied to spectrogram 
% For 'LogR', 'Rel' and 'Diff' this function computes power in the baseline
% p_b and outputs (i) p-p_b for 'Diff' (ii) 100*(p-p_b)/p_b for 'Rel' 
%                 (iii) log (p/p_b) for 'LogR'
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_eeg_tf_rescale.m 4316 2011-04-26 16:52:28Z vladimir $

SVNrev = '$Rev: 4316 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG Time-Frequency Rescale'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

try
    S.tf.method;
catch
    str  = {'LogR','Diff', 'Rel', 'Log', 'Sqrt', 'Zscore'};
    S.tf.method = spm_input('Rescale method','+1','m',str,char(str),1);
end

Din  = spm_eeg_load(D);
if ~isempty(strfind(Din.fnamedat,'ierp'))
    ierp=1;
else
    ierp=0;
end

tims = time(Din);
if ~isempty(strfind(Din.fname,'omegaf'))
    fF=1;
else
    fF=0;
end
Nf   = length(frequencies(Din));
if ~strcmpi(S.tf.method,'log')
  D    = clone(Din, ['r' Din.fname], [Din.nchannels Nf Din.nsamples Din.ntrials]);
else
    D    = clone(Din, ['rlg' Din.fname], [Din.nchannels Nf Din.nsamples Din.ntrials]);
end
try
    regu_param=S.tf.regu_param;
catch
    regu_param=1e-8;
end
try
    smooth=S.tf.smooth;
catch
    smooth=0;
end

q_pre=fspecial('gaussian', [5,5],0.5);
%q_pre2=fspecial('average', [3,3]);
if ierp
    q_pre2=fspecial('gaussian', [1,3],0.5);
else
    q_pre2=fspecial('gaussian', [1,7],0.75);
end
q_pre3=fspecial('gaussian', [3,3],0.6);
q_pre4=fspecial('gaussian', [3,7],0.75);
if smooth
    f1_ord=S.tf.f1_ord;
    f2_ord=S.tf.f2_ord;
    f3_ord=S.tf.f3_ord;
    f3_roi=S.tf.f3_roi;
    f4_ord=S.tf.f4_ord;
else
    f1_ord=0;
    f2_ord=0;
    f3_ord=0;
    f3_roi=[];
    f4_ord=0;
end
switch lower(S.tf.method)
    
    case {'logr','diff', 'rel', 'zscore'}
        regu_pwr=regu_param*nanmean(Din(:));
        try
            S.tf.Sbaseline;
        catch            
            if spm_input('Baseline dataset','+1','b',{'Same|Different'},[0 1],0)
                [Db, sts] = spm_select(1, 'mat', 'Select baseline M/EEG mat file');
                if ~sts, return; end
                S.tf.Db = Db;
            else
                S.tf.Db = [];
            end
            
            tmp_base = spm_input('Start and stop of baseline [ms]', '+1', 'i', '', 2);
            S.tf.Sbaseline = tmp_base/1000;
        end
        
        if isfield(S.tf, 'Db') && ~isempty(S.tf.Db)
            Db = spm_eeg_load(S.tf.Db);
        else
            Db = Din;
        end
        
        if any(abs(Din.frequencies-Db.frequencies)>0.1) || ~isequal(Db.chanlabels, Din.chanlabels) ||...
                (Db.ntrials>1 && (Db.ntrials~=Din.ntrials))
            error('The input dataset and the baseline dataset should have the same frequencies, channels and trial numbers');
        end        
        if ~fF
            for c=1:D.ntrials
                inds=find(tims>=S.tf.Sbaseline(1) & tims<=S.tf.Sbaseline(2));
                x=spm_squeeze(Din(:,:,:,c), 4);
                if Db.ntrials > 1
                    xbaseC=spm_squeeze(Db(:,:,:,c), 4);
                else
                    xbaseC=spm_squeeze(Db(:,:,:,1), 4);
                end
                if smooth
                    for ch=1: nchannels(D)
                        nts=shiftdim(x(ch,:,:),1);
                        ntb=shiftdim(xbaseC(ch,:,:),1);
                        if f1_ord>=1
                            for if1=1:f1_ord
                                nts=filter2(q_pre,nts);
                                ntb=filter2(q_pre,ntb);
                            end
                        end
                        if f2_ord>=1
                            for if2=1:f2_ord
                                nts=filter2(q_pre2,nts);
                                ntb=filter2(q_pre2,ntb);
                            end
                        end
                        if f3_ord>=1
                            if isempty(f3_roi)
                                for if1=1:f3_ord
                                    nts=filter2(q_pre3,nts);
                                    ntb=filter2(q_pre3,ntb);
                                end
                            else
                                if length(f3_roi)==1 && f3_roi>=1 && f3_roi<=size(nts,2)
                                    f3_roi=[1 round(f3_roi)];
                                elseif length(f3_roi)>1 && f3_roi(1)>=1 && f3_roi(2)<=size(nts,2) && f3_roi(1)<f3_roi(2)
                                    f3_roi=[f3_roi(1), f3_roi(2)];
                                    f3_roi=round(f3_roi);
                                else
                                    f3_roi=[1 size(nts,2)];
                                end
                                f3_roid=f3_roi(1):f3_roi(2);
                                for if1=1:f3_ord
                                    nts(:,f3_roid)=filter2(q_pre3,nts(:,f3_roid));
                                    ntb(:,f3_roid)=filter2(q_pre3,ntb(:,f3_roid));
                                end
                            end
                        end
                        if f4_ord>=1
                            for if2=1:f4_ord
                                nts=filter2(q_pre4,nts);
                                ntb=filter2(q_pre4,ntb);
                            end
                        end
                        x(ch,:,:)=shiftdim(nts,-1);
                        xbaseC(ch,:,:)=shiftdim(ntb,-1);
                    end
                end
                switch lower(S.tf.method)
                    case 'logr'
                        %xbase=mean(log10(xbase(:,:,inds)+regu_pwr),3);
                        xbase=log10(mean(xbaseC(:,:,inds)+regu_pwr,3));
                        D(:,:,:,c)= 10*(log10(x) - repmat(xbase,[1 1 D.nsamples 1]));
                        D = units(D, [], 'dB');
                    case 'diff'
                        xbase=mean(xbaseC(:,:,inds),3);
                        D(:,:,:,c)= (x - repmat(xbase,[1 1 D.nsamples 1]));
                    case 'zscore'
                        stdev = std(xbaseC(:,:,inds), [], 3)+regu_pwr;
                        xbase= mean(xbaseC(:,:,inds),3);                    
                        D(:,:,:,c)= (x - repmat(xbase,[1 1 D.nsamples 1]))./repmat(stdev,[1 1 D.nsamples 1]);
                    case 'rel'
                        xbase=mean(xbaseC(:,:,inds),3);
                        %D(:,:,:,c)= 100*((x./(repmat(xbase,[1 1 D.nsamples 1])+regu_pwr) - 1));
                        xbase_ex=repmat(xbase,[1 1 D.nsamples 1]);
                        D(:,:,:,c)= 100*(x-xbase_ex)./(xbase_ex+regu_pwr);
                        D = units(D, [], '%');
                end
            end
        else
            for c=1:D.ntrials
                    %inds=find(tims>=S.tf.Sbaseline(1) & tims<=S.tf.Sbaseline(2));
                    x=spm_squeeze(Din(:,:,:,c), 4);
                    if Db.ntrials > 1
                        xbaseC=spm_squeeze(Db(:,:,:,c), 4);
                    else
                        xbaseC=spm_squeeze(Db(:,:,:,1), 4);
                    end
                    if smooth
                        for ch=1: nchannels(D)
                            nts=shiftdim(x(ch,:,:),1);
                            ntb=shiftdim(xbaseC(ch,:,:),1);
                            if f1_ord>=1
                                for if1=1:f1_ord
                                    nts=filter2(q_pre,nts);
                                    ntb=filter2(q_pre,ntb);
                                end
                            end
                            if f2_ord>=1
                                for if2=1:f2_ord
                                    nts=filter2(q_pre2,nts);
                                    ntb=filter2(q_pre2,ntb);
                                end
                            end
                            if f3_ord>=1
                                if isempty(f3_roi)
                                    for if1=1:f3_ord
                                        nts=filter2(q_pre3,nts);
                                        ntb=filter2(q_pre3,ntb);
                                    end
                                else
                                    if length(f3_roi)==1 && f3_roi>=1 && f3_roi<=size(nts,2)
                                        f3_roi=[1 round(f3_roi)];
                                    elseif length(f3_roi)>1 && f3_roi(1)>=1 && f3_roi(2)<=size(nts,2) && f3_roi(1)<f3_roi(2)
                                        f3_roi=[f3_roi(1), f3_roi(2)];
                                        f3_roi=round(f3_roi);
                                    else
                                        f3_roi=[1 size(nts,2)];
                                    end
                                    f3_roid=f3_roi(1):f3_roi(2);
                                    for if1=1:f3_ord
                                        nts(:,f3_roid)=filter2(q_pre3,nts(:,f3_roid));
                                        ntb(:,f3_roid)=filter2(q_pre3,ntb(:,f3_roid));
                                    end
                                end
                            end
                            if f4_ord>=1
                                for if2=1:f4_ord
                                    nts=filter2(q_pre4,nts);
                                    ntb=filter2(q_pre4,ntb);
                                end
                            end
                            x(ch,:,:)=shiftdim(nts,-1);
                            xbaseC(ch,:,:)=shiftdim(ntb,-1);
                        end
                    end
                    switch lower(S.tf.method)
                        case 'logr'
                            xbase=log10(xbaseC + regu_pwr);
                            D(:,:,:,c)= 10*(log10(x) - xbase);
                            
                            D = units(D, [], 'dB');
                        case 'diff'
                            xbase=xbaseC(:,:,:);
                            D(:,:,:,c)= (x - xbase);
%                         case 'zscore'
%                             stdev = std(xbaseC(:,:,inds), [], 3);
%                             xbase= mean(xbaseC(:,:,inds),3);                    
%                             D(:,:,:,c)= (x - repmat(xbase,[1 1 D.nsamples 1]))./repmat(stdev,[1 1 D.nsamples 1]);
                        case 'rel'
                            xbase=xbaseC(:,:,:);
                            %D(:,:,:,c)= 100*(x./(xbase + regu_pwr) - 1);
                            D(:,:,:,c)= 100*(x - xbase)./(xbase + regu_pwr);
                            D = units(D, [], '%');
                    end
             end
            
        end
        
    case 'log'
            try
                rp=S.tf.regu_param;
            catch
                rp=1e-8;
            end
            rp_pwr=rp*nanmean(Din(:));
            for c=1:D.ntrials
                tr_data=Din(:,:,:,c);
                D(:,:,:,c) = log(tr_data+rp_pwr*ones(size(tr_data)));
            end
        
    case 'sqrt'
        for c=1:D.ntrials
            D(:,:,:,c) = sqrt(Din(:,:,:,c));
        end
        
    otherwise
        error('Unknown rescaling method.');
end

% Save
D = D.history(mfilename, S);
D.rescale_type=lower(S.tf.method);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG Time Frequency Rescale: done'); 
spm('Pointer','Arrow');
