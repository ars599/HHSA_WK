function D = spmm_eeg_average_LargeTF(S)
% Quick averaging each channel over trials or trial types, for time-frequency data
% The code is adapted from  spm_eeg_average_TF, allowing parallel averaging
% by blocks of channels
% WK Liang 2019
% Original FORMAT D = spm_eeg_average_TF(S)
%
% S         - optional input struct
% (optional) fields of S:
% S.D           - MEEG object or filename of M/EEG mat-file with epoched TF data
% S.circularise - flag that indicates whether average is straight (0) or
%                 vector (1) of phase angles.
% S.robust      - (optional) - use robust averaging (only for power)
%                 .savew  - save the weights in an additional dataset
%                 .bycondition - compute the weights by condition (1,
%                                default) or from all trials (0)
%                 .ks     - offset of the weighting function (default: 3)
%
% Output:
% D         - MEEG object (also written to disk)
%__________________________________________________________________________
%
% This function averages single trial time-frequency data within trial type.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_average_TF.m 3616 2009-12-08 15:16:39Z vladimir $

SVNrev = '$Rev: 3616 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG TF averaging'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end
if  ~isfield(S,'blocksize')
    nchan_blocksize=100;
else
    nchan_blocksize=S.blocksize;
end
D = spm_eeg_load(D);

%-Check data type
%--------------------------------------------------------------------------
if ~strcmp(D.type, 'single')
    error('This function can only be applied to single trial data');
elseif ~strncmp(D.transformtype, 'TF', 2) % TF and TFphase
    error('This function can only be applied to time-frequency data.');
end

Dnew = clone(D, ['m' fname(D)], [D.nchannels D.nfrequencies D.nsamples D.nconditions]);
    cl   = D.condlist;
    ni = zeros(1,D.nconditions);
    for i = 1:D.nconditions
        %w = pickconditions(D, deblank(cl{i}), 1)';
        w = indtrial(D, deblank(cl{i}), 'GOOD')';
        ni(i) = length(w);
        if ni(i) == 0
            warning('%s: No trials for trial type %d', D.fname, cl{i});
        end
    end

    %goodtrials  =  pickconditions(D, cl, 1);
    nc=D.nchannels;
    nblocks=ceil(nc/nchan_blocksize);
    for j = 1:nblocks
        if j< nblocks
           id_j=((j-1)*nchan_blocksize+1):j*nchan_blocksize;
        else
           id_j=((j-1)*nchan_blocksize+1): nc;
        end
        for i = 1:D.nconditions

            %w = pickconditions(D, deblank(cl{i}), 1)';
            w = indtrial(D, deblank(cl{i}), 'GOOD')';
            if isempty(w)
                continue;
            end

            %-Straight average
            %------------------------------------------------------------------
            if ~strcmp(D.transformtype, 'TFphase')
                    Dnew(id_j, :, :, i) = mean(D(id_j, :, :, w), 4);
                %-Vector average (eg PLV for phase)
                %------------------------------------------------------------------
            else
                tmp = D(id_j, :, :, w);
                tmp = exp(sqrt(-1)*tmp);
%                 if plv
                    Dnew(id_j, :, :, i) = abs(mean(tmp,4));
%                 else
%                     Dnew(j, :, :, i) = angle(mean(tmp,4));
%                 end
            end
        end
        fprintf('finish averaging for channel block %d\n',j);
    end


    Dnew = type(Dnew, 'evoked');

    %-Update some header information
    %--------------------------------------------------------------------------
    Dnew = conditions(Dnew, ':', cl);
    Dnew = repl(Dnew, ':', ni);

    %-Display averaging statistics
    %--------------------------------------------------------------------------
    disp(sprintf('%s: Number of replications per contrast:', Dnew.fname));  %-#
    s = [];
    for i = 1:D.nconditions
        s = [s sprintf('average %s: %d trials', cl{i}, ni(i))];
        if i < D.nconditions
            s = [s sprintf(', ')];
        else
            s = [s '\n'];
        end
    end
    disp(sprintf(s));                                                       %-#

    %-Save new evoked M/EEG dataset
    %--------------------------------------------------------------------------
    Dnew = Dnew.history(mfilename, S);
    save(Dnew);
D = Dnew;
%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG TF averaging: done'); spm('Pointer', 'Arrow');
