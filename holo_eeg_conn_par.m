function [Dplv] = holo_eeg_conn_par(S)
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


if nargin == 0
    S = [];
end
if ~isfield(S, 'parallel')
    S.parallel = 0;
end
if ~isfield(S, 'modu_type') || isempty(S.modu_type)
    modu_type = 1;
else
    modu_type=S.modu_type;
end
if ~isfield(S, 'conn_type') || isempty(S.conn_type)
    conn_type = 1;
else
    conn_type=S.conn_type;
end
if ~isfield(S, 'par_seed') || isempty(S.par_seed)
    par_seed = 0;
else
    par_seed=S.par_seed;
end
if ~isfield(S, 'seed')
    seed_conn=0;
else
    seed=S.seed;
    seed=seed(seed>0); % support multi seeds
%     if length(seed)>1
%         seed=seed(1);
%     end
    if isempty(seed)
        seed_conn=0;
    else
        seed_conn=1;
    end
end
if ~isfield(S, 'CoreN')
    S.CoreN = feature('numcores') - 1;
end
try
    sender=S.sender;
    receiver=S.receiver;
catch
    disp('this version require you enter a single IMF as the sender, as well as the receiver!');
    return
end
if length(sender)>1
    sender =sender(1);
end
% From 20190118, it can be multiple receivers, e.g., receiver=[4,5,6].
% However, sender can be only one till now.
% if length(receiver)>1
%     receiver =receiver(1);
% end
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
if seed_conn && par_seed
    if nd==3
        try
            parname=strrep(name,'src','par');
            Dfile_fm=strcat(pathstr,sep,'timfm_',name,ext);
            Dfile_fm2=strcat(pathstr,sep,'timfm2_',parname,ext); % second layer par file
            Dfile_ph=strcat(pathstr,sep,'timph_',name,ext);
            Dfile_ph2=strcat(pathstr,sep,'timph2_',parname,ext); % second layer par ph file
            dname=D.fname;
        catch
            error('selected file is not a source file or does not have par file!'); 
        end 
    elseif nd==4
        dname=D.fname;
        dname=dname(6:end);
        try
            parname=strrep(dname,'src','par');
            Dfile_fm=strcat(pathstr,sep,'timfm_',dname);
            Dfile_fm2=strcat(pathstr,sep,'timfm2_',parname); % second layer par file
            Dfile_ph=strcat(pathstr,sep,'timph_',dname);
            Dfile_ph2=strcat(pathstr,sep,'timph2_',parname); % second layer par ph file
        catch
            error('unknown filetype'); 
        end 
    else
       error('unknown filetype'); 
    end

else
    if nd==3   
        Dfile_fm=strcat(pathstr,sep,'timfm_',name,ext);
        Dfile_fm2=strcat(pathstr,sep,'timfm2_',name,ext);
        Dfile_ph=strcat(pathstr,sep,'timph_',name,ext);
        Dfile_ph2=strcat(pathstr,sep,'timph2_',name,ext);
        dname=D.fname;
    elseif nd==4
        dname=D.fname;
        dname=dname(6:end);
        Dfile_fm=strcat(pathstr,sep,'timfm_',dname);
        Dfile_fm2=strcat(pathstr,sep,'timfm2_',dname);
        Dfile_ph=strcat(pathstr,sep,'timph_',dname);
        Dfile_ph2=strcat(pathstr,sep,'timph2_',dname);
    else
       error('unknown filetype'); 
    end
end
%D=spm_eeg_load(Dfile);
try
    D_fm=spm_eeg_load(Dfile_fm);
    D_fm2=spm_eeg_load(Dfile_fm2);
    D_ph=spm_eeg_load(Dfile_ph);
    D_ph2=spm_eeg_load(Dfile_ph2);
catch
    disp('2 layer EEMD files missing'); 
    return;
end
try
    if D_fm2.sdownLV>0
        sdownLV=round(sdownLV);
    else
        sdownLV=0;
    end
catch
    sdownLV=0;
end
% if S.dyadic==1
%    dyadic=1;
% else
%     dyadic=0;
% end

%-Configure the analysis
%--------------------------------------------------------------------------
timesol=S.settings.timesol;
half_winsize=S.settings.half_winsize;
winsize=2*half_winsize+1;
winlen=1000*(winsize-1)/D_fm2.fsample;

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
% S.frequencies = fscale;
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
timeind = D_fm2.indsample(1e-3*min(S.timewin)):D_fm2.indsample(1e-3*max(S.timewin));
if sdownLV>0
    tind_all=1:sdownLV:nsamples(D_fm);
    timeindf=tind_all(timeind);
else
    timeindf=timeind;
end
Nchannels = length(chanind);
if seed_conn
    Nconn=Nchannels;
else
    Nconn=Nchannels^2;
end

%%%For holo conn
if ~seed_conn
  for irv=1:length(receiver)
        hi1=int2str(sender); hi2=int2str(receiver(irv));
        switch conn_type
           case 1
               switch modu_type
                   case 1
                           Dplv = clone(D, ['tconnplv_',hi1,hi2, 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_FM_','L',num2str(winlen),'_', dname], [Nconn 1 timesol D.ntrials]);
                   case 2
                           Dplv = clone(D, ['tconnplv_',hi1,hi2, 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_HM_','L',num2str(winlen),'_', dname], [Nconn 1 timesol D.ntrials]);
               end
               Dplv = transformtype(Dplv, 'TF');  % ISPC-time values are no longer "phase angle" 
           case 2
               switch modu_type
                   case 1
                           Dplv = clone(D, ['tconnplph_',hi1,hi2, 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_FM_','L',num2str(winlen),'_', dname], [Nconn 1 timesol D.ntrials]);
                   case 2
                           Dplv = clone(D, ['tconnplph_',hi1,hi2, 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_HM_','L',num2str(winlen),'_', dname], [Nconn 1 timesol D.ntrials]);
               end
               Dplv = transformtype(Dplv, 'TFphase'); % ISPC-trial values before average are "phase angle"
        end
        Dplv = Dplv.frequencies(:, sender);
        Dplv = timeonset(Dplv, 1e-3*min(S.timewin));
        Dplv = fsample(Dplv, timesol/(1e-3*max(S.timewin)-1e-3*min(S.timewin)));
        conn_idid=zeros(Nconn,2);
        conn_labels=cell(1,Nconn);
        for i=1:Nchannels
           for j=1:Nchannels
               conn_idid((i-1)*Nchannels+j,:)=[i,j];
               conn_labels{1,(i-1)*Nchannels+j}=strcat('conn','_',int2str(i),'_',int2str(j));
           end
        end
        Dplv = chanlabels(Dplv, 1:Nconn, conn_labels);
        Dplv = chantype(Dplv, 1:Nconn, 'LFP');
        Dplv.conn_idid=conn_idid;
        % vec_nimf=D_fm2.nimf;
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
        w=indtrial(D,deblank(cl{i}));
        %btrials=indtrial(D,deblank(cl{i}),'BAD');

        if sender== receiver(irv) % this is actually the conventional case that connections are within the same frequency range
            for k = w
                %fm_trialdata=D_fm(:,sender,timeindf,k);
                %fm_trialdata=permute(fm_trialdata,[3,1,2]);
                %pfm_trialdata=mat2cell(fm_trialdata,size(fm_trialdata,1),ones(1,size(fm_trialdata,2)));
%                 fm2_trialdata=D_fm2(:,:,timeind,k);
%                 fm2_trialdata=reshape(fm2_trialdata,size(fm2_trialdata,1),max_IMF,max_imf,[]);
%                 fm2_trialdata=permute(fm2_trialdata,[4,2,1,3]);
%                 fm2_receiver_trialdata=fm2_trialdata(:,:,:,receiver(irv));
%                 pfm2_trialdata=mat2cell(fm2_receiver_trialdata,size(fm2_receiver_trialdata,1),size(fm2_receiver_trialdata,2),ones(1,size(fm2_receiver_trialdata,3)));
                ph_trialdata=D_ph(:,sender,timeindf,k);
                ph_trialdata=permute(ph_trialdata,[3,1,2]);
                pph_trialdata=mat2cell(ph_trialdata,size(ph_trialdata,1),ones(1,size(ph_trialdata,2)));
%                 ph2_trialdata=D_ph2(:,:,timeind,k);
%                 ph2_trialdata=reshape(ph2_trialdata,size(ph2_trialdata,1),max_IMF,max_imf,[]);
%                 ph2_trialdata=permute(ph2_trialdata,[4,2,1,3]);
%                 ph2_receiver_trialdata=ph2_trialdata(:,:,:,receiver(irv));
%                 pph2_trialdata=mat2cell(ph2_receiver_trialdata,size(ph2_receiver_trialdata,1),size(ph2_receiver_trialdata,2),ones(1,size(ph2_receiver_trialdata,3)));
%                 pcurr_nIMF = vec_nIMF(k,:);
                
                   plv_trialdata=zeros(Nconn,1,timesol); 
                   %plph_trialdata=zeros(Nconn,1,timesol);
                 if modu_type==2 && conn_type==1
                     parfor icon=1:Nconn  
                        senrec=conn_idid(icon,:);
                        inich=senrec(1);
                        endch=senrec(2);
%                        curr_nIMF = pcurr_nIMF(endch);
                        [plv,~]=nspplotf3d_tres3_connLinear(2*pph_trialdata{inich}(:,1),2*pph_trialdata{endch}(:,1),timesol,half_winsize,S);
                        plv_trialdata(icon,1,:)=shiftdim(plv,-1);
                       %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                     end
                 elseif modu_type==1 && conn_type==1
                     parfor icon=1:Nconn  
                       senrec=conn_idid(icon,:);
                       inich=senrec(1);
                       endch=senrec(2);
%                        curr_nIMF = pcurr_nIMF(endch);
                       [plv,~]=nspplotf3d_tres3_connLinear(pph_trialdata{inich}(:,1),pph_trialdata{endch}(:,1),timesol,half_winsize,S);
                       plv_trialdata(icon,1,:)=shiftdim(plv,-1);
                       %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                     end
                 elseif modu_type==2 && conn_type==2
                     parfor icon=1:Nconn  
                        senrec=conn_idid(icon,:);
                        inich=senrec(1);
                        endch=senrec(2);
%                        curr_nIMF = pcurr_nIMF(endch);
                        [~,plph]=nspplotf3d_tres3_connLinear(2*pph_trialdata{inich}(:,1),2*pph_trialdata{endch}(:,1),timesol,half_winsize,S);
                        plv_trialdata(icon,1,:)=shiftdim(plph,-1);
                       %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                     end
                 elseif modu_type==1 && conn_type==2
                     parfor icon=1:Nconn  
                       senrec=conn_idid(icon,:);
                       inich=senrec(1);
                       endch=senrec(2);
%                        curr_nIMF = pcurr_nIMF(endch);
                       [~,plph]=nspplotf3d_tres3_connLinear(pph_trialdata{inich}(:,1),pph_trialdata{endch}(:,1),timesol,half_winsize,S);
                       plv_trialdata(icon,1,:)=shiftdim(plph,-1);
                       %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                     end
                 end

                   Dplv(:,:,:,k)=plv_trialdata;
                   %Dplph(:,:,:,k)=plph_trialdata;
            end % end k=w        
        else
            for k = w
                fm_trialdata=D_fm(:,sender,timeindf,k);
                fm_trialdata=permute(fm_trialdata,[3,1,2]);
                pfm_trialdata=mat2cell(fm_trialdata,size(fm_trialdata,1),ones(1,size(fm_trialdata,2)));
                fm2_trialdata=D_fm2(:,:,timeind,k);
                fm2_trialdata=reshape(fm2_trialdata,size(fm2_trialdata,1),max_IMF,max_imf,[]);
                fm2_trialdata=permute(fm2_trialdata,[4,2,1,3]);
                fm2_receiver_trialdata=fm2_trialdata(:,:,:,receiver(irv));
                pfm2_trialdata=mat2cell(fm2_receiver_trialdata,size(fm2_receiver_trialdata,1),size(fm2_receiver_trialdata,2),ones(1,size(fm2_receiver_trialdata,3)));
                ph_trialdata=D_ph(:,sender,timeindf,k);
                ph_trialdata=permute(ph_trialdata,[3,1,2]);
                pph_trialdata=mat2cell(ph_trialdata,size(ph_trialdata,1),ones(1,size(ph_trialdata,2)));
                ph2_trialdata=D_ph2(:,:,timeind,k);
                ph2_trialdata=reshape(ph2_trialdata,size(ph2_trialdata,1),max_IMF,max_imf,[]);
                ph2_trialdata=permute(ph2_trialdata,[4,2,1,3]);
                ph2_receiver_trialdata=ph2_trialdata(:,:,:,receiver(irv));
                pph2_trialdata=mat2cell(ph2_receiver_trialdata,size(ph2_receiver_trialdata,1),size(ph2_receiver_trialdata,2),ones(1,size(ph2_receiver_trialdata,3)));
                pcurr_nIMF = vec_nIMF(k,:);
                plv_trialdata=zeros(Nconn,1,timesol); 
                   %plph_trialdata=zeros(Nconn,1,timesol);
                if modu_type==2 && conn_type==1
                   parfor icon=1:Nconn  
                       senrec=conn_idid(icon,:);
                       inich=senrec(1);
                       endch=senrec(2);
                       curr_nIMF = pcurr_nIMF(endch);
                       [plv,~]=nspplotf3d_tres3_conn(2*pfm_trialdata{inich}(:,1),2*pph_trialdata{inich}(:,1),pfm2_trialdata{endch}(:,1:curr_nIMF),pph2_trialdata{endch}(:,1:curr_nIMF),timesol,half_winsize,S);
                       plv_trialdata(icon,1,:)=shiftdim(plv,-1);
                       %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                   end
                elseif modu_type==1 && conn_type==1
                     parfor icon=1:Nconn  
                       senrec=conn_idid(icon,:);
                       inich=senrec(1);
                       endch=senrec(2);
                       curr_nIMF = pcurr_nIMF(endch);
                       [plv,~]=nspplotf3d_tres3_conn(pfm_trialdata{inich}(:,1),pph_trialdata{inich}(:,1),pfm2_trialdata{endch}(:,1:curr_nIMF),pph2_trialdata{endch}(:,1:curr_nIMF),timesol,half_winsize,S);
                       plv_trialdata(icon,1,:)=shiftdim(plv,-1);
                       %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                     end
                elseif modu_type==2 && conn_type==2
                   parfor icon=1:Nconn  
                       senrec=conn_idid(icon,:);
                       inich=senrec(1);
                       endch=senrec(2);
                       curr_nIMF = pcurr_nIMF(endch);
                       [~,plph]=nspplotf3d_tres3_conn(2*pfm_trialdata{inich}(:,1),2*pph_trialdata{inich}(:,1),pfm2_trialdata{endch}(:,1:curr_nIMF),pph2_trialdata{endch}(:,1:curr_nIMF),timesol,half_winsize,S);
                       plv_trialdata(icon,1,:)=shiftdim(plph,-1);
                       %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                   end
                elseif modu_type==1 && conn_type==2
                     parfor icon=1:Nconn  
                       senrec=conn_idid(icon,:);
                       inich=senrec(1);
                       endch=senrec(2);
                       curr_nIMF = pcurr_nIMF(endch);
                       [~,plph]=nspplotf3d_tres3_conn(pfm_trialdata{inich}(:,1),pph_trialdata{inich}(:,1),pfm2_trialdata{endch}(:,1:curr_nIMF),pph2_trialdata{endch}(:,1:curr_nIMF),timesol,half_winsize,S);
                       plv_trialdata(icon,1,:)=shiftdim(plph,-1);
                       %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                     end
                end

                   Dplv(:,:,:,k)=plv_trialdata;
                   %Dplph(:,:,:,k)=plph_trialdata;
            end % end k=w
        end
    end

    %-Save new M/EEG dataset(s)
    Dplv = Dplv.history(mfilename, S);
    save(Dplv);
    %save(Dplph);

  end  % for multiple receivers
else % for seed_conn
   for s=1:length(seed) 
         for irv=1:length(receiver)
               hi1=int2str(sender); hi2=int2str(receiver(irv));
               switch conn_type
                   case 1
                       switch modu_type
                           case 1
                                   Dplv = clone(D, ['tconnplv_',hi1,hi2,'seed',int2str(seed(s)),'_', 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_FM_','L',num2str(winlen),'_', dname], [Nconn 1 timesol D.ntrials]);
                           case 2
                                   Dplv = clone(D, ['tconnplv_',hi1,hi2,'seed',int2str(seed(s)),'_', 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_HM_','L',num2str(winlen),'_', dname], [Nconn 1 timesol D.ntrials]);
                       end
                       Dplv = transformtype(Dplv, 'TF');  % ISPC-time values are no longer "phase angle" 
                   case 2
                       switch modu_type
                           case 1
                                   Dplv = clone(D, ['tconnplph_',hi1,hi2,'seed',int2str(seed(s)),'_', 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_FM_','L',num2str(winlen),'_', dname], [Nconn 1 timesol D.ntrials]);
                           case 2
                                   Dplv = clone(D, ['tconnplph_',hi1,hi2,'seed',int2str(seed(s)),'_', 't',int2str(min(S.timewin)),'_',int2str(max(S.timewin)),'_HM_','L',num2str(winlen),'_', dname], [Nconn 1 timesol D.ntrials]);
                       end
                       Dplv = transformtype(Dplv, 'TFphase'); % ISPC-trial values before average are "phase angle"
               end
               Dplv = Dplv.frequencies(:, sender);
               Dplv = timeonset(Dplv, 1e-3*min(S.timewin));
               Dplv = fsample(Dplv, timesol/(1e-3*max(S.timewin)-1e-3*min(S.timewin)));
               Dplv = chanlabels(Dplv, 1:Nconn, D.chanlabels(chanind));
               Dplv = chantype(Dplv, 1:Nconn, D.chantype(chanind));
               Dplv.seed=seed(s);
               % vec_nimf=D_fm2.nimf;
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
                w=indtrial(D,deblank(cl{i}));
                %btrials=indtrial(D,deblank(cl{i}),'BAD');
                if sender== receiver(irv)
                    for k = w
%                         fm_trialdata=D_fm(:,sender,timeindf,k);
%                         fm_trialdata=permute(fm_trialdata,[3,1,2]);
%                         pfm_trialdata=mat2cell(fm_trialdata,size(fm_trialdata,1),ones(1,size(fm_trialdata,2)));
%                         fm2_trialdata=D_fm2(seed(s),:,timeind,k);
%                         fm2_trialdata=reshape(fm2_trialdata,size(fm2_trialdata,1),max_IMF,max_imf,[]);
%                         fm2_trialdata=permute(fm2_trialdata,[4,2,1,3]);
%                         fm2_receiver_trialdata=fm2_trialdata(:,:,:,receiver(irv));
%                         pfm2_trialdata=mat2cell(fm2_receiver_trialdata,size(fm2_receiver_trialdata,1),size(fm2_receiver_trialdata,2));
                        ph_trialdata=D_ph(:,sender,timeindf,k);
                        ph_trialdata=permute(ph_trialdata,[3,1,2]);
                        pph_trialdata=mat2cell(ph_trialdata,size(ph_trialdata,1),ones(1,size(ph_trialdata,2)));
%                         ph2_trialdata=D_ph2(seed(s),:,timeind,k);
%                         ph2_trialdata=reshape(ph2_trialdata,size(ph2_trialdata,1),max_IMF,max_imf,[]);
%                         ph2_trialdata=permute(ph2_trialdata,[4,2,1,3]);
%                         ph2_receiver_trialdata=ph2_trialdata(:,:,:,receiver(irv));
%                         pph2_trialdata=mat2cell(ph2_receiver_trialdata,size(ph2_receiver_trialdata,1),size(ph2_receiver_trialdata,2));
%                         pcurr_nIMF = vec_nIMF(k,:);
%                         curr_nIMF = pcurr_nIMF(seed(s));
%                         seedfm2_trialdata=pfm2_trialdata{1}(:,1:curr_nIMF);
%                         seedph2_trialdata=pph2_trialdata{1}(:,1:curr_nIMF);
                        plv_trialdata=zeros(Nchannels,1,timesol); 
                        %plph_trialdata=zeros(Nchannels,1,timesol);
                        curr_seed=seed(s);
                        if modu_type==2 && conn_type==1
                               parfor icon=1:Nchannels  
                                   %[inich,endch]=conn_idid(icon,:);
                                   %endch=seed;
                                   %inich=icon;
                                   %curr_nIMF = pcurr_nIMF(seed);
                                   %[plv,~]=nspplotf3d_tres3_conn(2*pfm_trialdata{icon}(:,1),2*pph_trialdata{icon}(:,1),seedfm2_trialdata,seedph2_trialdata,timesol,half_winsize,S);
                                   [plv,~]=nspplotf3d_tres3_connLinear(2*pph_trialdata{icon}(:,1),2*pph_trialdata{curr_seed}(:,1),timesol,half_winsize,S);
                                   plv_trialdata(icon,1,:)=shiftdim(plv,-1);
                                   %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                               end
                         elseif modu_type==1 && conn_type==1
                               parfor icon=1:Nchannels  
                                   %[inich,endch]=conn_idid(icon,:);
                                   %endch=seed;
                                   %inich=icon;
                                   %curr_nIMF = pcurr_nIMF(seed);
                                   %[plv,~]=nspplotf3d_tres3_conn(pfm_trialdata{icon}(:,1),pph_trialdata{icon}(:,1),seedfm2_trialdata,seedph2_trialdata,timesol,half_winsize,S);
                                   [plv,~]=nspplotf3d_tres3_connLinear(pph_trialdata{icon}(:,1),pph_trialdata{curr_seed}(:,1),timesol,half_winsize,S);
                                   plv_trialdata(icon,1,:)=shiftdim(plv,-1);
                                   %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                               end
                         elseif modu_type==2 && conn_type==2
                               parfor icon=1:Nchannels  
                                   %[inich,endch]=conn_idid(icon,:);
                                   %endch=seed;
                                   %inich=icon;
                                   %curr_nIMF = pcurr_nIMF(seed);
                                   %[~,plph]=nspplotf3d_tres3_conn(2*pfm_trialdata{icon}(:,1),2*pph_trialdata{icon}(:,1),seedfm2_trialdata,seedph2_trialdata,timesol,half_winsize,S);
                                   [~,plph]=nspplotf3d_tres3_connLinear(2*pph_trialdata{icon}(:,1),2*pph_trialdata{curr_seed}(:,1),timesol,half_winsize,S);
                                   plv_trialdata(icon,1,:)=shiftdim(plph,-1);
                                   %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                               end
                         elseif modu_type==1 && conn_type==2
                               parfor icon=1:Nchannels  
                                   %[inich,endch]=conn_idid(icon,:);
                                   %endch=seed;
                                   %inich=icon;
                                   %curr_nIMF = pcurr_nIMF(seed);
                                   %[~,plph]=nspplotf3d_tres3_conn(pfm_trialdata{icon}(:,1),pph_trialdata{icon}(:,1),seedfm2_trialdata,seedph2_trialdata,timesol,half_winsize,S);
                                   [~,plph]=nspplotf3d_tres3_connLinear(pph_trialdata{icon}(:,1),pph_trialdata{curr_seed}(:,1),timesol,half_winsize,S);
                                   plv_trialdata(icon,1,:)=shiftdim(plph,-1);
                                   %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                               end
                         end
                           Dplv(:,:,:,k)=plv_trialdata;
                           %Dplph(:,:,:,k)=plph_trialdata;
                    end % end k=w
                else
                    for k = w
                        fm_trialdata=D_fm(:,sender,timeindf,k);
                        fm_trialdata=permute(fm_trialdata,[3,1,2]);
                        pfm_trialdata=mat2cell(fm_trialdata,size(fm_trialdata,1),ones(1,size(fm_trialdata,2)));
                        fm2_trialdata=D_fm2(seed(s),:,timeind,k);
                        fm2_trialdata=reshape(fm2_trialdata,size(fm2_trialdata,1),max_IMF,max_imf,[]);
                        fm2_trialdata=permute(fm2_trialdata,[4,2,1,3]);
                        fm2_receiver_trialdata=fm2_trialdata(:,:,:,receiver(irv));
                        pfm2_trialdata=mat2cell(fm2_receiver_trialdata,size(fm2_receiver_trialdata,1),size(fm2_receiver_trialdata,2));
                        ph_trialdata=D_ph(:,sender,timeindf,k);
                        ph_trialdata=permute(ph_trialdata,[3,1,2]);
                        pph_trialdata=mat2cell(ph_trialdata,size(ph_trialdata,1),ones(1,size(ph_trialdata,2)));
                        ph2_trialdata=D_ph2(seed(s),:,timeind,k);
                        ph2_trialdata=reshape(ph2_trialdata,size(ph2_trialdata,1),max_IMF,max_imf,[]);
                        ph2_trialdata=permute(ph2_trialdata,[4,2,1,3]);
                        ph2_receiver_trialdata=ph2_trialdata(:,:,:,receiver(irv));
                        pph2_trialdata=mat2cell(ph2_receiver_trialdata,size(ph2_receiver_trialdata,1),size(ph2_receiver_trialdata,2));
                        pcurr_nIMF = vec_nIMF(k,:);
                        curr_nIMF = pcurr_nIMF(seed(s));
                        seedfm2_trialdata=pfm2_trialdata{1}(:,1:curr_nIMF);
                        seedph2_trialdata=pph2_trialdata{1}(:,1:curr_nIMF);
                        plv_trialdata=zeros(Nchannels,1,timesol); 
                        %plph_trialdata=zeros(Nchannels,1,timesol);
                        if modu_type==2 && conn_type==1
                               parfor icon=1:Nchannels  
                                   %[inich,endch]=conn_idid(icon,:);
                                   %endch=seed;
                                   %inich=icon;
                                   %curr_nIMF = pcurr_nIMF(seed);
                                   [plv,~]=nspplotf3d_tres3_conn(2*pfm_trialdata{icon}(:,1),2*pph_trialdata{icon}(:,1),seedfm2_trialdata,seedph2_trialdata,timesol,half_winsize,S);
                                   plv_trialdata(icon,1,:)=shiftdim(plv,-1);
                                   %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                               end
                         elseif modu_type==1 && conn_type==1
                               parfor icon=1:Nchannels  
                                   %[inich,endch]=conn_idid(icon,:);
                                   %endch=seed;
                                   %inich=icon;
                                   %curr_nIMF = pcurr_nIMF(seed);
                                   [plv,~]=nspplotf3d_tres3_conn(pfm_trialdata{icon}(:,1),pph_trialdata{icon}(:,1),seedfm2_trialdata,seedph2_trialdata,timesol,half_winsize,S);
                                   plv_trialdata(icon,1,:)=shiftdim(plv,-1);
                                   %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                               end
                         elseif modu_type==2 && conn_type==2
                               parfor icon=1:Nchannels  
                                   %[inich,endch]=conn_idid(icon,:);
                                   %endch=seed;
                                   %inich=icon;
                                   %curr_nIMF = pcurr_nIMF(seed);
                                   [~,plph]=nspplotf3d_tres3_conn(2*pfm_trialdata{icon}(:,1),2*pph_trialdata{icon}(:,1),seedfm2_trialdata,seedph2_trialdata,timesol,half_winsize,S);
                                   plv_trialdata(icon,1,:)=shiftdim(plph,-1);
                                   %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                               end
                         elseif modu_type==1 && conn_type==2
                               parfor icon=1:Nchannels  
                                   %[inich,endch]=conn_idid(icon,:);
                                   %endch=seed;
                                   %inich=icon;
                                   %curr_nIMF = pcurr_nIMF(seed);
                                   [~,plph]=nspplotf3d_tres3_conn(pfm_trialdata{icon}(:,1),pph_trialdata{icon}(:,1),seedfm2_trialdata,seedph2_trialdata,timesol,half_winsize,S);
                                   plv_trialdata(icon,1,:)=shiftdim(plph,-1);
                                   %plph_trialdata(icon,1,:)=shiftdim(plph,-1);
                               end
                         end
                           Dplv(:,:,:,k)=plv_trialdata;
                           %Dplph(:,:,:,k)=plph_trialdata;
                    end % end k=w
                end
            end

        %-Save new M/EEG dataset(s)
          Dplv = Dplv.history(mfilename, S);
          save(Dplv);
          %save(Dplph);

        end  % for multiple receivers
   end % for multi seeds
end

