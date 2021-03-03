function [T,chan_hood, loc2d, label]=perform_onesample_t_Topo(S, f_ranges, chan_hood, cond_cell,  time_roi, max_dist, sigbool, tail)
if nargin<8
    tail=0;
end
if ~isfield(S, 'weighted') || isempty(S.weighted) 
    weighted=1;
else
    weighted=S.weighted;
end
if ~isfield(S,'timebin')
    S.timebin=[];
end
if isempty(S.timebin)
    timebin=0;
else
    timebin=S.timebin;
end
if ~isfield(S,'chs')
    S.chs=[];
    chs=[];
else
    chs=S.chs;
end
if ~isfield(S, 'alpha')
    alpha=0.05;
else
    alpha=S.alpha;
end
D  = spm_eeg_load(S.D);
if isempty(chs) % all chs are included
   conn_idid=D.conn_idid;
   loc2d=conn_idid(:,[2, 1]);
   label=chanlabels(D,1:length(loc2d))';
else
    A=ismember(D.conn_idid(:,1),chs);
    B=ismember(D.conn_idid(:,2),chs);
    subconn_chs=find(A & B);
    subconn_idid=D.conn_idid(subconn_chs,:);
    [~,locb]=ismember(subconn_idid,chs);
    loc2d=locb(:,[2, 1]);
    label=chanlabels(D,subconn_chs)';
end
conn=1;
if isempty(chan_hood)
    mychan_hood=0;
else
    mychan_hood=1;
end
S.Din=D;
if length(cond_cell(:))>1
    choice='concatenation';
end
    ep_initime=D.timeonset*1000;
    roi=time_roi;
    abs_roi=[D.time(1) D.time(end)]*1000;
    if roi(1)<abs_roi(1)
       roi(1)=abs_roi(1) ;
    end
    if roi(2)>abs_roi(2)
       roi(2)=abs_roi(2) ;
    end
    samp_interval=1000/D.fsample;
    slice_roi=round((roi-ep_initime)/samp_interval)+1;
%slice_times=ep_initime+samp_interval*((slice_roi(1)-1):(slice_roi(2)-1));
if conn==1 
    if isempty(chs)
        channels=meegchannels(D,{'EEG','MEG','LFP'});
    else
        channels=subconn_chs;
    end
else
    channels=meegchannels(D,{'EEG','MEG','LFP'});
end
structD=struct(D);
Freqs=frequencies(D);
h=1;
 idF=find(Freqs>=f_ranges{h}(1) & Freqs<=f_ranges{h}(2));
 Ddata=zeros(nchannels(D), 1, nsamples(D), ntrials(D));
 for c=1:D.ntrials
    x=spm_squeeze(D(:,idF,:,c), 4);
    x=mean(x,2);
    %D(:,1,:,c)= x;
    Ddata(:,1,:,c)= x;
 end
 if length(cond_cell(:))>1   
       switch choice
            case 'first'
                 idx1=indtrial(D,cond_cell{1});
                 tmpdata1=Ddata(channels,1,slice_roi(1):slice_roi(2),idx1);
                 tmpdata1=permute(tmpdata1,[2,1,3,4]);
                 tmpdata1=shiftdim(tmpdata1);
            case 'concatenation'
                idx1=indtrial(D,cond_cell{1});
                %repl_1=getRepl(structD, idx1);
                if weighted
                    repl_1=getRepl(structD, idx1);
                else
                    repl_1=ones(1,length(idx1));
                end
                repl_1_4d(1,1,:)=repl_1;
                tmp=Ddata(channels,1,slice_roi(1):slice_roi(2),idx1);
                tmp=permute(tmp,[2,1,3,4]);
                tmp=shiftdim(tmp);
                tmpdata1=tmp.*repmat(repl_1_4d,[length(channels),length(slice_roi(1):slice_roi(2)),1]);
                c1=repl_1;
                for k=2:size(cond_cell,1)*size(cond_cell,2)
                    idxh1=indtrial(D,cond_cell{k});
                    %repl_1=getRepl(structD, idxh1);
                    if weighted
                        repl_1=getRepl(structD, idxh1);
                    else
                        repl_1=ones(1,length(idxh1));
                    end

                    repl_1_4d(1,1,:)=repl_1;
                    tmp=Ddata(channels,1,slice_roi(1):slice_roi(2),idxh1);
                    tmp=permute(tmp,[2,1,3,4]);
                    tmp=shiftdim(tmp);
                    tmpdata1=tmpdata1+tmp.*repmat(repl_1_4d,[length(channels),length(slice_roi(1):slice_roi(2)),1]);
                    c1=c1+repl_1;
                end
                c1_4d(1,1,:)=c1;
                tmpdata1=tmpdata1./repmat(c1_4d,[length(channels),length(slice_roi(1):slice_roi(2)),1]);

            case 'no test'
                return
       end
 elseif length(cond_cell(:))==1
        idx1=indtrial(D,cond_cell{1});
        tmpdata1=Ddata(channels,1,slice_roi(1):slice_roi(2),idx1);
        tmpdata1=permute(tmpdata1,[2,1,3,4]);
        tmpdata1=shiftdim(tmpdata1);
 else
        return
 end
 if timebin==0   
     avgdata1=tmpdata1;
     avgdata2=zeros(size(avgdata1));
 else
     bno=abs(round(timebin));
     if bno>1   
        avgdata1=zeros(size(tmpdata1,1),bno,size(tmpdata1,3));
        %avgdata2=zeros(size(tmpdata2,1),bno,size(tmpdata2,3));
        b_sep=linspace(1,size(tmpdata1,2),bno+1);
        b_sep_ceil=ceil(b_sep);
        b_sep_floor=floor(b_sep);
        b_sep_ids=[b_sep_ceil(1:bno);b_sep_floor(2:(bno+1))]';
        for ib=1:bno
            avgdata1(:,ib,:)=mean(tmpdata1(:,b_sep_ids(ib,1):b_sep_ids(ib,2),:),2);
            %avgdata2(:,ib,:)=mean(tmpdata2(:,b_sep_ids(ib,1):b_sep_ids(ib,2),:),2);
        end
        avgdata2=zeros(size(avgdata1));
     elseif bno==1
        avgdata1=zeros(size(tmpdata1,1),1,size(tmpdata1,3));
        %avgdata2=zeros(size(tmpdata2,1),bno,size(tmpdata2,3));
            avgdata1(:,1,:)=mean(tmpdata1(:,:,:),2);
            %avgdata2(:,ib,:)=mean(tmpdata2(:,b_sep_ids(ib,1):b_sep_ids(ib,2),:),2);
        avgdata2=zeros(size(avgdata1));
     end
 end
 if S.stats
    if sum(sigbool(1,1:5))>0
        [t df] = ttest_cell({ avgdata1 avgdata2 });
        T(1,h).t=t;
        T(1,h).df=df;
        if sigbool(1)==1
            if tail==1
               t_p=tinv(0.95,df);
            elseif tail==0
                t_p=tinv(0.975,df);
            end
            s=abs(t)>t_p;
        end
        if sigbool(2)==1
            if tail==1
               t_p=tinv(0.99,df);
            elseif tail==0
                t_p=tinv(0.995,df);
            end
            s=abs(t)>t_p;
        end
        if sigbool(3)==1
            if tail==1
               t_p=tinv(0.999,df);
            elseif tail==0
                t_p=tinv(0.9995,df);
            end
            s=abs(t)>t_p;
        end
        if sigbool(4)==1
            if tail==1
               P=1-tcdf(abs(t),df);
            elseif tail==0
               P=2*(1-tcdf(abs(t),df));
            end
            [p_fdr,P_masked]=fdr(P,0.05);
            s=P_masked;
        end
        if sigbool(5)==1
            if tail==1
               P=1-tcdf(abs(t),df);
            elseif tail==0
               P=2*(1-tcdf(abs(t),df));
            end
            [p_fdr,P_masked]=fdr(P,0.01);
            s=P_masked;
        end
        T(1,h).s=s;
        T(1,h).tail=tail;
    end
    if sum(sigbool(1,8:10))>0
        if ~mychan_hood
            chan_hood=spatial_neighbors_spm(chanlocs,max_dist);
            if conn==2 && ~isempty(chs)
                chan_hood=chan_hood(chs,chs);
            end
        end
    end
    if sigbool(8)==1
       %chan_hood=spatial_neighbors_spm(chanlocs,max_dist);
       [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm1(avgdata1,chan_hood,2000,0.05,0,alpha);
       T(1,h).t=t_orig;
       T(1,h).pval=pval;
       T(1,h).tail=0;
    end
    if sigbool(9)==1
       %chan_hood=spatial_neighbors_spm(chanlocs,max_dist); 
       [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm1(avgdata1,chan_hood,2000,0.05,1,alpha);
       T(1,h).t=t_orig;
       T(1,h).pval=pval;
       T(1,h).tail=1;
    end
    if sigbool(10)==1
       %chan_hood=spatial_neighbors_spm(chanlocs,max_dist); 
       [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm1(avgdata1,chan_hood,2000,0.05,-1,alpha);
       T(1,h).t=t_orig;
       T(1,h).pval=pval;
       T(1,h).tail=-1;
    end
    T(1,h).m_data=shiftdim(mean(avgdata1,ndims(avgdata1)));
    T(1,h).t=T(1,h).m_data;
 else
    T(1,h).m_data=shiftdim(mean(avgdata1,ndims(avgdata1)));
    T(1,h).t=T(1,h).m_data;
    T(1,h).pval=ones(size(T(1,h).t));
    T(1,h).s=zeros(size(T(1,h).t));
    chan_hood=[];
 end

  %end
%end
% plot_bands_tmap_1D_fdr( T, slice_times,sigbool,cond_cell);
%view_topography_T
function repl=getRepl(structD, idx)
    repl=zeros(1,length(idx));
    for k=1:length(repl)
        repl(k)=structD.trials(1,idx(k)).repl;
    end


