function [T, chan_hood, loc2d, label, S]=perform_paired_tw_Topov(S, f_ranges, chan_hood, cond_cell,  time_roi, max_dist, sigbool, tail)
if size(cond_cell,2)~=2
    errordlg('Paired t-test must select exactly two conditions!','Conditions selected error','modal');
    return;
end
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
if ~isfield(S,'chs')
    S.chs=[];
    chs=[];
else
    chs=S.chs;
end
if isempty(S.timebin)
    timebin=0;
else
    timebin=S.timebin;
end
if ~isfield(S,'logpwr')
    S.logpwr=[];
end
if isempty(S.logpwr)
    logpwr=0;
else
    logpwr=S.logpwr;
end

if ~isfield(S, 'tfce') % only for nonparam
    tfce=0;
else
    tfce=S.tfce;
end

if ~isfield(S, 'nonparam')
    nonparam=0;
else
    nonparam=S.nonparam;
    if nonparam && ~tfce
        if ~isfield(S, 'nonparam_tail')
            nonparam_tail=0;
        else
            nonparam_tail=S.nonparam_tail;
        end
    end
end
if ~isfield(S, 'alpha')
    alpha=0.05;
else
    alpha=S.alpha;
end
if ~isfield(S, 'n_perm')
    n_perm=2000;
else
    n_perm=S.n_perm;
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
    S.subconn_chs=subconn_chs;
end
conn=1;
if isempty(chan_hood)
    mychan_hood=0;
else
    mychan_hood=1;
end
S.Din=D;
if size(cond_cell,1)>1
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
if conn==1 
    if isempty(chs)
        channels=meegchannels(D,{'EEG','MEG','LFP'});
    else
        channels=subconn_chs;
    end
elseif conn==2
    if isempty(chs)
        channels=meegchannels(D,{'EEG','MEG','LFP'});
    else
        channels=chs;
    end
else
    channels=meegchannels(D,{'EEG','MEG','LFP'});
end

structD=struct(D);
Freqs=frequencies(D);
    
%for h=1:length(bandbool)
h=1;
idF=find(Freqs>=f_ranges{h}(1) & Freqs<=f_ranges{h}(2));
Ddata=zeros(nchannels(D), 1, nsamples(D), ntrials(D));
for c=1:D.ntrials
    x=spm_squeeze(D(:,idF,:,c), 4);
    x=mean(x,2);
    %D(:,1,:,c)= x;
    Ddata(:,1,:,c)= x;
end
if size(cond_cell,1)>1
   switch choice
        case 'first'
            idx1=indtrial(D,cond_cell{1,1});
            idx2=indtrial(D,cond_cell{1,2});
            tmpdata1=Ddata(channels,1,slice_roi(1):slice_roi(2),idx1);
            tmpdata2=Ddata(channels,1,slice_roi(1):slice_roi(2),idx2);
            tmpdata1=permute(tmpdata1,[2,1,3,4]);
            tmpdata1=shiftdim(tmpdata1);
            tmpdata2=permute(tmpdata2,[2,1,3,4]);
            tmpdata2=shiftdim(tmpdata2);
        case 'concatenation'
            idx1=indtrial(D,cond_cell{1,1});
            idx2=indtrial(D,cond_cell{1,2});
%                     repl_1=getRepl(structD, idx1);
%                     repl_2=getRepl(structD, idx2);
            if weighted
                repl_1=getRepl(structD, idx1);
                repl_2=getRepl(structD, idx2);
            else
                repl_1=ones(1,length(idx1));
                repl_2=ones(1,length(idx2));
            end

            repl_1_4d(1,1,:)=repl_1;
            repl_2_4d(1,1,:)=repl_2;
            tmp1=Ddata(channels,1,slice_roi(1):slice_roi(2),idx1);
            tmp2=Ddata(channels,1,slice_roi(1):slice_roi(2),idx2);
            tmp1=permute(tmp1,[2,1,3,4]);
            tmp1=shiftdim(tmp1);
            tmp2=permute(tmp2,[2,1,3,4]);
            tmp2=shiftdim(tmp2);
            tmpdata1=tmp1.*repmat(repl_1_4d,[length(channels),length(slice_roi(1):slice_roi(2)),1]);
            tmpdata2=tmp2.*repmat(repl_2_4d,[length(channels),length(slice_roi(1):slice_roi(2)),1]);
            c1=repl_1;
            c2=repl_2;
            for k=2:size(cond_cell,1)
                idxh1=indtrial(D,cond_cell{k,1});
%                 idx1=[idx1 idxh1];
                idxh2=indtrial(D,cond_cell{k,2});
%                 idx2=[idx2 idxh2];
%                         repl_1=getRepl(structD, idxh1);
%                         repl_2=getRepl(structD, idxh2);
                if weighted
                    repl_1=getRepl(structD, idxh1);
                    repl_2=getRepl(structD, idxh2);
                else
                    repl_1=ones(1,length(idxh1));
                    repl_2=ones(1,length(idxh2));
                end

                repl_1_4d(1,1,:)=repl_1;
                repl_2_4d(1,1,:)=repl_2;
                tmp1=Ddata(channels,1,slice_roi(1):slice_roi(2),idxh1);
                tmp2=Ddata(channels,1,slice_roi(1):slice_roi(2),idxh2);
                tmp1=permute(tmp1,[2,1,3,4]);
                tmp1=shiftdim(tmp1);
                tmp2=permute(tmp2,[2,1,3,4]);
                tmp2=shiftdim(tmp2);
                tmpdata1=tmpdata1+tmp1.*repmat(repl_1_4d,[length(channels),length(slice_roi(1):slice_roi(2)),1]);
                tmpdata2=tmpdata2+tmp2.*repmat(repl_2_4d,[length(channels),length(slice_roi(1):slice_roi(2)),1]);
                c1=c1+repl_1;
                c2=c2+repl_2;
            end
            c1_4d(1,1,:)=c1;
            c2_4d(1,1,:)=c2;
            tmpdata1=tmpdata1./repmat(c1_4d,[length(channels),length(slice_roi(1):slice_roi(2)),1]);
            tmpdata2=tmpdata2./repmat(c2_4d,[length(channels),length(slice_roi(1):slice_roi(2)),1]);

        case 'no test'
            return
   end
elseif size(cond_cell,1)==1
    idx1=indtrial(D,cond_cell{1,1});
    idx2=indtrial(D,cond_cell{1,2});
    tmpdata1=Ddata(channels,1,slice_roi(1):slice_roi(2),idx1);
    tmpdata2=Ddata(channels,1,slice_roi(1):slice_roi(2),idx2);
    tmpdata1=permute(tmpdata1,[2,1,3,4]);
    tmpdata1=shiftdim(tmpdata1);
    tmpdata2=permute(tmpdata2,[2,1,3,4]);
    tmpdata2=shiftdim(tmpdata2);
else
    return
end
if timebin==0
   if ~ logpwr 
       avgdata1=tmpdata1;
       avgdata2=tmpdata2;
   else
       avgdata1=tmpdata1;
       avgdata2=tmpdata2;
       avgdata1=log10(avgdata1+0.00001);
       avgdata2=log10(avgdata2+0.00001);
   end
else
    bno=abs(round(timebin));
    if bno>1
        avgdata1=zeros(size(tmpdata1,1),bno,size(tmpdata1,3));
        avgdata2=zeros(size(tmpdata2,1),bno,size(tmpdata2,3));
        b_sep=linspace(1,size(tmpdata1,2),bno+1);
        b_sep_ceil=ceil(b_sep);
        b_sep_floor=floor(b_sep);
        b_sep_ids=[b_sep_ceil(1:bno);b_sep_floor(2:(bno+1))]';
        for ib=1:bno
            avgdata1(:,ib,:)=mean(tmpdata1(:,b_sep_ids(ib,1):b_sep_ids(ib,2),:),2);
            avgdata2(:,ib,:)=mean(tmpdata2(:,b_sep_ids(ib,1):b_sep_ids(ib,2),:),2);
        end
        if logpwr
        avgdata1=log10(avgdata1+0.00001);
        avgdata2=log10(avgdata2+0.00001);
        end
    elseif bno==1
        avgdata1=zeros(size(tmpdata1,1),1,size(tmpdata1,3));
        avgdata2=zeros(size(tmpdata2,1),1,size(tmpdata2,3));
        avgdata1(:,1,:)=mean(tmpdata1(:,:,:),2);
        avgdata2(:,1,:)=mean(tmpdata2(:,:,:),2);
    end
end
if ~nonparam
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
end
if nonparam
    if ~mychan_hood
        chan_hood=spatial_neighbors_spm(chanlocs,max_dist);
        if conn==2 && ~isempty(chs)
            chan_hood=chan_hood(chs,chs);
        end
    else
        if conn==1 && ~isempty(chs)
            chan_hood=chan_hood(subconn_chs,subconn_chs);
        end
    end
   if ~tfce 
       [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm1x(avgdata1-avgdata2,chan_hood,n_perm,0.05,nonparam_tail,alpha);
       T(1,h).t=t_orig;
       T(1,h).pval=pval;
       T(1,h).tail=nonparam_tail;
   else
        maxN_nbrs=max(sum(chan_hood,2));
        tfce_nbrs=zeros(size(chan_hood,1),maxN_nbrs);
        for gch=1:size(chan_hood,1)
            tmp_nbrs=find(chan_hood(gch,:));
            if ~isempty(tmp_nbrs)
              tfce_nbrs(gch,1:length(tmp_nbrs))=tmp_nbrs';
            end
        end

        Results = ept_TFCE(permute(avgdata1,[3,1,2]), permute(avgdata2,[3,1,2]), [], 'nPerm',n_perm,'type','d','chn',tfce_nbrs);
        T(1,h).t=Results.Obs;
        T(1,h).pval=Results.P_Values;
        T(1,h).tail=0;
   end
end
T(1,h).m_data={mean(avgdata1,ndims(avgdata1)) , mean(avgdata2,ndims(avgdata2))};
T(1,h).data={avgdata1 , avgdata2};
function repl=getRepl(structD, idx)
    repl=zeros(1,length(idx));
    for k=1:length(repl)
        repl(k)=structD.trials(1,idx(k)).repl;
    end


