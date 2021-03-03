function [T,P,cfg,df,testdata]=perform_paired_tw2n(D,cond_cell,channels,time_roi,idFreqs,S)
% sigbool must be a 3x1 boolean vector corresponding to the boundary of
% p value equal to [0.05, 0.01, 0.001]
if size(cond_cell,2)~=2
    errordlg('Paired t-test must select exactly two conditions!','Conditions selected error','modal');
    return;
end
if iscell(channels) && size(channels,2)==2
    n2pc=1;
else
    n2pc=0;
end

if ~isfield(S, 'weighted') || isempty(S.weighted) 
    weighted=1;
else
    weighted=S.weighted;
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
        if ~isfield(S, 'nonparam_dcmax')
            nonparam_dcmax=0;
        else
            nonparam_dcmax=S.nonparam_dcmax;
        end
    end
end
if ~isfield(S, 'alpha')
    alpha=0.05;
else
    alpha=S.alpha;
end
if isempty(alpha)
    alpha=0.05;
end
if ~isfield(S, 'n_perm')
    n_perm=2000;
else
    n_perm=S.n_perm;
end

if ~isempty(strfind(D.fnamedat,'omegaf'))
    fF=1;
else
    fF=0;
end
ph_time=~isempty(strfind(D.fnamedat,'hh1fp'))||~isempty(strfind(D.fnamedat,'hh2fp'))||~isempty(strfind(D.fnamedat,'hh1p'))||~isempty(strfind(D.fnamedat,'hh2p'));

if ~isempty(strfind(D.fnamedat,'pdf'))||~isempty(strfind(D.fnamedat,'mse'))
    pdf=1;
else
    pdf=0;
end
if ~isempty(strfind(D.fnamedat,'ierp'))
    ierp=1;
else
    ierp=0;
end
if ~isfield(D, 'dyadic')
    dyadic=0;
else
    dyadic=D.dyadic;
end
if ~isfield(S, 'collapse')
    collapse=0;
else
    collapse=S.collapse;
end
if collapse && idFreqs(1)==2
   dc_inc=1;
   idFreqs=[1;idFreqs];
else
   dc_inc=0; 
end
if fF
    ep_initime=D.timeonset;
    if dyadic
       roi=log2(time_roi); 
    else
       roi=time_roi;
    end
    abs_roi=[D.time(1) D.time(end)];
    if roi(1)<abs_roi(1)
       roi(1)=abs_roi(1) ;
    end
    if roi(2)>abs_roi(2)
       roi(2)=abs_roi(2) ;
    end
    samp_interval=1/D.fsample;
    slice_roi=round((roi-ep_initime)/samp_interval)+1;
elseif pdf|| ph_time
    ep_initime=D.timeonset;
    roi=time_roi;
    abs_roi=[D.time(1) D.time(end)];
    if roi(1)<abs_roi(1)
       roi(1)=abs_roi(1) ;
    end
    if roi(2)>abs_roi(2)
       roi(2)=abs_roi(2) ;
    end
    samp_interval=1/D.fsample;
    slice_roi=round((roi-ep_initime)/samp_interval)+1;
else
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
end
structD=struct(D);
if n2pc 
else
        if size(cond_cell,1)>1
            %choice = questdlg('for paired t-test, you can enter only two conditions. Test first row or concatenating to one row','t-test limit','first','concatenation','no test','concatenation');
           choice='concatenation'; 
           switch choice
                case 'first'
                    idx1=indtrial(D,cond_cell{1,1},'GOOD');
                    idx2=indtrial(D,cond_cell{1,2},'GOOD');
                    tmpdata1=D(channels,idFreqs,slice_roi(1):slice_roi(2),idx1);
                    tmpdata2=D(channels,idFreqs,slice_roi(1):slice_roi(2),idx2);
                case 'concatenation'
                    idx1=indtrial(D,cond_cell{1,1},'GOOD');
                    idx2=indtrial(D,cond_cell{1,2},'GOOD'); 
                    if weighted
                        repl_1=getRepl(structD, idx1);
                        repl_2=getRepl(structD, idx2);
                    else
                        repl_1=ones(1,length(idx1));
                        repl_2=ones(1,length(idx2));
                    end
                    repl_1_4d(1,1,1,:)=repl_1;
                    repl_2_4d(1,1,1,:)=repl_2;
                    tmpdata1=D(channels,idFreqs,slice_roi(1):slice_roi(2),idx1).*repmat(repl_1_4d,[length(channels),length(idFreqs),length(slice_roi(1):slice_roi(2)),1]);
                    tmpdata2=D(channels,idFreqs,slice_roi(1):slice_roi(2),idx2).*repmat(repl_2_4d,[length(channels),length(idFreqs),length(slice_roi(1):slice_roi(2)),1]);
                    c1=repl_1;
                    c2=repl_2;
                    for h=2:size(cond_cell,1)
                        idxh1=indtrial(D,cond_cell{h,1},'GOOD');
                        idxh2=indtrial(D,cond_cell{h,2},'GOOD');
                        if weighted
                            repl_1=getRepl(structD, idxh1);
                            repl_2=getRepl(structD, idxh2);
                        else
                            repl_1=ones(1,length(idxh1));
                            repl_2=ones(1,length(idxh2));
                        end
                        repl_1_4d(1,1,1,:)=repl_1;
                        repl_2_4d(1,1,1,:)=repl_2;
                        tmpdata1=tmpdata1+D(channels,idFreqs,slice_roi(1):slice_roi(2),idxh1).*repmat(repl_1_4d,[length(channels),length(idFreqs),length(slice_roi(1):slice_roi(2)),1]);
                        tmpdata2=tmpdata2+D(channels,idFreqs,slice_roi(1):slice_roi(2),idxh2).*repmat(repl_2_4d,[length(channels),length(idFreqs),length(slice_roi(1):slice_roi(2)),1]);
                        c1=c1+repl_1;
                        c2=c2+repl_2;
                    end
                    c1_4d(1,1,1,:)=c1;
                    c2_4d(1,1,1,:)=c2;
                    tmpdata1=tmpdata1./repmat(c1_4d,[length(channels),length(idFreqs),length(slice_roi(1):slice_roi(2)),1]);
                    tmpdata2=tmpdata2./repmat(c2_4d,[length(channels),length(idFreqs),length(slice_roi(1):slice_roi(2)),1]);

                case 'no test'
                    return
           end
        elseif size(cond_cell,1)==1
            idx1=indtrial(D,cond_cell{1,1},'GOOD');
            idx2=indtrial(D,cond_cell{1,2},'GOOD');
            tmpdata1=D(channels,idFreqs,slice_roi(1):slice_roi(2),idx1);
            tmpdata2=D(channels,idFreqs,slice_roi(1):slice_roi(2),idx2);
        end
end
if length(idx1)~=length(idx2)
    errordlg('The number of trials in the two conditions in paired t-test must be equal!','Trial selected error','modal');
    return;
end
if length(idx1)==1
    errordlg('The number of trials for paired t-test must be larger than 1!','Trial selected error','modal');
    return;
end

 
if size(tmpdata1,1)==1
    avgdata1=tmpdata1;
    avgdata2=tmpdata2;
else
    avgdata1=mean(tmpdata1);
    avgdata2=mean(tmpdata2);
end
avgdata1=shiftdim(avgdata1);
avgdata2=shiftdim(avgdata2);
if ph_time
   avgdata1=repmat(avgdata1(:,2:(end-1),:),[1 2 1]); 
   avgdata2=repmat(avgdata2(:,2:(end-1),:),[1 2 1]); 
end

if fF
    xmin=D.time(slice_roi(1));
    xmax=D.time(slice_roi(2));
elseif ph_time
%     xmin=roi(1);
%     xmax=roi(2)+2*pi;
    xmin=-pi;
    xmax=3*pi;
else
%     xmin=roi(1);
%     xmax=roi(2);
    xmin=1000*D.time(slice_roi(1));
    xmax=1000*D.time(slice_roi(2));
end
Freqs=frequencies(D);
if dc_inc
    ymin=Freqs(idFreqs(2));
else
    ymin=Freqs(idFreqs(1));
end
ymax=Freqs(idFreqs(end));
if fF
    yy=repmat(Freqs(idFreqs),1,size(avgdata1,2));
    xx=repmat(D.time(slice_roi(1):slice_roi(2)),size(avgdata1,1),1); 
    delbool=(yy>(xx-0.25));
    delbool_3d=repmat(delbool,1,1,size(avgdata1,3));
    avgdata1(delbool_3d)=0;
    avgdata2(delbool_3d)=0;
end
testdata=avgdata1-avgdata2;
if ~ nonparam
    [T df] = ttest_cell({ avgdata1 avgdata2 });
    P=1-tcdf(abs(T),df);
    test1=isnan(T);
    i_nan = find(test1);
    T(i_nan)=0;
    test2=isnan(P);
    i_nan = find(test2);
    P(i_nan)=1;
else
    if ~tfce
        fsize=size(testdata,1);
        chan_hood=zeros(fsize,fsize);
            if ierp
                cIMFs=D.cIMFs(idFreqs,:);
                for j=1:fsize
                    set_j=cIMFs(j,1):cIMFs(j,2);
                    for k=1:fsize
                        set_k=cIMFs(k,1):cIMFs(k,2);
                        tmp_xor=setxor(set_j,set_k);
                        if length(tmp_xor)<=1
                           chan_hood(j,k)=1;
                        end
                    end
                end
            else
                for j=1:fsize
                    if j==1
                        chan_hood(1,1)=1;
                        chan_hood(1,2)=1;
                    elseif j==fsize
                        chan_hood(fsize,fsize-1)=1;
                        chan_hood(fsize,fsize)=1;
                    else
                        chan_hood(j,j-1)=1;
                        chan_hood(j,j)=1;
                        chan_hood(j,j+1)=1;
                    end
                end
            end
        try
           cluster_roi=S.cluster_roi;
        catch
            cluster_roi=[];
        end
        if isempty(cluster_roi)
            if dc_inc && nonparam_dcmax
                [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm1dcam(testdata,chan_hood,n_perm,0.05,nonparam_tail,alpha);
            else
                [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm1x(testdata,chan_hood,n_perm,0.05,nonparam_tail,alpha);
            end
            [T df] = ttest_cell({ avgdata1 avgdata2 });
            %T=t_orig;
            P=pval;
           test1=isnan(T);
           i_nan = find(test1);
           T(i_nan)=0;
        else
                if cluster_roi(1)< xmin %roi(1)
                   cluster_roi(1)=xmin; %roi(1) ;
                end
                if cluster_roi(2)> xmax %roi(2)
                   cluster_roi(2)=xmax; % roi(2) ;
                end
                ini_id=ceil((cluster_roi(1)-xmin)/samp_interval)+1;
                end_id=floor((cluster_roi(2)-xmin)/samp_interval)+1;
            cluster_id=ini_id:end_id;
            testdata_roi=testdata(:,cluster_id,:);
            [T df] = ttest_cell({ avgdata1 avgdata2 });
            P=ones(size(T));
            if dc_inc && nonparam_dcmax
                [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm1dcam(testdata_roi,chan_hood,n_perm,0.05,nonparam_tail,alpha);
            else
                [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm1(testdata_roi,chan_hood,n_perm,0.05,nonparam_tail,alpha);
            end
            %T=t_orig;
            test1=isnan(T);
            i_nan = find(test1);
            T(i_nan)=0;
            P(:,cluster_id)=pval;
        end
    else
        try
           cluster_roi=S.cluster_roi;
        catch
            cluster_roi=[];
        end
        fsize=size(testdata,1);
        tfce_nbrs=zeros(fsize,fsize);
        if ierp
            cIMFs=D.cIMFs(idFreqs,:);
            for j=1:fsize
                set_j=cIMFs(j,1):cIMFs(j,2);
                count_j=0;
                for k=1:fsize
                    set_k=cIMFs(k,1):cIMFs(k,2);
                    tmp_xor=setxor(set_j,set_k);
                    if length(tmp_xor)<=1
                        count_j=count_j+1;
                        tfce_nbrs(j,count_j)=k;
                    end
                end
            end
        else
            for j=1:fsize
                if j==1
                    tfce_nbrs(1,1)=1;
                    tfce_nbrs(1,2)=2;
                elseif j==fsize
                    tfce_nbrs(fsize,1)=fsize-1;
                    tfce_nbrs(fsize,2)=fsize;
                else
                    tfce_nbrs(j,1)=j-1;
                    tfce_nbrs(j,2)=j;
                    tfce_nbrs(j,3)=j+1;
                end
            end
        end
        if isempty(cluster_roi)
%             avgdata1_roi=shiftdim(avgdata1,-1);
%             avgdata2_roi=shiftdim(avgdata2,-1);
%             avgdata1_imx=permute(avgdata1_roi, [4,1,2,3]);
%             avgdata2_imx=permute(avgdata2_roi, [4,1,2,3]);
%             avgdata1_imx=repmat(avgdata1_imx,[1,2,1,1]);
%             avgdata2_imx=repmat(avgdata2_imx,[1,2,1,1]);
%             tfce_nbrs=[1,2;1,2]; % two psudo channels
%             Results = ept_TFCE(avgdata1_imx, avgdata2_imx,[] , 'nPerm',n_perm,'type','d','chn',tfce_nbrs);
            avgdata1_roi=avgdata1;
            avgdata2_roi=avgdata2;
            avgdata1_imx=permute(avgdata1_roi, [3,1,2]);
            avgdata2_imx=permute(avgdata2_roi, [3,1,2]);
            Results = ept_TFCE(avgdata1_imx, avgdata2_imx,[] , 'nPerm',n_perm,'type','d','chn',tfce_nbrs);

            [~, df] = ttest_cell({ avgdata1 avgdata2 });
            T=Results.Obs(:,:);
            test1=isnan(T);
            i_nan = find(test1);
            T(i_nan)=0;
            P=Results.P_Values(:,:);
        else
                if cluster_roi(1)< xmin %roi(1)
                   cluster_roi(1)=xmin; %roi(1) ;
                end
                if cluster_roi(2)> xmax %roi(2)
                   cluster_roi(2)=xmax; % roi(2) ;
                end
                ini_id=ceil((cluster_roi(1)-xmin)/samp_interval)+1;
                end_id=floor((cluster_roi(2)-xmin)/samp_interval)+1;
            cluster_id=ini_id:end_id;
%             avgdata1_roi=shiftdim(avgdata1(:,cluster_id,:),-1);
%             avgdata2_roi=shiftdim(avgdata2(:,cluster_id,:),-1);
%             avgdata1_imx=permute(avgdata1_roi, [4,1,2,3]);
%             avgdata2_imx=permute(avgdata2_roi, [4,1,2,3]);
%             avgdata1_imx=repmat(avgdata1_imx,[1,2,1,1]);
%             avgdata2_imx=repmat(avgdata2_imx,[1,2,1,1]);
%             tfce_nbrs=[1,2;1,2]; % two psudo channels
            avgdata1_roi=avgdata1(:,cluster_id,:);
            avgdata2_roi=avgdata2(:,cluster_id,:);
            avgdata1_imx=permute(avgdata1_roi, [3,1,2]);
            avgdata2_imx=permute(avgdata2_roi, [3,1,2]);
            Results = ept_TFCE(avgdata1_imx, avgdata2_imx,[] , 'nPerm',n_perm,'type','d','chn',tfce_nbrs);
            [T df] = ttest_cell({ avgdata1 avgdata2 });
            P=ones(size(T));
            %T=t_orig;
            %T=Results.Obs';
            test1=isnan(T);
            i_nan = find(test1);
            T(i_nan)=0;
            P(:,cluster_id)=Results.P_Values(:,:);
        end
    end
end

t_max = tinv(0.995,df);
t_min=-1*t_max;
clim=[t_min,t_max];

cfg.clim=clim;
cfg.xmin=xmin;
cfg.xmax=xmax;
cfg.ymin=ymin;
cfg.ymax=ymax;
cfg.fF=fF;
cfg.dyadic=dyadic;

%colorbar('FontSize',14);
function repl=getRepl(structD, idx)
    repl=zeros(1,length(idx));
    for k=1:length(repl)
        repl(k)=structD.trials(1,idx(k)).repl;
    end
