function [T,P, cfg]=perform_onesample_correlation2(avgdata1, behavior, D, time_roi, idFreqs, S)
% T: correlation, P: p value of correlation
% sigbool must be a 3x1 boolean vector corresponding to the boundary of
% p value equal to [0.05, 0.01, 0.001]


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
        if ~isfield(S, 'nonparam_clusmod')
            nonparam_clusmod=0;
        else
            nonparam_clusmod=S.nonparam_clusmod;
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
elseif pdf || ph_time
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
if strcmp(D.transformtype, 'TF') && isfield(D,'wtc')
    idFreqs=1:D.nfrequencies;
end

% fq_roi=[1,30];
% clim=[-100, 100];
% slice_roi=round((roi-ep_initime)/samp_interval)+1;


% tmpdata1=D(channels,idFreqs,slice_roi(1):slice_roi(2),idx1);
% % tmpdata2=D(channels,fq_roi(1):fq_roi(2),slice_roi(1):slice_roi(2),idx2); 
% if size(tmpdata1,1)==1
%     avgdata1=tmpdata1;
% %     avgdata2=tmpdata2;
% else
%     avgdata1=mean(tmpdata1);
% %     avgdata2=mean(tmpdata2);
% end
% avgdata1=squeeze(avgdata1);
% avgdata2=zeros(size(avgdata1));
% avgdata2=squeeze(avgdata2);
% [T df] = ttest_cell({ avgdata1 avgdata2 });
% t_max = tinv(0.999,df);
% t_min=-1*t_max;
% clim=[t_min,t_max];
% clim=[-100,100];
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
    delboolp=(yy>(xx-1));
    %delboolp_3d=repmat(delboolp,1,1,size(avgdata1,3));
    %P(delboolp_3d)=1;

end
if ~ nonparam
    T=zeros(size(avgdata1,1),size(avgdata1,2));
    P=ones(size(avgdata1,1),size(avgdata1,2));
    for j=1:size(avgdata1,1)
        tmpdata=squeeze(avgdata1(j,:,:));
        tmpdata=permute(tmpdata,[2,1]);
        [rho,pval] = corr(behavior,tmpdata);
        T(j,:)=rho;
        P(j,:)=pval;
    end
    test1=isnan(T);
    i_nan = find(test1);
    T(i_nan)=0;
    test2=isnan(P);
    i_nan = find(test2);
    P(i_nan)=1;
    if fF
       P(delboolp)=1; 
    end
else
    if ~tfce
        %testdata=avgdata1-avgdata2;
        fsize=size(avgdata1,1);
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
           behaviorEx=repmat(shiftdim(behavior,-2),[size(avgdata1,1), size(avgdata1,2)]);
           if fF && nonparam_clusmod
               [pval,t_orig, corr_obs, est_alpha, seed_state]=clust_perm_corrmodu(avgdata1,behaviorEx,chan_hood,n_perm,0,alpha,xx,yy);
           elseif dc_inc && nonparam_dcmax
               [pval,t_orig, corr_obs, est_alpha, seed_state]=clust_perm_corrdcam(avgdata1,behaviorEx,chan_hood,n_perm,0,alpha);
           else
               [pval,t_orig, corr_obs, est_alpha, seed_state]=clust_perm_corrx(avgdata1,behaviorEx,chan_hood,n_perm,0,alpha);
           end
           T=corr_obs;
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
            avgdata1_roi=avgdata1(:,cluster_id,:);
           behaviorEx=repmat(shiftdim(behavior,-2),[size(avgdata1_roi,1), size(avgdata1_roi,2)]);
           if dc_inc && nonparam_dcmax
               [pval,t_orig, corr_obs, est_alpha, seed_state]=clust_perm_corrdcam(avgdata1_roi,behaviorEx,chan_hood,n_perm,0,alpha);
           else
               [pval,t_orig, corr_obs, est_alpha, seed_state]=clust_perm_corrx(avgdata1_roi,behaviorEx,chan_hood,n_perm,0,alpha);
           end
           T=zeros(size(avgdata1,1),size(avgdata1,2));
           P=ones(size(avgdata1,1),size(avgdata1,2));
           for j=1:size(avgdata1,1)
                tmpdata=squeeze(avgdata1(j,:,:));
                tmpdata=permute(tmpdata,[2,1]);
                [rho,~] = corr(behavior,tmpdata);
                T(j,:)=rho;
           end
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
        fsize=size(avgdata1,1);
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
            avgdata1_roi=avgdata1;
            avgdata1_imx=permute(avgdata1_roi, [3,1,2]);
            Results = ept_TFCE_CD(avgdata1_imx, behavior,[] , 'nPerm',n_perm,'type','c','chn',tfce_nbrs);
            T=Results.Obs(:,:);
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
            avgdata1_roi=avgdata1(:,cluster_id,:);
            avgdata1_imx=permute(avgdata1_roi, [3,1,2]);
            Results = ept_TFCE_CD(avgdata1_imx, behavior,[] , 'nPerm',n_perm,'type','c','chn',tfce_nbrs);
            %[T df] = ttest_cell({ avgdata1 avgdata2 });
            T=zeros(size(avgdata1,1),size(avgdata1,2));
            P=ones(size(avgdata1,1),size(avgdata1,2));
            for j=1:size(avgdata1,1)
                tmpdata=squeeze(avgdata1(j,:,:));
                tmpdata=permute(tmpdata,[2,1]);
                [rho,~] = corr(behavior,tmpdata);
                T(j,:)=rho;
            end
            test1=isnan(T);
            i_nan = find(test1);
            T(i_nan)=0;
            P(:,cluster_id)=Results.P_Values(:,:);
        end
    end
end

clim=[-0.6,0.6];

cfg.clim=clim;
cfg.xmin=xmin;
cfg.xmax=xmax;
cfg.ymin=ymin;
cfg.ymax=ymax;
cfg.fF=fF;
cfg.dyadic=dyadic;
% colorbar('FontSize',14,'YTick',-100:50:100,'YTickLabel',{'-100%','-50%','base','50%','100%'});