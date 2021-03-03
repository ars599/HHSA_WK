function [T,P,cfg,df,avgdata1]=perform_onesample_tw_power2(D,cond_cell,channels,time_roi,idFreqs,S)
% sigbool must be a 3x1 boolean vector corresponding to the boundary of
% p value equal to [0.05, 0.01, 0.001]
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
structD=struct(D);
if strcmp(D.transformtype, 'TF') && isfield(D,'wtc')
    idFreqs=1:D.nfrequencies;
end


if length(cond_cell)>1
   %choice = questdlg('for one sample t-test, you can only test one condition. test first or concatenation','One sample t-test limit','first','concatenation','no test','concatenation');
   choice= 'concatenation';
   switch choice
        case 'first'
             %idx1=D.pickconditions(cond_cell{1});
             idx1=indtrial(D,cond_cell{1},'GOOD');
             tmpdata1=D(channels,idFreqs,slice_roi(1):slice_roi(2),idx1);
        case 'concatenation'
              %idx1=D.pickconditions(cond_cell{1});
            idx1=indtrial(D,cond_cell{1},'GOOD');  
            if weighted
                repl_1=getRepl(structD, idx1);
            else
                repl_1=ones(1,length(idx1));
            end
            repl_1_4d(1,1,1,:)=repl_1;
            tmpdata1=D(channels,idFreqs,slice_roi(1):slice_roi(2),idx1).*repmat(repl_1_4d,[length(channels),length(idFreqs),length(slice_roi(1):slice_roi(2)),1]);
            c1=repl_1;
            for h=2:size(cond_cell,1)*size(cond_cell,2)
                %idxh1=D.pickconditions(cond_cell{h});
                idxh1=indtrial(D,cond_cell{h},'GOOD');
%                 idx1=[idx1 idxh1];
                if weighted
                    repl_1=getRepl(structD, idxh1);
                else
                    repl_1=ones(1,length(idxh1));
                end
                repl_1_4d(1,1,1,:)=repl_1;
                tmpdata1=tmpdata1+D(channels,idFreqs,slice_roi(1):slice_roi(2),idxh1).*repmat(repl_1_4d,[length(channels),length(idFreqs),length(slice_roi(1):slice_roi(2)),1]);
                c1=c1+repl_1;
            end
            c1_4d(1,1,1,:)=c1;
            tmpdata1=tmpdata1./repmat(c1_4d,[length(channels),length(idFreqs),length(slice_roi(1):slice_roi(2)),1]);
           
        case 'no test'
            return
   end
elseif length(cond_cell)==1
    %idx1=D.pickconditions(cond_cell{1});
    idx1=indtrial(D,cond_cell{1},'GOOD');
    tmpdata1=D(channels,idFreqs,slice_roi(1):slice_roi(2),idx1);
else
    return
end



if size(tmpdata1,1)==1
    avgdata1=tmpdata1;
else
    avgdata1=mean(tmpdata1);
end
avgdata1=shiftdim(avgdata1);
if ph_time
   avgdata1=repmat(avgdata1(:,2:(end-1),:),[1 2 1]); 
end
avgdata2=zeros(size(avgdata1));
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
if fF && S.stats
    yy=repmat(Freqs(idFreqs),1,size(avgdata1,2));
    xx=repmat(D.time(slice_roi(1):slice_roi(2)),size(avgdata1,1),1); 
    delbool=(yy>(xx-0.25));
    delbool_3d=repmat(delbool,1,1,size(avgdata1,3));
    avgdata1(delbool_3d)=0;
end
testdata=avgdata1-avgdata2;
df=size(testdata,ndims(testdata))-1;
if S.stats
    if ~ nonparam
            [T df] = ttest_cell({ avgdata1 avgdata2 });
            %P_masked=zeros(size(T));
            P=1-tcdf(abs(T),df);
    else
        if ~tfce
            fsize=size(testdata,1);
            chan_hood=zeros(fsize,fsize);
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
                    [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm1x(testdata_roi,chan_hood,n_perm,0.05,nonparam_tail,alpha);
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
            if isempty(cluster_roi)
    %             avgdata1_imx=permute(avgdata1, [3,2,1]);
    %             avgdata2_imx=permute(avgdata2, [3,2,1]);
                avgdata1_roi=shiftdim(avgdata1,-1);
                avgdata2_roi=shiftdim(avgdata2,-1);
                avgdata1_imx=permute(avgdata1_roi, [4,1,2,3]);
                avgdata2_imx=permute(avgdata2_roi, [4,1,2,3]);
                avgdata1_imx=repmat(avgdata1_imx,[1,2,1,1]);
                avgdata2_imx=repmat(avgdata2_imx,[1,2,1,1]);
                tfce_nbrs=[1,2;1,2]; % two psudo channels
                Results = ept_TFCE(avgdata1_imx, avgdata2_imx,[] , 'nPerm',n_perm,'type','d','chn',tfce_nbrs);
                [~, df] = ttest_cell({ avgdata1 avgdata2 });
                T=shiftdim(Results.Obs(1,:,:),1);
                test1=isnan(T);
                i_nan = find(test1);
                T(i_nan)=0;
                P=shiftdim(Results.P_Values(1,:,:),1);
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
    %             avgdata1_roi=avgdata1(:,cluster_id,:);
    %             avgdata2_roi=avgdata2(:,cluster_id,:);
    %             avgdata1_imx=permute(avgdata1_roi, [3,2,1]);
    %             avgdata2_imx=permute(avgdata2_roi, [3,2,1]);
                avgdata1_roi=shiftdim(avgdata1(:,cluster_id,:),-1);
                avgdata2_roi=shiftdim(avgdata2(:,cluster_id,:),-1);
                avgdata1_imx=permute(avgdata1_roi, [4,1,2,3]);
                avgdata2_imx=permute(avgdata2_roi, [4,1,2,3]);
                avgdata1_imx=repmat(avgdata1_imx,[1,2,1,1]);
                avgdata2_imx=repmat(avgdata2_imx,[1,2,1,1]);
                tfce_nbrs=[1,2;1,2]; % two psudo channels
                Results = ept_TFCE(avgdata1_imx, avgdata2_imx,[] , 'nPerm',n_perm,'type','d','chn',tfce_nbrs);
                [T df] = ttest_cell({ avgdata1 avgdata2 });
                P=ones(size(T));
                %T=t_orig;
                %T=Results.Obs';
                test1=isnan(T);
                i_nan = find(test1);
                T(i_nan)=0;
                P(:,cluster_id)=shiftdim(Results.P_Values(1,:,:),1);
            end
        end
    end
else
    T=[];
    P=[];
end
t_max = tinv(0.999,df);
t_min=-1*t_max;
%clim=[t_min,t_max];
if S.stats
    if isfield(D,'rescale_type')
        switch D.rescale_type
            case 'rel'
                clim=[-30,30];
            case 'logr'
                clim=[-2,2];
            otherwise
                clim=[-3,3];
        end
    else
        clim=[-30,30]; 
    end
else
    clim=0.5*[(-1)*max(avgdata1(:)), max(avgdata1(:))];
end
cfg.clim=clim;
cfg.xmin=xmin;
cfg.xmax=xmax;
cfg.ymin=ymin;
cfg.ymax=ymax;
cfg.fF=fF;
cfg.dyadic=dyadic;




function repl=getRepl(structD, idx)
    repl=zeros(1,length(idx));
    for k=1:length(repl)
        repl(k)=structD.trials(1,idx(k)).repl;
    end

% colorbar('FontSize',14,'YTick',-100:50:100,'YTickLabel',{'-100%','-50%','base','50%','100%'});