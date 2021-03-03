function [Omega_time_spec]=produce_Omega_time(All_nt)
% All_nt: fres x Rreq x samples x imfs
All_size=size(All_nt);
Omega_time_spec=zeros(All_size(2),All_size(3));
% for i_imf=1:All_size(4)
%     for i_fre=1:All_size(1)
tmpE1=squeeze(sum(All_nt,1)); % squeeze(tmpE1);
if numel(All_size)==4 && All_size(4)>1
   tmpE1=squeeze(sum(tmpE1,3));
end
%     end
% end
Omega_time_spec=tmpE1;
%log_Omega_time_spec=log10(Omega_time_spec);