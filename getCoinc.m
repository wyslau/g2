function [N, t, dif] = getCoinc(v1, v2, W, depth)
%% Given two vectors v1, v2 and window W, getCoinc returns the number of 
%coincidences within the given window. It returns the
%elements of v2 that are 'coincidences' as t (this is asymmetric) 

it=0;
t_ha = [];
t_hb = [];
t1_ha = [];
t1_hb = [];
d_ha = [];
d_hb = [];

while it <= depth

    % find indices (tind) of closest element of v2 for every v1, above and
    % below v1. Record the absolute difference (diff) each time.
    [tinda, diffa] = nearestpointDBH(v1,v2,'previous',it);
    [tindb, diffb] = nearestpointDBH(v1,v2,'next',it);

    % return only heralded events indices from ch2 within W with 
    % corresponding diff
    t_ha = [t_ha; tinda(diffa <= W)];
    t_hb = [t_hb; tindb(diffb <= W)];
        
    d_ha = [d_ha; diffa(diffa <= W)];
    d_hb = [d_hb; diffb(diffb <= W)];
    
    % find corresponding timestamps for heralds from ch1
    t1_ha = [t1_ha; find(diffa <= W)];    
    t1_hb = [t1_hb; find(diffb <= W)];    

    %iterate iterator
    it=it+1;

end

% % unique heralded ch2 events
% [t, IA, ~] = unique([t_ha ; t_hb]);
% %corresponding time difference from herald (multiply 'befores' by -1
% %because nearestpoint returns only absolute differences
% tst = [d_ha ; -1*d_hb];
% dif = tst(IA);
% 
% %sort output (this is important for further nearest finding)
% [t, IA] = sort(t);
% dif = dif(IA);
% 
% %two ways of counting coincidences in this channel
% N = length(t_ha) + length(t_hb);
% %N = length(t);

t(:,1) = [t1_ha; t1_hb];
t(:,2) = [t_ha; t_hb];
dif = [d_ha; -1*d_hb];

N = length(dif);
end