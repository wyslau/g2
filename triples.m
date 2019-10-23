function [comvalue, iv1, iv2] = triples(v1, v2)

[comvalue, iv1, iv2] = intersect(v1, v2,'stable');
% remove first set of common t2 times from t23 and look again for
% reptitions - first iteration
v2a = v2;
v2a(iv2,:) = [];
[comvaluea, iv1a, iv2a] = intersect(v1, v2a,'stable');      
% combine with initial iteration
comvalue = [comvalue; comvaluea];
iv1 = [iv1; iv1a];
iv2 = [iv2; iv2a];

% second iteration
v2b = v2a;
v2b(iv2a,:) = [];
[comvalueb, iv1b, iv2b] = intersect(v1, v2b,'stable');      
% combine with initial iteration
comvalue = [comvalue; comvalueb];
iv1 = [iv1; iv1b];
iv2 = [iv2; iv2b];
clear comvaluea iv1a iv2a

% third iteration - should be more than enough (check for high singles)
v2c = v2b;
v2c(iv2b,:) = [];
[comvaluec, iv1c, iv2c] = intersect(v1, v2c,'stable');      
% combine with initial iteration
comvalue = [comvalue; comvaluec];
iv1 = [iv1; iv1c];
iv2 = [iv2; iv2c];
clear comvalueb iv1b iv2b

% fourth iteration - should be more than enough (check for high singles)
v2d = v2c;
v2d(iv2c,:) = [];
[comvalued, iv1d, iv2d] = intersect(v1, v2d,'stable');      
% combine with initial iteration
comvalue = [comvalue; comvalued];
iv1 = [iv1; iv1d];
iv2 = [iv2; iv2d];
clear comvaluec iv1c iv2c

% fifth iteration - should be more than enough (check for high singles)
v2e = v2d;
v2e(iv2d,:) = [];
[comvaluee, iv1e, iv2e] = intersect(v1, v2e,'stable');      
% combine with initial iteration
comvalue = [comvalue; comvaluee];
iv1 = [iv1; iv1e];
iv2 = [iv2; iv2e];
clear comvalued iv1d iv2d
clear comvaluee iv1e iv2e

% this doesn't work as it doesn't keep track of the triple events
% individually:
% t12com = ismember(t12(2,:),comt2);          % find all common t2 times in t12/t23 (ie. look for repetitions)
% t23com = ismember(t23(1,:),comt2);          % note: 'ismember' returns logic
% t12comi = find(t12com);                 % finds index of all non-zero elements of logic 
% t23comi = find(t23com);
% clear t12com t23com
% t1223a = [tt_1i(t12(1,t12comi)), tt_2i(t12(2,t12comi))]';
% Diff1223a = t1223a(1,:) - t1223a(2,:);

% this doesn't work: overcounts by a fair bit...:
% [n1223ai, t1223ai, Diff1223ai] = getCoinc(tt_1i(t12(1,t12comi)), tt_2i(t12(2,t12comi)), W, dep);
% [~, t1223b, Diff1223b] = getCoinc(tt_1i(t12(1,t12comi)), tt_3i(t23(2,t23comi)), W, dep);
% [n1223old, t1223old, Diff1223old] = getCoinc(tt_1i(t12(1,:)), tt_3i(t23(2,:)), W, dep);

% % original N1332 time differences between heralded ch2 and ch3 events
% [N1332, t1332, Diff1332] = getCoinc(tt_1i(t13(1,:)), tt_2i(t32(2,:)), W, dep);
% disp(['N1332 done: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])
