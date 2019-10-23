clear all 
%Sarah's code modified by DBH

disp(['Start: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

%% Import data

% No filter for N1; 16kHz data from Rambach
fn = 'TimeTagsOut-2016-05-30-12-59-20.dat';         
window = 5000; 
filtwin = 250;
winh = window/2;

% depth to perform 'nearest' search
dep = 2;

% bin size in ns
binsize = 0.0823045; 
% Window half in bins, used one to the right and one to the left for full window
W = floor(window/binsize/2); 
    
% Reading out channel and timetag information
[channel,~,tt] = ReadTimetags(fn);          
disp('Loading done')    

%%  Separating channels 

i = channel == 0; tt_1 = tt(i);         %% CH 0 timetags, trigger
j = channel == 1; tt_ch2 = tt(j);      %% CH 1 timetags
k = channel == 2; tt_ch3 = tt(k);      %% CH 2 timetags

% Smaller dataset
% dsize = 5000000;
tt_1 = tt_1(1:end);
tt_ch2 = tt_ch2(1:end);
tt_ch3 = tt_ch3(1:end);

% Convert from uint64 to double
tt_1i = double(tt_1);
tt_2i = double(tt_ch2);
tt_3i = double(tt_ch3);

% Singles counts
N1 = length(tt_1i);
N2 = length(tt_2i);
N3 = length(tt_3i);

clear tt_1 tt_ch2 tt_ch3 tt i j k channel

%% Raw coincidence analysis with 'Nearestpoint'

disp(['Starting coincidences: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

% t = indices of tt corresponding to coincidence event
% Diff = time difference of coincidence event (within W)

% N12
[N12, t12, Diff12] = getCoinc(tt_1i, tt_2i, W, dep);
disp(['N12 done: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

% N13
[N13, t13, Diff13] = getCoinc(tt_1i, tt_3i, W, dep);
disp(['N13 done: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

% N23
[N23, t23, Diff23] = getCoinc(tt_2i, tt_3i, W, dep);
disp(['N23 done: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

% N32
[N32, t32, Diff32] = getCoinc(tt_3i, tt_2i, W, dep);
disp(['N32 done: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

%% N1223 time differences between heralded ch2 and ch3 events
 
% comt = value of common element
% comt_i = index of common element 

[~, comt2i12, comt2i23] = triples(t12(:,2), t23(:,1));

[~, tripi1223b, Diff1223b] = getCoinc(tt_1i(t12(comt2i12,1)), tt_3i(t23(comt2i23,2)), W, dep);
% t1223b contains indicies of tt_1i(t12(comt2i12,1)), essentially comt2i12,
% that were in coincidence with a t3 (ie. not a simple index translation)
t1223b = [t12(comt2i12(tripi1223b(:,1)),1), t23(comt2i23(tripi1223b(:,2)),2)];
t1223a = [t12(comt2i12(tripi1223b(:,1)),1), t12(comt2i12(tripi1223b(:,1)),2)];
% so t1223 contains relevant indices of tt for triples
Diff1223a = Diff12(comt2i12(tripi1223b(:,1)));
N1223 = length(Diff1223a);
clear comt2i12 comt2i23
disp(['N1223 done: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

% N1332 time differences between heralded ch2 and ch3 events
[~, comt3i13, comt3i32] = triples(t13(:,2), t32(:,1));

[~, tripi1332b, Diff1332b] = getCoinc(tt_1i(t13(comt3i13,1)), tt_2i(t32(comt3i32,2)), W, dep);
t1332b = [t13(comt3i13(tripi1332b(:,1)),1), t32(comt3i32(tripi1332b(:,2)),2)];
t1332a = [t13(comt3i13(tripi1332b(:,1)),1), t13(comt3i13(tripi1332b(:,1)),2)];
Diff1332a = Diff13(comt3i13(tripi1332b(:,1)));
N1332 = length(Diff1332a);
clear comt3i13 comt3i32
disp(['N1332 done: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

clear tt_1i tt_2i tt_3i

% % histogram these time differences for approximate heralded g2
% figure(1); histogram(dif1213, 100)

% % unheralded source
% [N23, t23, dif23] = getCoinc(tt_2i, tt_3i, W, dep);
% figure(2); histogram(dif23, 100)

%% Histograms

% N12
[Hist12, offst12, dt12] = histoffset(Diff12, binsize);

% N13
[Hist13, offst13, dt13] = histoffset(Diff13, binsize);

% N23
[Hist23, offst23, dt23] = histoffset(Diff23, binsize);

% N32
[Hist32, offst32, dt32] = histoffset(Diff32, binsize);

% N1223
[Hist1223a, offst1223a, dt1223a] = histoffset(Diff1223a, binsize);
[Hist1223b, offst1223b, dt1223b] = histoffset(Diff1223b, binsize);

% N1332
[Hist1332a, offst1332a, dt1332a] = histoffset(Diff1332a, binsize);
[Hist1332b, offst1332b, dt1332b] = histoffset(Diff1332b, binsize);

% clearvars -except binsize filtwin dt12 dt13 dt23 dt32 dt1223 dt1332 ...
%     offst12 offst13 offst23 offst32 offst1223 offst1332 ...
%     Diff12 Diff13 Diff23 Diff32 Diff1223 Diff1332

%% Filtering step: only take coincidences within integer multiples of 1/FSR

FSR = 8.28/binsize;     % FSR in bins
ds = 12;                % Bin width over peaks (13 = 1.07ns)
fwin = floor(filtwin/(2*FSR*binsize));    % Coincidence window, halved for +/-, binsize to cancel FSR

% Filtered N12(0)
nmax12 = floor(dt12/2/FSR);

m = -nmax12;
peakbins12 = (round(m*FSR-ds/2)+offst12:round(m*FSR+ds/2)+offst12)';
% column vector of time bins around left-most peak m

for m = -nmax12+1:nmax12
    peakbins12 = [peakbins12 ; ...
        (round(m*FSR-ds/2)+offst12:round(m*FSR+ds/2)+offst12)']; 
                                  % vector with all 'good' bins
end

% N12f = sum(sum(peakbins12 == Diff12)); 
% sums the logical for whenever Diff12 falls in a peakbin
N12f=0;
    for J = 1:length(peakbins12)
        N12f = sum(peakbins12(J) == Diff12)+N12f;
    end

% clear peakbins12
disp(['N12f done: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

% Filtered N13(0)
nmax13 = floor(dt13/2/FSR);

m = -nmax13;
peakbins13 = (round(m*FSR-ds/2)+offst13:round(m*FSR+ds/2)+offst13)';

for m = -nmax13+1:nmax13
    peakbins13 = [peakbins13 ; ...
        (round(m*FSR-ds/2)+offst13:round(m*FSR+ds/2)+offst13)'];
end

% N13f = sum(sum(peakbins13 == Diff13)); 
N13f=0;
    for J = 1:length(peakbins13)
        N13f = sum(peakbins13(J) == Diff13)+N13f;
    end

% clear peakbins13
disp(['N13f done: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

% Filtered N1223(0)
parfor i = 1:length(Diff1223a)
    if (sum(Diff1223a(i) == peakbins12) > 0) && (sum(Diff1223b(i) == peakbins13) > 0)
        n1223f(i) = 1;
    else 
        n1223f(i) = 0;
    end
end
N1223f = sum(n1223f);
clear n1223f

disp(['N1223f done: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

% Filtered N1332(0)
parfor i = 1:length(Diff1332a)
    if (sum(Diff1332a(i) == peakbins13) > 0) && (sum(Diff1332b(i) == peakbins12) > 0)
        n1332f(i) = 1;
    else 
        n1332f(i) = 0;
    end
end
N1332f = sum(n1332f);
clear n1332f

disp(['N1332f done: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

%% Filtered N13(t)

% Zero = t in ns, zero point around which coincidences are taken, moving in FSRs  
Zero = (-181*8.28:8.28:1500);          
% Zero = (-120*8.28:8.28:1000);  % smaller range for less data
N13t = zeros(length(Zero),1);

for Zer = 1:length(Zero)    
    m = -(fwin);
    T3{Zer} = (round((m+Zero(Zer)-1)*FSR-ds/2)+offst13:round((m+Zero(Zer)-1)*FSR+ds/2)+offst13)';
           
    for m = -fwin+1:fwin
        T3{Zer} = [ T3{Zer} ; ...
            (round((m+Zero(Zer)-1)*FSR-ds/2)+offst13:round((m+Zero(Zer)-1)*FSR+ds/2)+offst13)' ];   % Vector filled with GOOD time bins around Zero(in units of FSR)              
    end
            
parfor i = 1:length(Diff13)
    if sum(Diff13(i) == T3{Zer}) > 0
        n13t(i) = 1;
    else 
        n13t(i) = 0;
    end
end

N13t(Zer) = sum(n13t);
clear n13t

%     if mod(Zer,10) == 0
%         Zer
%     end
            
end

disp(['N13t done: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

%% Filtered N1223(t)

for Zer = 1:length(Zero)    
    m = -(fwin);
    T2{Zer} = (round((m+Zero(Zer)-1)*FSR-ds/2)+offst12:round((m+Zero(Zer)-1)*FSR+ds/2)+offst12)';
           
    for m = -fwin+1:fwin
        T2{Zer} = [ T2{Zer} ; ...
            (round((m+Zero(Zer)-1)*FSR-ds/2)+offst12:round((m+Zero(Zer)-1)*FSR+ds/2)+offst12)' ];   % Vector filled with GOOD time bins around Zero(in units of FSR)              
    end
            
parfor i = 1:length(Diff1223a)
    if (sum(Diff1223a(i) == T2{Zer}) > 0) && (sum(Diff1223b(i) == T3{Zer}) > 0)
        n1223t(i) = 1;
    else 
        n1223t(i) = 0;
    end
end

N1223t(Zer) = sum(n1223t);
clear n1223t

end

clear T

disp(['N1223t done: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

%% Filtered N1332(t)

for Zer = 1:length(Zero)    
                
parfor i = 1:length(Diff1332a)
    if (sum(Diff1332a(i) == T3{Zer}) > 0) && (sum(Diff1332b(i) == T2{Zer}) > 0)
        n1332t(i) = 1;
    else 
        n1332t(i) = 0;
    end
end

N1332t(Zer) = sum(n1332t);
clear n1332t
    
end

disp(['N1332t done: ',datestr(now,'dd-mm-yyyy HH:MM:SS')])

%% g2(0)_a - N1223

g2_a = N1*N1223f/(N12f*N13f);
g2err_a = g2_a*((N1)^-1+(N1223f)^-1+(N12f)^-1+(N13f)^-1)^0.5;
g2raw_a = N1*N1223/(N12*N13);
g2rawe_a = g2raw_a*((N1)^-1+(N1223)^-1+(N12)^-1+(N13)^-1)^0.5;

%% g2(t)_a - N1223

g2t_a = N1*N1223t./(N12f.*N13t);
g2terr_a = g2t_a.*((N1).^-1+(N1223t).^-1+(N12f).^-1+(N13t).^-1).^0.5;

%% g2(0)_b - N1332

g2_b = N1*N1332f/(N12f*N13f);
g2err_b = g2_b*((N1)^-1+(N1332f)^-1+(N12f)^-1+(N13f)^-1)^0.5;
g2raw_b = N1*N1332/(N12*N13);
g2rawe_b = g2raw_b*((N1)^-1+(N1332)^-1+(N12)^-1+(N13)^-1)^0.5;

%% g2(t)_b - N1332

g2t_b = N1*N1332t./(N12f.*N13t);
g2terr_b = g2t_b.*((N1).^-1+(N1332t).^-1+(N12f).^-1+(N13t).^-1).^0.5;

%% Saving the workspace

save('g2v2_16kHz_5000_f250')

%%
p = profile('info')
save myprofiledata p

