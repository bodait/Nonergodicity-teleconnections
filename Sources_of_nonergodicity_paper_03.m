% This code is shared as supporting material for the paper: 
% 
% "Sources of nonergodicity for teleconnections as cross-correlations" by
% Bodai, Lee, Aneesh
%
% specifically, for Figures S1 and other results reported in Part 1 of the
% Supplementary Information inline.

% Author: Tamas Bodai (bodai at pusan.ac.kr, bodait at yahoo.com)
% Date: 07.02.22. 


clear 
close all

% Number of Monte Carlo experiments for establishing H0. Go for 1e4 on a
% big computer, 1e3 on a laptop.
MC = 1e3; 

AorC = 1; % AISMR (1) or CRU (0) 
yf_aismr = 2014; % final year for AISMR data

if AorC
 load Reproduce_Kumar_1999_02_AISMR_ismr_nino
 y0 = 1880-1;
else
 load Reproduce_Kumar_1999_02_CRU_ismr_nino
 y0 = 1901-1;
end

% Switch for restricting the study period to match (Krishna-Kumar et al.
% 1999) 
KK = 0;

nino = detrend(nino',2); 

% Window size etc. used just for smoothing/anomalising 
ws = 21; 
b = (1/ws)*ones(1,ws);
sg_ord = 1;
% Window for a short period where r_tau is suspected to be significantly
% different. 
ws2_vec = 3:2:37;%25; % 28.07.21. These need to be odd now because of sgloay in 'anomalise'
nws2 = length(ws2_vec);

% The null-hipothesis H0 is a stationary r derived from obs'... but that r
% is uncertain!!                                                           
% For now, just get this with a largish ws
% nino_ano = anomalise(nino,b,sg_ord,ws);
% ismr_ano = anomalise(ismr,b,sg_ord,ws);
% Involving the ends too will yield a bit larger p-values
  nino_ano = anomalise_2(nino,sg_ord,ws); 
  ismr_ano = anomalise_2(ismr,sg_ord,ws);
[r_na,P,RL,RU] = corrcoef(nino_ano,ismr_ano); r_na = r_na(1,2);
temp = atanh(RU) - atanh(RL);
% Standard deviation of the Fisher transform
dz = temp(1,2)/4; % can just assign 0 here to see what's the difference

% Data period for the raw data
% AISMR: 1880-2014 % Note: data available untill 2019 now
% CRU:   1901-2010

% KK99 worked with data till 1997; let's remove data after that to see if
% KK could have obtained significance. Do only with the AISMR data.
if KK && AorC
ismr_ano = ismr_ano(1:end-(yf_aismr-1997));
nino_ano = nino_ano(1:end-(yf_aismr-1997));
nino     = nino    (1:end-(yf_aismr-1997)); 
end

tsl0 = length(nino); % 'tsl' for 'time series length'
tsl  = length(nino_ano);

yr = (1:tsl)+y0;

% First, consider the "one parameter test", T0, the parameter being the
% limiting year for the unequal partition

p_c20thc_K_vec = zeros(1,tsl)*NaN;
c20thc_K_1     = zeros(1,tsl)*NaN;
c20thc_K_2     = zeros(1,tsl)*NaN;
std_sum        = zeros(1,tsl)*NaN;
% The Fisher law is not accurate for too small tau, so, apply a minimum
tau_min = 10;
for i1 = tau_min:tsl-tau_min
    temp = corrcoef(ismr_ano(   1:i1 ),nino_ano(   1:i1 ));
    c20thc_K_1(i1) = temp(1,2);
    temp = corrcoef(ismr_ano(i1+1:end),nino_ano(i1+1:end));
    c20thc_K_2(i1) = temp(1,2);
    c20thc_K_h1_h2_diff = atanh(c20thc_K_1(i1)) - atanh(c20thc_K_2(i1));
    std_sum(i1) = sqrt(1/(i1-3) + 1/(tsl - i1 -3));
    p_c20thc_K_vec(i1) = 2*(1 - normcdf(abs(c20thc_K_h1_h2_diff),0,std_sum(i1))); 
end

% Figure S1 d
figure; plot(yr,p_c20thc_K_vec,yr,c20thc_K_1,yr,c20thc_K_2,yr,std_sum)

p_min_3_obs = nanmin(p_c20thc_K_vec);


% Second, check the detection rate under NOT H0, but a hockey-stick-like
% forced change, with two FIXED partitions, the case of equal partitions,
% and another case when a partition is restricted to the head of the hockey
% stick. This is to show that the second test T2 is more powerful than the
% first T1 under the said conditions... In the same time i'm showing also
% what my intention was with detectability_by_YT.m: even if there is a
% considerable forced change in the end of the 20th c., the detectability
% measured by the detection rate is meager.

% (Actually, i do it after the Third point below, which is a similar piece
% of code.)

x = nino/std(nino);
x_ano = nino_ano;
tsl = length(x_ano);

% Hockey stick corr coeff time series. Like in detectability_by_YT.m
sotw = 30; % start of the teleconnection weakening: years before the end
dr = 0.2; % delta r in a period of [sotw] years
r_hs = -0.5*ones(1,tsl);
r_hs(end-sotw+1:end) = r_hs(end-sotw+1:end) + (1:sotw)/sotw*dr;

dtsc = zeros(MC,2); % detection of significant change
for i2 =  1:MC
    if mod(i2,100) == 0
        i2
    end
y = x .* sqrt(r_hs.^2./(1-r_hs.^2)) .* sign(r_hs) + randn(1,tsl0);

y_ano = anomalise_2(y,sg_ord,ws);

pdl = [floor(tsl0/2) tsl0-sotw]; % partition delimiting year
p_c20thc_K_vec = zeros(1,2)*NaN;
for i1 = 1:2
    temp = corrcoef(y_ano(   1:pdl(i1) ),x_ano(   1:pdl(i1) ));
    c20thc_K_1 = temp(1,2);
    temp = corrcoef(y_ano(pdl(i1)+1:end),x_ano(pdl(i1)+1:end));
    c20thc_K_2 = temp(1,2);
    c20thc_K_h1_h2_diff = atanh(c20thc_K_1) - atanh(c20thc_K_2);
    std_sum = sqrt(1/(pdl(i1)-3) + 1/(tsl - pdl(i1) -3));
    p_c20thc_K_vec(i1) = 2*(1 - normcdf(abs(c20thc_K_h1_h2_diff),0,std_sum)); 
end
dtsc(i2,:) = p_c20thc_K_vec<0.05; 
end

% Detection rates
dtr_2 = sum(dtsc)/MC; % >5% indeed 


% Third, check the detection rate under H0 when allowing any partition,
% and plot also the distribution of p-values

dtsc = zeros(MC,1); % detection of significant change
p_min_3 = dtsc; 
p_3 = zeros(MC,tsl)*NaN; % Check again detection rate with fixed sample (test T0)
fh1 = figure;
fh2 = figure;
for i2 =  1:MC
    if mod(i2,100) == 0
        i2
    end
%x = randn(1,tsl0);
y = x * sqrt(r_na^2/(1-r_na^2)) * sign(r_na) + randn(1,tsl0);

%x_ano = anomalise_2(x,sg_ord,ws); 
y_ano = anomalise_2(y,sg_ord,ws);

%r_tau = movmean(x_ano.*y_ano,ws)./movstd(x_ano,ws)./movstd(y_ano,ws);
%figure(fh1); plot(r_tau)

tsl = length(x_ano);

p_c20thc_K_vec = zeros(1,tsl)*NaN;
for i1 = tau_min:tsl-tau_min
    temp = corrcoef(y_ano(   1:i1 ),x_ano(   1:i1 ));
    c20thc_K_1 = temp(1,2);
    temp = corrcoef(y_ano(i1+1:end),x_ano(i1+1:end));
    c20thc_K_2 = temp(1,2);
    c20thc_K_h1_h2_diff = atanh(c20thc_K_1) - atanh(c20thc_K_2);
    std_sum = sqrt(1/(i1-3) + 1/(tsl - i1 -3));
    p_c20thc_K_vec(i1) = 2*(1 - normcdf(abs(c20thc_K_h1_h2_diff),0,std_sum)); 
end
p_min_3(i2) = nanmin(p_c20thc_K_vec);
%figure(fh2); plot(p_c20thc_K_vec)
dtsc(i2) = sum(p_c20thc_K_vec<0.05)>0; % problem with multiple testing
    % This could be also: dtsc(i2) = p_min_3(i2)<0.05;
p_3(i2,:) = p_c20thc_K_vec; 
end

% Detection rate when choosing the test parameter freely, i.e., cheekily 
dtr_3 = sum(dtsc)/MC; % >5% indeed 

% Distribution of max'l p-values... Skewed towards small values indeed
% % Figure S1 e
figure; histogram(p_min_3,50)
title(['detection rate: ' num2str(dtr_3)])

% Then the revised p-value is:
pv_3 = 2*min(nansum(p_min_3 > p_min_3_obs),nansum(p_min_3 < p_min_3_obs))/MC;

% How does the detection rate for T0 depend on sample size; does it get
% more liberal as it decreases? ... Doesn't really look like that.
figure; plot(sum(p_3<0.05,1)/MC)


% Fourth, allow for a partition when one of the intervals is not compact
% but the other one is an interlude: test T2. This way we have TWO
% parameters of the test. What is not the adjusted p-value?! I expect it to
% be smaller because it is a more powerful test in the supposed situation
% that the forced change is nonmonotonic.

% The code is a modified version of the third point above

p_min_4 = zeros(MC,1); 
tau_opt_4 = zeros(MC,1);
tic
for i2 =  1:MC
%parfor i2 =  1:MC % no parfor here if we use parfor in min_p
    % IF we use parfor in min_p, then we can follow progress here properly,
    % being also make prediction when the code would finish, desirable for
    % large MC like 1e4.
    if mod(i2,100) == 0
        i2
        toc
    end
    y = x * sqrt(r_na^2/(1-r_na^2)) * sign(r_na) + randn(1,tsl0);

    y_ano = anomalise_2(y,sg_ord,ws);

    temp = min_p(x_ano,y_ano,tau_min);
    p_min_4(i2) = temp(1); 
    tau_opt_4(i2) = temp(2);
end

% Detection rate when choosing the test parameter freely, i.e., cheekily 
dtr_4 = sum(p_min_4<0.05)/MC; % >5% indeed 

% Distribution of max'l p-values... Skewed towards small values indeed
% Figure S1 e
figure; histogram(p_min_4,50)
title(['detection rate: ' num2str(dtr_4)])

% Them the revised p-value is:
temp = min_p(nino_ano,ismr_ano,tau_min);
p_min_4_obs = temp(1); 
tau_opt = temp(2); 
mwy_opt = temp(3)+y0; % mid window year
pv_4 = 2*min(nansum(p_min_4 > p_min_4_obs),nansum(p_min_4 < p_min_4_obs))/MC;


% Figure S1 b. As an illustration for the paper, plot the approx distr of
% r_tau for the tau's (ws2) used manually in min_p to plot r_tau time
% series.
tau = 20;
sig_F = 1/sqrt(tau-3); % "Fisher sigma"
z_na = atanh(r_na);
z_samp = linspace(z_na-5*sig_F,z_na+5*sig_F,1e3);
r_samp = tanh(z_samp);

p_of_z_tau = normpdf(z_samp,z_na,sig_F);
p_of_r_tau = p_of_z_tau ./ (1-r_samp.^2);

% Note for the transformation of rv
% syms r z
% z=atanh(r);
% diff(z)
%  
% ans =
%  
% -1/(r^2 - 1)

figure; semilogx(p_min_4,tau_opt_4,'.')
ylabel('\tau_{opt}')
xlabel('p_{T0}')
set(findall(gcf,'-property','FontSize'),'FontSize',16)

% Figure S1 f. An additional figure for the paper: shows that...
figure; plot(p_of_r_tau,r_samp)

save Sources_of_nonergodicity_paper_03_x

function out = min_p(x_ano,y_ano,ws2_min)

tsl = length(x_ano);
t_vec = 1:tsl;

% There is a max of the window length so that the length of t_vec_rst would
% also not be smaller than ws2_min
ws2_vec = ws2_min:tsl-ws2_min;

p_c20thc_K_vec = zeros(tsl,tsl)*NaN; % NaN is important because of nanmax below
% 04.08.21. For checking and figure making purposes
r_tau_mid_vec = p_c20thc_K_vec;
parfor i2 = ws2_vec
    ws2 = i2;
p_c20thc_K_i1  = zeros(1  ,tsl)*NaN; % auxiliary var' for parfor
r_tau_mid_i1 = p_c20thc_K_i1;
for i1 = 1:tsl
    try
    % Try out all periods of length ws2 standing apart from the rest.
    t_vec_mid = i1:i1+ws2-1;
    t_vec_rst = setdiff(t_vec,t_vec_mid);
    r_tau_mid = corrcoef(x_ano(t_vec_mid),y_ano(t_vec_mid)); 
    r_tau_mid = r_tau_mid(1,2);
    r_tau_rst = corrcoef(x_ano(t_vec_rst),y_ano(t_vec_rst)); 
    r_tau_rst = r_tau_rst(1,2);
    
    dr_tau = atanh(r_tau_mid) - atanh(r_tau_rst);
    
    std_sum = sqrt(1/(ws2-3) + 1/(tsl - ws2 -3));
    
    %p_c20thc_K_vec(i2,i1) = 2*(1 - normcdf(abs(dr_tau),0,std_sum)); 
    % 04.08.21. For parfor
     p_c20thc_K_i1(  :,i1) = 2*(1 - normcdf(abs(dr_tau),0,std_sum)); 
     
      r_tau_mid_i1(  :,i1) = r_tau_mid;
    catch
    end
end
p_c20thc_K_vec(i2,:) = p_c20thc_K_i1; % 04.08.21. For parfor
r_tau_mid_vec(i2,:) = r_tau_mid_i1;
end
%p_min = nanmin(p_c20thc_K_vec(:));
% 04.08.21. Determine also the optimal values of the two test parameters
[p_min,i2_min] = nanmin(p_c20thc_K_vec); % min in COLUMNS are found 
[p_min,i1_min] = nanmin(p_min);
i2_min = i2_min(i1_min);
out = [p_min, i2_min, ...
         ceil(i2_min/2) + t_vec(i1_min)]; % midwindow

% % Interesting to look at this for the obs realisation 
% %i2 = i2_min; % obviously
% i2 = 20; % and another interesting choice
% y0 = 1880;
% figure; plot(t_vec+y0-1+ceil(i2/2), r_tau_mid_vec(i2,:)) % Figure S1 a
% figure; plot(t_vec+y0-1+ceil(i2/2),p_c20thc_K_vec(i2,:)) % Figure S1 c
end

% Make this bit of code also a function so that the code is more efficient
% and easier to debug
function x_ano = anomalise(x,b,sg_ord,ws)
    x_ma = filter(b,1,x);
    x_cubicMA = sgolayfilt(x_ma(ws:end), sg_ord, ws);
    x_ano = x(ceil(ws/2):end-floor(ws/2)) - x_cubicMA;
end

% Alternatively, DO NOT DISCARD the end bits of the data of half the window
% length ws. But bear in mind that these are not as good quality
% anomalies!!
function x_ano = anomalise_2(x,sg_ord,ws)
    x_SG = sgolayfilt(x, sg_ord, ws);
    x_ano = x - x_SG;
end