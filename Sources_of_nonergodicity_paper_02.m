% This code is shared as supporting material for the paper: 
% 
% "Sources of nonergodicity for teleconnections as cross-correlations" by
% Bodai, Lee, Aneesh
%
% specifically, for Figure 3 b. 

% Author: Tamas Bodai (bodai at pusan.ac.kr, bodait at yahoo.com)
% Date: 07.02.22. 


clear
close all

load ENSO_IM_teleco_11_aismr_nino3

% Sample values for tau
tau_vec = 3:2:31;
ntau = length(tau_vec);

doplot = 0;

wsm = 0; % which smoothing technique to apply; 1: movav alone, 0: movav & sgolay, or, sgolay meant with sg_ord > 1

% Realisations to keep 
r2k = [ 1:63]; % keep all

aismr = aismr(:,r2k);
nino3 = nino3(:,r2k);

sg_ord = 1;

nr = size(nino3,2); % number of realisations
ny = size(nino3,1); % * years
yr = 1880 + (1:ny); 

% Ensemble means
nino_ea = mean(nino3,2);
ismr_ea = mean(aismr,2);

% 20.07.21. The new loop
eps_n = zeros(1,ntau);
eps_a = zeros(1,ntau);
for i_tau = 1:ntau

ws = tau_vec(i_tau); 
b = (1/ws)*ones(1,ws);
ws2 = ws;
ws3 = ws2;

nino_sm = nino3*0;
ismr_sm = aismr*0;
if doplot
fh_ni = figure; hold on
fh_ai = figure; hold on
end
for i_r = 1:nr % when doplot=1, pick a realisation for looking at the result of smoothing
    i_r
% Copied from Reproduce_Kumar_1999_02.m 
nino = nino3(:,i_r)';
ismr = aismr(:,i_r)';

% Is there a need for moving average? Maybe not, but filtering twice gives
% a smoother signal.
nino_sg = sgolayfilt(nino, sg_ord, ws); 
ismr_sg = sgolayfilt(ismr, sg_ord, ws);
    % With sg_ord=1 e.g. nino_sg is the same as nino_ma, except that it
    % goes the full data span!!

nino_sg2 = sgolayfilt(nino_sg, sg_ord, ws2); % 22.01.21. ws -> ws2
ismr_sg2 = sgolayfilt(ismr_sg, sg_ord, ws2); % *

% Let's not use movav done by filter, only sgolay, so that we don't
% contract the data span
if     wsm 
    nino_sm(:,i_r) = nino_sg;
    ismr_sm(:,i_r) = ismr_sg;
else
    nino_sm(:,i_r) = nino_sg2;
    ismr_sm(:,i_r) = ismr_sg2;
end

end

% Low and high freq anomalies
nino_ano_lf = nino_sm - mean(nino_sm,2); 
nino_ano_hf = nino3   - nino_sm; 

ismr_ano_lf = ismr_sm - mean(ismr_sm,2); 
ismr_ano_hf = aismr   - ismr_sm;

% And the full anomalies
nino_ano_fu = nino3   - nino_ea; 
ismr_ano_fu = aismr   - ismr_ea;

eps_n_e = var(nino_ano_lf,[],2)./var(nino_ano_hf,[],2);
eps_a_e = var(ismr_ano_lf,[],2)./var(ismr_ano_hf,[],2);
eps_n(i_tau) = mean(eps_n_e(ws2:end-ws2));
eps_a(i_tau) = mean(eps_a_e(ws2:end-ws2));

if ws == 21
    figure; plot(yr,eps_n_e,yr,eps_a_e); 
    xlabel('year')
end
end

figure; hold on
%subplot(1,2,1); 
plot(tau_vec,eps_n,'o-')
%subplot(1,2,2); 
plot(tau_vec,eps_a,'o-')
plot([10 100],[0.1 0.01],'Color','k')
legend({'Nino3';'AISMR';'slope -1'})
xlabel('\tau')
ylabel('\epsilon^2')
set(findall(gcf,'-property','FontSize'),'FontSize',16)