% This code is shared as supporting material for the paper: 
% 
% "Sources of nonergodicity for teleconnections as cross-correlations" by
% Bodai, Lee, Aneesh
%
% specifically, for Figures 1, 2, 3a, 4, S2, Table 1. 

% Author: Tamas Bodai (bodai at pusan.ac.kr, bodait at yahoo.com)
% Date: 07.02.22. 

% Obtain function 'movingslope' from "Mathworks File Exchange":
% https://ch.mathworks.com/matlabcentral/fileexchange/16997-movingslope

clear
close all

load ENSO_IM_teleco_11_aismr_nino3

% Set doplot = 1 to look at the result of decomposition for a single
% realisation. In this case, use a break point after the plotting command.
doplot = 0;
% Which smoothing technique to apply; 1: movav alone, 0: movav & sgolay,
% or, sgolay meant with sg_ord > 1 
wsm = 0; 

% Realisations to keep - to check if r_lf is completely different
% for different "ensemble realisations"
r2k = [ 1:63]; % keep all
%r2k = [ 1:31];
%r2k = [32:62];

aismr = aismr(:,r2k);
nino3 = nino3(:,r2k);

% Run a moving average "filter" to identify a low-frequency variability
% which will serve as a basis for identifying anomalies 
% Parameters for "filter" to facilitate moving averaging
ws = 21; 
b = (1/ws)*ones(1,ws);
ws2 = ws; 
ws3 = ws2; % tau for running window temporal cc r_YT and r_KK

sg_ord = 1;

nr = size(nino3,2); % number of realisations
ny = size(nino3,1); % * years
yr = 1880 + (1:ny); % time vector

% Ensemble means
nino_ea = mean(nino3,2);
ismr_ea = mean(aismr,2);

% Smooth single realisations in order to operatively define LF variability
nino_sm = nino3*0;
ismr_sm = aismr*0;
if doplot
fh_ni = figure; hold on
fh_ai = figure; hold on
end
for i_r = 1:nr % when doplot=1, pick a realisation for looking at the result of smoothing
    i_r
nino = nino3(:,i_r)';
ismr = aismr(:,i_r)';

nino_ma = filter(b,1,nino);
nino_cubicMA = sgolayfilt(nino_ma(ws:end), sg_ord, ws2); 
nino_ano = nino(ceil(ws/2):end-floor(ws/2)) - nino_cubicMA;
    % Actually, i never used this var.

ismr_ma = filter(b,1,ismr);
ismr_cubicMA = sgolayfilt(ismr_ma(ws:end), sg_ord, ws2); 
ismr_ano = ismr(ceil(ws/2):end-floor(ws/2)) - ismr_cubicMA;

% Is there a need for moving average? Maybe not, but filtering twice gives
% a smoother signal.
nino_sg = sgolayfilt(nino, sg_ord, ws); 
ismr_sg = sgolayfilt(ismr, sg_ord, ws);
    % With sg_ord=1 e.g. nino_sg is the same as nino_ma, except that it
    % goes the full data span!!

nino_sg2 = sgolayfilt(nino_sg, sg_ord, ws2); 
ismr_sg2 = sgolayfilt(ismr_sg, sg_ord, ws2); 

% Let's not use movav done by filter, only sgolay, so that we don't
% contract the data span
if     wsm 
    nino_sm(:,i_r) = nino_sg;
    ismr_sm(:,i_r) = ismr_sg;
else
    nino_sm(:,i_r) = nino_sg2;
    ismr_sm(:,i_r) = ismr_sg2;
end

% Have a look at the currently treated realisation
if doplot
figure(fh_ni); plot(yr(ws:end)-floor(ws/2),nino_ma(ws:end)) 
figure(fh_ai); plot(yr(ws:end)-floor(ws/2),ismr_ma(ws:end))

figure(fh_ni); plot(yr(ws:end)-floor(ws/2),nino_cubicMA) 
figure(fh_ai); plot(yr(ws:end)-floor(ws/2),ismr_cubicMA)

figure(fh_ni); plot(yr,nino_sg)
figure(fh_ai); plot(yr,ismr_sg)

figure(fh_ni); plot(yr,nino_sg2); hold off
figure(fh_ai); plot(yr,ismr_sg2); hold off
end

end

% Low and high freq anomalies
nino_ano_lf = nino_sm - mean(nino_sm,2); 
% This is what we should insist on, i.e., that the HF anomaly is
% what people do and can only do in practice, not relying on an ensemble,
% just a single realisation:
nino_ano_hf = nino3   - nino_sm; 

ismr_ano_lf = ismr_sm - mean(ismr_sm,2); 
ismr_ano_hf = aismr   - ismr_sm;

% And the true full anomalies by definition
nino_ano_fu = nino3   - nino_ea; 
ismr_ano_fu = aismr   - ismr_ea;

% Check if the anomalies have zero mean
figure; plot(yr,mean(nino_ano_fu,2),yr,mean(nino_ano_lf,2),yr,-mean(nino_ano_hf,2))
    % Only the HF anomalies have no zero mean indeed, as only that is based
    % on an operational definition which is not actually a natural
    % definition 

    
% Ensemble-wise Pearson corr coeff's
% Vector of all the corr coeff's; apply the def of Pearson cc
r_fu = zeros(1,ny);
r_lf = zeros(1,ny);
r_hf = zeros(1,ny);
r_hl = zeros(1,ny); % of the same quantity: nino
r_lh = zeros(1,ny); % low freq nino to hf rain
for i1 = 1:ny
    temp = corrcoef(nino_ano_fu(i1,:),ismr_ano_fu(i1,:));
    r_fu(i1) = temp(1,2); 
    temp = corrcoef(nino_ano_lf(i1,:),ismr_ano_lf(i1,:));
    r_lf(i1) = temp(1,2); 
    temp = corrcoef(nino_ano_hf(i1,:),ismr_ano_hf(i1,:));
    r_hf(i1) = temp(1,2);
    temp = corrcoef(nino_ano_hf(i1,:),nino_ano_lf(i1,:));
    r_hl(i1) = temp(1,2);
    % 22.01.21.
    temp = corrcoef(nino_ano_lf(i1,:),ismr_ano_hf(i1,:));
    r_lh(i1) = temp(1,2);
end

% r_fu and r_hf looks very similar
temp = corrcoef(r_fu,r_hf);
r_fuhf = temp(1,2);


eps_n_e = var(nino_ano_lf,[],2)./var(nino_ano_hf,[],2);
eps_a_e = var(ismr_ano_lf,[],2)./var(ismr_ano_hf,[],2);
eps_n = mean(eps_n_e(ws2:end-ws2));
eps_a = mean(eps_a_e(ws2:end-ws2));

r_hf_sm = sgolayfilt(r_hf, sg_ord, ws2);    
r_lf_re = (r_fu'.*sqrt((1+eps_n_e).*(1+eps_a_e)) - r_hf')./sqrt(eps_n_e.*eps_a_e);
r_lf_re_sm = sgolayfilt(r_lf_re, sg_ord, ws2); 

% For Table 1 in paper
std(r_hf_sm)
% Take the term that appears added to r_hf_sm, and smooth it the very same
% way as r_hf_sm for comparability. Only that ...
r_lf_re2_sm = sgolayfilt(r_fu'.*sqrt((1+eps_n_e).*(1+eps_a_e)) - r_hf', sg_ord, ws2);
% ...this inherits the massive "end outliers" of eps_n_e, eps_a_e...
% figure; plot(yr,r_lf_re2_sm,yr,r_lf.*sqrt(eps_n_e.*eps_a_e)')
% ...so, remove before calculating std
std(r_lf_re2_sm(ceil(ws2/2):end-ceil(ws2/2))) 
std(r_lf)


% Figure 2 or Figure S2
fh = figure; 
subplot(2,2,2); plot(yr,eps_n_e,yr,eps_a_e) 
legend({'$\epsilon_x^2$'; '$\epsilon_y^2$'},'Interpreter','latex','Orientation','horizontal')
xlim([1880,2100])
subplot(2,2,1); plot(yr,r_fu,yr,r_hf) 
xlim([1880,2100])
legend({'$r$'; '$r_H$'},'Interpreter','latex','Orientation','horizontal')
subplot(2,2,3); plot(yr,r_lf_re,yr,r_lf_re_sm,yr,r_lf)
xlim([1880,2100])
xlabel('year')
legend({'$r_L^{\bullet}$'; ...
    '$r_L^{\bullet*}$'; '$r_L$'},'Interpreter','latex','Orientation','horizontal')
subplot(2,2,4); plot(yr,r_hf_sm-mean(r_hf_sm),yr,r_lf_re2_sm,...
    yr,r_lf'.*sqrt(eps_n_e.*eps_a_e))
xlim([1880,2100])
legend({'$r_H^*-\overline{r_H^*}$'; '$\rho_L^{\bullet*}$'; ...
    '$\rho_L$'},'Interpreter','latex','Orientation','horizontal')
xlabel('year')
set(findall(gcf,'-property','FontSize'),'FontSize',16)


% Moving correlation by direct calculation
r_YT = zeros(length(ismr)-ws3+1,nr);
% Krishna-Kumar's *, again, as far as we can guess 
r_KK = zeros(length(ismr)-ws3+1,nr);
for i2 = 1:nr
for i1 = 1:size(r_YT,1)
    temp = corrcoef(nino3      (i1:i1+ws3-1,i2),aismr      (i1:i1+ws3-1,i2)); 
    r_YT(i1,i2) = temp(1,2);
    temp = corrcoef(nino_ano_hf(i1:i1+ws3-1,i2),ismr_ano_hf(i1:i1+ws3-1,i2));
    r_KK(i1,i2) = temp(1,2);
end
end
  

% Figure 1 and Figure 4
fh = figure; subplot(2,1,1)
plot(yr,r_fu,'Color',[1 1 1]*0.7)
hold on
r_fu_sm = sgolayfilt(r_fu, sg_ord, ws2);
plot(yr,r_fu_sm) 
plot(yr(ceil(ws3/2):ceil(ws3/2)+size(r_YT,1)-1),mean(r_YT,2),...
     yr(ceil(ws3/2):ceil(ws3/2)+size(r_YT,1)-1),mean(r_KK,2),yr,r_hf_sm) 
xlim([1880,2100])
xlabel('year')
legend({'$r$'; '$r^*$'; '$\langle r_{\tau,d} \rangle$'; ...
    '$\langle r_{\tau,p} \rangle$';...
    '$r_H^*$'},'Interpreter','latex')
   

% Moving window slopes as "regression" coefficients
alph = movingslope(nino_ea,ws3);
beta = movingslope(ismr_ea,ws3);

biasfact = alph.*beta./sqrt(var(nino_ano_hf,[],2).*var(ismr_ano_hf,[],2))*ws3^2/12;

figure(fh)
subplot(2,1,2)
plot(yr(ceil(ws3/2):ceil(ws3/2)+size(r_YT,1)-1),mean(r_YT,2)-mean(r_KK,2),...
    yr,biasfact)

  
epsx2 = alph.^2*ws3^2/12./var(nino_ano_hf,[],2);
epsy2 = beta.^2*ws3^2/12./var(ismr_ano_hf,[],2);

r_YT_pad = zeros(size(epsx2))*NaN;
r_YT_pad(ceil(ws3/2):ceil(ws3/2)+size(r_YT,1)-1) = mean(r_YT,2);
r_KK_pad = zeros(size(epsx2))*NaN;
r_KK_pad(ceil(ws3/2):ceil(ws3/2)+size(r_YT,1)-1) = mean(r_KK,2);

r_YT_re = (r_KK_pad + sqrt(epsx2.*epsy2))./sqrt((1+epsx2).*(1+epsy2));

%figure; plot(yr,r_YT_pad,yr,r_YT_re)

figure(fh)
subplot(2,1,2); hold on; plot(yr,r_YT_re-r_KK_pad)
xlim([1880,2100])
xlabel('year')
legend({'$\langle r_{\tau,d} \rangle - \langle r_{\tau,p} \rangle$'; '$\epsilon_x\epsilon_y$';...
    '$\langle r_{\tau,d} \rangle^{\bullet} - \langle r_{\tau,p} \rangle$'},'Interpreter','latex')
set(findall(gcf,'-property','FontSize'),'FontSize',16)


% Calculate power spectra to see if Nino3 and AISMR are really
% more like white noise, such that eps~1/tau. Consider only
% 20th c.
% Figure 3a
pxx_nino = zeros(nr,1e3);
pxx_ismr = zeros(nr,1e3);
for i1 = 1:nr
    [pxx,f] = periodogram(nino_ano_fu(1:120,i1),[],[],1);
    if i1 == 1
        n_f = length(f);
        pxx_nino = pxx_nino(:,1:n_f);
        pxx_ismr = pxx_ismr(:,1:n_f);
    end
    pxx_nino(i1,:) = pxx;
    [pxx,f] = periodogram(ismr_ano_fu(1:120,i1),[],[],1);
    pxx_ismr(i1,:) = pxx;
end
figure; 
%subplot(1,2,1); plot(f,mean(pxx_nino))
%subplot(1,2,2); plot(f,mean(pxx_ismr))
loglog(f,mean(pxx_nino),f,mean(pxx_ismr)*10e9,'LineWidth',1)
legend({'Nino3';'AISMR'})
xlabel('f [yr^-1]')
ylabel('power [arbitrary unit]')
set(findall(gcf,'-property','FontSize'),'FontSize',16)