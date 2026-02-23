%[text] # Linear Inverse Model - Bering Sea Ice
%[text] This live script is the principal workflow for Cox and Penland (2026).
%[text] **Objective:** Build a linear inverse model (LIM) to evaluate the predictibility of wind-forced expansion and retreat of the wintertime Bering Sea ice edge.
%[text] **Purpose:** Update what is known about predictablity of this system; benchmark expectations for minimum skill of NOAA's fully-coupled UFS-Arctic system, currently in development; potential to develop into future statistical forecasst guidance tool
%[text] **Timescales of interest:** Synoptic to subseasonal scales, mostly 48 hours to 3 weeks.
%[text] **Code Dependencies:** 
%[text] MATLAB R2025b, but should be backward compatable to ~R2019b. 
%[text] Torrence and Compo (1998) wavelet codes [https://github.com/ct6502/wavelets](https://github.com/ct6502/wavelets) 
%[text] Colorbrewer for colormaps. The exact code is depreciated, I think, so I suggest replacing cbrewer.m calls with your preference.
%[text] There is a section near the end that uses m\_map, [https://www-old.eoas.ubc.ca/~rich/map.html](https://www-old.eoas.ubc.ca/~rich/map.html) and another that uses frank-pk-DataViz-3/daboxplot, [https://www.mathworks.com/matlabcentral/fileexchange/74851-daboxplot](https://www.mathworks.com/matlabcentral/fileexchange/74851-daboxplot) 
%[text] **References:**
%[text] Cox, C.J. and C. Penland (2026) Medium-range predictability of the wintertime Bering Sea ice edge using Linear Inverse Modeling. submitted to *Journal of Geophysical Research - Oceans*.
%[text] Torrence, C. and G.P. Compo (1998) A practical guide to wavelets. Bulletin of the American Meteorological Society, 79, 61-78. [https://doi.org/10.1175/1520-0477(1998)079\<0061:APGTWA\>2.0.CO;2](https://doi.org/10.1175/1520-0477(1998)079%3C0061:APGTWA%3E2.0.CO;2)
%[text] \_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_
%%
%[text] ## I. Options and Loading
%[text] The primary data is the 25x25 km sea ice index based on SSM/I&SSMIS satellite sensors produced by the National Snow and Ice Data Center (NSIDC). We use Sea Ice Index (SII) v3 (Fetterer et al., 2017). SII begins 26 October 1978, but it is provided only every 2 days until 13 January 1988. Here we will select the period of continuous daily granularity.
%[text] LIM "stations" for the present model are meridians of longitude. Valid selections are 180°E  to 200°E. We have preselected 181°E  to 198°E for reasons described in the manuscript, though the we show the work ahead of time.
%[text] The Aleutian Low Beaufort Sea Anticyclone (ALBSA) is a daily index of 850 hPa geopotential height that is sensitive to the juxtaposition of the Aleutian Low and the Beaufort High. It captures the regional atmospheric dynamical patterns we hypothesize to be important to the tendencies of the wind forcing on the sea ice. See Cox et al. (2019).
%[text] **New References:**
%[text] Cox, C. J., Stone, R. S., Douglas, D. C., Stanitski, D. M., and Gallagher, M. R. (2019) The Aleutian Low – Beaufort Sea Anticyclone: A climate index correlated with the timing of springtime melt in the Pacific Arctic. *Geophysical Research Letters*, 46(13), 7464-7473. [https://doi.org/10.1029/2019GL083306](https://doi.org/10.1029/2019GL083306) 
%[text] Cox, C. J. (2025). Aleutian Low – Beaufort Sea Anticyclone (ALBSA) index, 1940-2025. \[Dataset\]. Arctic Data Center. [https://doi.org/10.18739/A2WP9T85G](https://doi.org/10.18739/A2WP9T85G)
%[text] Fetterer, F. Knowles, K., Meier, W. N., Savoie, M., and Windnagel, A. K. (2017) Sea Ice Index (G02135, Version 3). \[Dataset\]. National Snow and Ice Data Center. [https://doi.org/10.7265/N5K072F8](https://doi.org/10.7265/N5K072F8) 
% Set you paths
code_path  = '/path' ; % access to project mfiles
fig_path   = '/path' ; % store figures here
data_path  = '/path' ; % NSIDC matfile data directory 
albsa_path = '/path' ; % Location of ALBSA data
nsidc_path = '/path' ; % NSIDC native data directory

% Load NSDIC?
%   1 = Reload from file (takes a while)
%   2 = Load preprocessed matfile
load_type = 2;

% Save Figures?
%   1 = Save 
%   2 = Skip
sav_fig = 2;

% Save Model Weights?
%   1 = Save
%   2 = Skip
save_weights = 2;

% Set meridians [degrees longitude]
stas = 175:200; 

% Begin date of NSIDC data.
date_beg = datenum([1988 1 13 0 0 0]);
date_end = datenum([2025 7 1 0 0 0]); % make this 90 days after March 31 in last year

% De-seasonalizing options
% deseas 
%   = 1 remove 1st component
%   = 2 remove 1st, 2nd components
%   = 3 remove 1st, 2nd, 3rd  components
deseas = 2;

% What months are included in the analysis?
mo_inc = 1:3;

% Data Split, training and validation sets
train_years = 1988:2019;
validation_years = 2020:2025;

% Access some stuff
addpath(code_path); 
load([data_path,'nsidc_data']);

% % Pad the time series with NaN to give longer forecasts somewhere to be
add_dates = max(iceedge_dn)+1:date_end;
iceedge_dn = [iceedge_dn add_dates];
iceedge = transpose(padarray(transpose(iceedge),length(add_dates),NaN, 'post'));

% Subset period of analysis
ibeg = find(iceedge_dn == date_beg);
iend = find(iceedge_dn == date_end);
iceedge = iceedge(:,ibeg:iend);
iceedge_dn = iceedge_dn(ibeg:iend);

% Date info
[yyyy,mm,dd,~,~,~] = datevec(iceedge_dn); % year, month, day
iceedge_doy = iceedge_dn-datenum([yyyy;yyyy.*0+1;yyyy.*0;yyyy.*0;yyyy.*0;yyyy.*0]')'; % day of year

% Collect indexes for training and validation sets
train_inds = find(ismember(yyyy,train_years));
validation_inds = find(ismember(yyyy,validation_years));

% Loads daily ALBSA index values (Cox et al., 2019)
% - Resamples to the NSDIC time series, removes April-December.
albsa = rd_netcdf([albsa_path,'albsa_index.nc']);
albsa.dn = datenum([1940 1 1 0 0 0]) + double(albsa.valid_time);
i_resamp = find(ismember(albsa.dn,iceedge_dn));
albsalim = albsa.index(i_resamp); albsalim(~ismember(mm,mo_inc)) = NaN;
amap = albsa.z(:,:,i_resamp); amap(:,:,~ismember(mm,mo_inc)) = NaN;
era_uv = rd_netcdf([data_path,'winds4albsa.nc']);
era_uv.dn = datenum([1988 1 1 0 0 0]) + double(era_uv.valid_time);
i_resamp = find(ismember(era_uv.dn,iceedge_dn));
era_uv.u = era_uv.u(:,:,i_resamp);
era_uv.v = era_uv.v(:,:,i_resamp);

% Number of "stations", which for this LIM are longitudes. This is ALL possible station (170-200)
nDato = size(iceedge,1);
%%
%[text] ## II. De-seasonalizing and outlier removal
%[text] This analysis will focus on variability at timescales smaller than the seasonal cycle, as well as contributions from interannual variablity in the first three components of the seasonal cycle. At the pan-Arctic scale, Lopes et al. (2023) showed that the first three components of the seasonal cycle in the sea ice extent are all signficant, representing 95.6%, 2.69%, and 0.43% of the variance. The remaining 1.3% was considered residual "noise". We will explore the predctability of this residual given the hypothesis that persistence from the polar jet (as expressed by ALBSA) results in mean tendencies for on- or off-ice floe a 1-2 weeks (Cox et al. 2019). We will also repeat the analysis being inclusive of component three, and again with components two and three. Therefore, we define the seasonal cycle continously, not climatologically (it varies interannually). We use a wavelet deseasonalizing method (e.g., Szolgayová et al., 2014). Specfically, we reconstruct the seasonal cycle time series from the transform using pre-defined thresholds defined by the seasonal harmonic periods and adjusted for the Fourier period of the wavelet.
%[text] The wavelet procedure is Torrence and Compo (1998). Output printed to screen is produced by wavelet.m, called by the waveletwrapper. The output refers to closure equations that firstly provide a check that a sufficiently large po2 was selected for reconstruction (of the full time series) and secondly that the reconstruction accuracy is acceptable. Past work has suggested that reconstructions should have losses \< 2% (Cox et al., 2014), which we achieve here.
%[text] **New References:**
%[text] Cox, C.J., V.P. Walden, G.P. Compo, P.M. Rowe, M.D. Shupe, and K. Steffen (2014) Downwelling longwave flux over Summit, Greenland, 2010-2012: Analysis of surface-based observations and evaluation of ERA-Interim using wavelets. *Journal of Geophysical Research - Atmospheres*, 119, 12317-12337. [https://doi.org/10.1002/2014JD021975](https://doi.org/10.1002/2014JD021975) 
%[text] Lopes, F., Couterillot, V., Gibert, D., and Le Mouël, J.-L. (2023) On the annual and semi-annual components of variations in extent of Arctic and Antarctic sea-ice. *Geosciences*, 13, 21. [https://doi.org/10.3390/geosciences13010021](https://doi.org/10.3390/geosciences13010021) 
%[text] Szolgayová, E., Arlt, J., Blöschl, G., and Szolgay, J. (2014) Wavelet based deseasonalization for modelling and forecasting of daily discharge series considering long range dependence. *Journal of Hydrology and Hydromechanics*, 62, 24-32. [https://doi.org/10.2478/johh-2014-0011](https://doi.org/10.2478/johh-2014-0011) 
iceedgedt = iceedge; % dt is the de-seasonalized time series
seas = []; % store the season time series
ndspk = 0; % de-spiking counter
gooddata = NaN(1,nDato);
for k = 1:nDato
    
    fprintf(['Deseasonalizing Meridian = ',num2str(stas(k)),'\n'])

    % Assign vector, interpolating over nans
    merid_k = iceedge(k,:);
    merid_k_gaps = interp1( iceedge_dn(~isnan(merid_k)), merid_k(~isnan(merid_k)) , iceedge_dn, 'linear'); % linear interp gaps
    merid_k_nonan = interp1( iceedge_dn(~isnan(merid_k_gaps)), merid_k_gaps(~isnan(merid_k_gaps)) , iceedge_dn, 'nearest','extrap'); % remainder is padding, so nearest

    % first despike by comparison to a 10-day running mean
    y1 = movmean(merid_k_nonan,10);
    tmp = iceedge(k,:)-y1;
    ibad = find(abs(tmp) > std(tmp)*8);
    ndspk = ndspk + length(ibad);
    iceedge(k,ibad) = NaN;
    gooddata(k) = length(find(~isnan(iceedge(k,ismember(mm,mo_inc)))));
    
    % fill gaps. this is just for fitting the seasonal cycle
    % have tried also gap filling with neighboring stations and climatology, but this has little effect
    merid_k = iceedge(k,:);
    merid_k_gaps = interp1( iceedge_dn(~isnan(merid_k)), merid_k(~isnan(merid_k)) , iceedge_dn, 'linear'); % linear interp gaps
    merid_k_nonan = interp1( iceedge_dn(~isnan(merid_k_gaps)), merid_k_gaps(~isnan(merid_k_gaps)) , iceedge_dn, 'nearest','extrap'); % remainder is padding, so nearest

    % then deseasonalize
    data = waveletwrapper(merid_k_nonan,iceedge_dn);

    % the equivalent Fourier period of the wavelet (Table 1, T&C1998)
    lambda = (4*pi*data.period) / (6 + sqrt(2+6^2)); 

    % now filter:
    if deseas == 1
        ind = find(lambda > exp((log(365/1)+log(365/2))/2)); % remove 1st component
    elseif deseas == 2
        ind = find(lambda > exp((log(365/2)+log(365/3))/2)); % remove 1st, 2nd components
    elseif deseas == 3
        ind = find(lambda > exp((log(365/3)+log(365/4))/2)); % remove 1st, 2nd, 3rd  components
    end
   
    y1 = sum(data.x(:,ind),2) + mean(merid_k_nonan);
    seas = [seas; y1'];
    iceedgedt(k,:) = iceedge(k,:)-y1';
end
% Display the variance of the original timeseries, seasonal cycle
disp(['Variance of raw timeseries: ',num2str(var(mean(iceedge,'omitnan'),'omitnan')),' deg^2']);
disp(['Variance of seasonal cycle: ',num2str(var(mean(seas,'omitnan'),'omitnan')),' deg^2']);

% We are interested in wintertime data, though the seasonal cycle is, in
% practice, most easily separable for JFM. Therefore, we set all other
% months, April-December, to NaN.
iceedgedt(:,~ismember(mm,mo_inc))=NaN;

% despike again, this time by year
yyyyu = unique(yyyy);
for k = 1:length(yyyyu)
    ind = find(yyyy == yyyyu(k));
    a1 = [iceedgedt(:,ind)'-mean(iceedgedt(:,ind),'omitnan')']';
    a2 = reshape(a1,size(a1,1)*size(a1,2),1);
    ind2 = find(abs(a2) > std(a2,'omitnan').*5);
    ndspk = ndspk + length(ind2);
    tmp = iceedgedt(:,ind);
    tmp(ind2) = NaN;
    iceedgedt(:,ind) = tmp;
end

disp(['Despiking removed ', num2str(ndspk ./ (length(find(~isnan(iceedgedt)))+ndspk) * 100),'%']);
%%
%[text] ## III. Subset stations
%[text] So far, we have been inclusive of all meridians between 175 and 200 °E, as were selected by Wohl (1991). However, meridians near the edges of the domain have less data and less variance. They may not be very representative of the dynamics and could artificially inflate the calculated skill. This section shows how we came up with the subset 181 and 198 °E.
%[text] **New References:**
%[text] Wohl, G.M. (1991) Sea ice edge forecast verifcation program for the Bering Sea. *National Weather Digest*, 16(4), 6-12.
% Calculate the variance and "uptime" for each station
nlags = 90;
sta_uptime = NaN(nDato,90)       ; % fraction of time the station has valid data, JFM
sta_var    = NaN(nDato,1)        ; % variance of the station distribution, JFM
sta_af     = NaN(nDato,nlags+1)  ; % station autocorrelation function, JFM
sta_eF     = NaN(nDato,1)        ; % station eFolding time based on sta_af
for h = 1:nDato

    % calculate variance of the station
    sta_var(h) = var(iceedgedt(h,:),'omitnan');

    % calculate the autocorrelation function, save it
    % calculate the eFolding time from the AF
    datatmp = iceedgedt(h,:);
    output = [0:nlags].*NaN;
    for q = 0:nlags
        if q == 0
            output(h+1) = 1;
        else
            x = transpose(datatmp(1:end-q)) - mean(datatmp,'omitnan');
            y = transpose(datatmp(q+1:end)) - mean(datatmp,'omitnan');
            [xx pp] = nancorrcoef(x,y);
            if length(xx) > 1
                output(q+1) = xx(2);
            end
        end
    end
    output(1)=1;
    sta_af(h,:) = output; 
    nl2 = 0:0.01:nlags;
    highres = interp1(0:nlags,output,nl2);
    sta_eF(h) = nl2(find(highres < 1/exp(1),1,'first'));
    
    % calculate uptime
    for k = 1:90
        ind1 = find(iceedge_doy == k);
        sta_uptime(h,k) = length(find(~isnan(iceedgedt(h,ind1)))) / length(ind1);
    end

end

% Criteria for subsetting stations: 181-198 
%stind = find(mean(sta_uptime,2,'omitnan') >= 0.75 & sta_eF <= 16 & sta_var >= 0.32); % 0.32 is ~35 km using deseas = 2
stind = 7:24; % hard-code to facilitate deseas

% Save the unconsolidated values for Figure 2
ie = iceedgedt;

% Now with all the statistics compiled, subset the data and send to LIM
iceedgedt = iceedgedt(stind,:);
nDat = size(iceedgedt,1);
%%
%[text] ## IV. Figure 1, Map
%[text] Map of the Bering Sea region (Albers Equal-Area Conic projection, basemap made with Natural Earth). Sea ice concentrations (SIC) are an example day, 14 February 2023, using NSIDC Sea Ice Index data. White dots are the identified ice edge latitudes at each of the 18 meridians used for the LIM, which are displayed as solid gray lines. Dashed gray lines are other meridians considered previously. The Bering Sea shelf break is visible in the basemap southwest of the ice edge near the center of the map.
% get the grid
fid = fopen([nsidc_path,'psn25lats_v3.dat']);
grlats = fread(fid,[304 448],'int',0,'n')./100000;
fclose(fid);
fid = fopen([nsidc_path,'psn25lons_v3.dat']);
grlons = fread(fid,[304 448],'int',0,'n')./100000;
grlons = mod(grlons,360);
fclose(fid);
lon = 0:1:360;
lat = 50:0.1:90;
latflip = fliplr(lat);
[xq,yq] = meshgrid(lon,lat);

beg_day = datenum([2023 2 14 0 0 0]);
beg_day_ind = find(iceedge_dn == beg_day);
tmp1 = importdata([nsidc_path,'2023/02_Feb/N_20230214_concentration_v3.0.tif']); 
tmp1 = double(tmp1.cdata);
tmp1(tmp1 > 1000) = NaN; tmp1 = tmp1./10;
sic = griddata(grlons,grlats,double(tmp1'),double(xq),double(yq))./10;

latlim = [50 72];
lonlim = [175 205];

[A,RA,attrib] = readBasemapImage("landcover",latlim,lonlim,5);

clf;
[xGrid,yGrid] = worldGrid(RA);
[latGrid,lonGrid] = projinv(RA.ProjectedCRS,xGrid,yGrid);
lonGrid=wrapTo360(lonGrid);
h = worldmap(latlim,lonlim);
setm(h,'MapProjection','eqaconicstd')
geoshow(latGrid,lonGrid,A)
tmp2 = tmp1; tmp2(isnan(tmp1)) = 0;
tmp3 = tmp2; tmp3(tmp3 > 0) = 1;
e1 = pcolorm(grlats,grlons,tmp2'); shading flat
set(e1, 'FaceAlpha', 'flat', 'AlphaData', tmp3','edgecolor','none');
hh=colorbar; hh.Label.String = 'SIC [%]'; clim([0 100]); hh.FontSize = 8;
hh.Position(4) = 0.6; hh.Position(2) = 0.25; hh.Position(1) = 0.8;
lon1 = stas(stind);
lat1 = iceedge(lon1-(min(stas)-1),beg_day_ind);
for k = 1:length(lon1)
    plotm(lat1(k),lon1(k),'marker','.','markersize',12,'color',[0.99 0.99 0.99]);
end
for k = stas
    plotm([50 72],[k k],'k--','color',[0.5 0.5 0.5])
end
for k = 1:length(stind)
    plotm([50 72],[stas(stind(k)) stas(stind(k))],'color',[0.25 0.25 0.25])
end
setm(h,'mlabellocation',5,'plabellocation',5,'grid','off','FontSize',8);
mlabelzero22pi; % not obvious, but switches label from degW to degE 
colormap(parula);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 18 14]);
set(gca,'fontsize',14); if sav_fig == 1 & deseas == 2; print('-dpng','-r300',[fig_path,'fig1.png']); end
%%
%[text] ## V. Figure 2a,b,c,d; Figure S2
%[text] *Figure 2:* 
%[text] (a) Autocorrelation functions for each meridian. Solid lines are merdians included in the LIM and dashed lines are those at the edge of the domain included by Wohl (1991), but not considered here; (b) Contemporaneous cross-correlations amongst the meridians; (c) lagged (𝜏 = 6 days) cross-correlations amongst the meridians; (d) difference between (b) and (c) (lagged – contemporaneous). The boxes in panels (b)-(d) surround the meridians included in the LIM.
%[text] *Figure S2:*
%[text] Frequency distribution functions of xi anomaly time series. All were deseasonalized by removing the first and second components of the seasonal cycle. Bin size is 0.1°.
% Calculate matrices

lag = 6; % use a 6-day lag

% Figure 2b: Contemporaneous correlation matrixm, Rc
%   slopes code calculates how quickly the covariance decreases from with increasing distance from a meridian
Rc = NaN(nDato,nDato);
Rcslopes1 = NaN(nDato,nDato);
Rcslopes2 = NaN(nDato,nDato);
for k = 1:nDato
    for h = 1:nDato
        [r p] = corrcoef(ie(h,:)',ie(k,:)','Rows','pairwise');
        Rc(h,k) = r(2);
        if h >= k
            Rcslopes1(h,k) = r(2)^2;
        end
        if k > h
            Rcslopes2(h,k) = r(2)^2;
        end
    end
end
slopes = NaN(2,nDato);
for k = 1:nDato
    thex = stas;
    they = Rcslopes1(k,:);
    thex(isnan(they)) = [];
    they(isnan(they)) = [];
    if length(thex) >= 11
        tmp = polyfit(thex,they,1);
        slopes(1,k) = tmp(1);
    end
    thex = stas;
    they = Rcslopes2(k,:);
    thex(isnan(they)) = [];
    they(isnan(they)) = [];
    if length(thex) >= 11
        tmp = polyfit(thex,they,1);
        slopes(2,k) = tmp(1);
    end
end

% Figure 2b: Lagged correlation matrix
R = NaN(nDato,nDato);
for k = 1:nDato
    for h = 1:nDato
        [r p] = corrcoef(ie(h,1:end-(lag-1))',ie(k,lag:end)','Rows','pairwise');
        R(h,k) = r(2);
    end
end

% Create a figure of the autocorrelation functions for each station.
% Save the eFolding times, we need them later for the AR1 fx.
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
tiledlayout(2, 2, "TileSpacing", "compact");
cmapax1 = jet(nDato);
cmapax4 = flipud(cbrewer('div', 'RdBu',nDato, 'linear'));

ax1=nexttile;
for k = 1:nDato
    if ismember(stas(k),stas(stind))
        plot(0:nlags,sta_af(k,:),'color',cmapax1(k,:),'linewidth',2); hold on;
    else
        plot(0:nlags,sta_af(k,:),'k--','color',cmapax1(k,:),'linewidth',1); hold on;
    end
end
xlabel('\tau [days]');
ylabel('autocorr [r]');
grid on;
xlim([0 60]); ylim([-0.6 1]);
lats = stas;
set(gca,'ytick',-0.6:0.1:1,'xtick',0:10:100);
%plot([0:60],exp((-1/mean(sta_eF(stind)))*[0:60]),'k','linewidth',2)
h1=legend(num2str(transpose(lats)),'location','northeast'); 
h1.FontSize = 12;
h1.NumColumns = 4;
pa=text(-9,1,'a)'); pa.FontSize=20; pa.FontWeight = 'bold';
set(gca,'fontsize',14);

ax2=nexttile;
imagesc(stas,stas,Rc); set(gca,'ydir','normal');
xlabel('Longitude [\circ]'); ylabel('Longitude [\circ]'); 
clim([0 1]); h=colorbar; h.Label.String = 'r';
colormap(parula);
rectangle('Position',[stas(stind(1))-0.5 stas(stind(1))-0.5 stas(stind(end))-stas(stind(1))+1 stas(stind(end))-stas(stind(1))+1],'linewidth',2)
set(gca,'xtick',174:2:200,'ytick',174:2:200);
colormap(ax2,'parula');
pb=text(170,200,'b)'); pb.FontSize=20; pb.FontWeight = 'bold';
set(gca,'fontsize',14);

ax3=nexttile;
imagesc(stas,stas,R); set(gca,'ydir','normal');
xlabel('Longitude [\circ]'); ylabel('Longitude [\circ]'); 
clim([0 1]); h=colorbar; h.Label.String = 'r';
colormap(parula);
rectangle('Position',[stas(stind(1))-0.5 stas(stind(1))-0.5 stas(stind(end))-stas(stind(1))+1 stas(stind(end))-stas(stind(1))+1],'linewidth',2)
set(gca,'xtick',174:2:200,'ytick',174:2:200);
colormap(ax3,'parula');
pc=text(170,200,'c)'); pc.FontSize=20; pc.FontWeight = 'bold';
set(gca,'fontsize',14);

ax4=nexttile;
imagesc(stas,stas,R-Rc); set(gca,'ydir','normal');
xlabel('Longitude [\circ]'); ylabel('Longitude [\circ]'); 
clim([-0.4 0.4]); h=colorbar; h.Label.String = 'r';
rectangle('Position',[stas(stind(1))-0.5 stas(stind(1))-0.5 stas(stind(end))-stas(stind(1))+1 stas(stind(end))-stas(stind(1))+1],'linewidth',2)
set(gca,'xtick',174:2:200,'ytick',174:2:200);
colormap(ax4,cmapax4);
pd=text(170,200,'d)'); pd.FontSize=20; pd.FontWeight = 'bold';
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 40 32]);
set(gca,'fontsize',14); if sav_fig == 1 & deseas == 2; print('-dpng','-r300',[fig_path,'fig2_raw.png']); end

% Make the histogram supplements plot
% % See distributions of each station.
clf;
thevars = NaN(1,nDat);
for k = 1:nDat
    subplot(4,5,k);
    thevars(k) = std(iceedgedt(k,:));
    histN(iceedgedt(k,:),-15:0.1:15,'line','linewidth',1); hold on;
    title(([num2str(stas(stind(k))),'\circ']));
    set(gca,'xtick',-15:5:15); xlim([-6 6]); grid on;
    if k == nDat; xlabel('[\circ]'); end
    if k == 11; ylabel('                                [%]'); end
end
subplot(4,5,nDat+1);
histN(mean(iceedgedt(:,:)),-15:0.1:15,'line','linewidth',1);
set(gca,'xtick',-15:5:15); xlim([-6 6]); grid on; title('Mean');
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 32 16]);
set(gca,'fontsize',12); if sav_fig == 1 & deseas == 2; print('-dpng','-r300',[fig_path,'figS2.png']); end
%%
%[text] ## VI. Do LIM
%[text] The code block is a wrapper for dolim.m. LIMs are built for a series of proposed 𝜏0 (tau0\_arr) for leads, 𝜏 (Tau\_arr), of 1-90 days; forecasts, expected errors (espilon, ε), and diagnostics (e.g., Nyquist frequency detections) are stored in arrays. We will use a 𝜏0 = 4. The nominal criterion for selecting 𝜏0 is that it doesn't matter! Actually, there are some guidelines. The assumed lack of dependency of results on 𝜏0 will be tested in Section VII. Let's aim for the smallest acceptable values. You cannot use a 𝜏0 \>= any 𝜏 having a Nyquist frequency at stated 𝜏0 (Penland, 2019), and best perhaps to buffer that criterion too. For this application, 𝜏0 \< 3 are suspect due to the nature of the input data (see Section VII). Later, we will also discover that the observed lack of dependency of forecast skill on 𝜏0 actually corresponds to a tradeoff between the maximum amplitude of the stochastic growth and the correlation between the optimal and final structures. Nominally, this code block would choose the largest value of 𝜏0 that is -2 𝜏 from the first detected Nyquist frequency. In practice, I find that the appearance of these Nyquists is quite sensitive to even small adjustments to feature engineering of the input data (e.g., handling of outliers). The present configuration does not find any Nyquists until rather large 𝜏 so I have simply chose 𝜏0 = 4. 
%[text] There are several credits I would like to provide here, which also apply to Section VI and VII.
%[text] (1) Penland, C., M.D. Fowler, D.L. Jackson, and R. Cifelli (2021) Forecasts of opportunity in Northern California soil moisture. *Land*, 10, 713. [https://doi.org/10.3390/land10070713](https://doi.org/10.3390/land10070713) 
%[text] (2) Megan Fowler's LIM routines developed for Penland et al. (2021), [https://github.com/NOAA-PSL/Linear\_Inverse\_Modeling](https://github.com/NOAA-PSL/Linear_Inverse_Modeling) 
%[text] (3) Cecile Penland's *Introduction to Stochastic Processes & Linear Modeling* coursework, originally developed as a guest lecturer series at MIT.
%[text] **New References:**
%[text] Penland, C. (2019) The Nyquist issue in Linear Inverse Modeling. Monthly Weather Review, 147, 1341-1349. [https://doi.org/10.1175/MWR-D-18-0104.1](https://doi.org/10.1175/MWR-D-18-0104.1)
% For shuffling test
% % Data Split, training and validation sets
% ayrs = 1988:2025;
% vinds = randi(length(ayrs),6,1);
% train_years = ayrs(~ismember(1:length(ayrs),vinds)); %1988:2019;
% validation_years = ayrs(vinds); %2020:2025;
% % Collect indexes for training and validation sets
% train_inds = find(ismember(yyyy,train_years));
% validation_inds = find(ismember(yyyy,validation_years));
%
% For fixed-origin rolling cross-validation
%train_years = 1988:1999;
%validation_years = 2000:2005;
%train_inds = find(ismember(yyyy,train_years));
%validation_inds = find(ismember(yyyy,validation_years));


% We will test the first 9 tau0s, each lags out to 90 days.
x = transpose(iceedgedt); % get oriented.
disp(['Number of training samples = ',num2str(length(find(~isnan(x(train_inds,:)))))]);
disp(['Number of validation samples = ',num2str(length(find(~isnan(x(validation_inds,:)))))]);
% fx 
tau0_arr = 1:9                                                           ; % Test these tau0s 
Tau_arr  = 1:90                                                          ; % Create array of Tau's to get error at (out to ~3 weeks)
    
fx_all            = NaN(length(tau0_arr),length(Tau_arr),nDat,length(x)) ; % Fx stored 4d array, [tau0,leadtime,station,time]
epsilon_sta_val   = NaN(length(tau0_arr),length(Tau_arr),nDat)           ; % Expected error per station
epsilon_val       = NaN(length(tau0_arr),length(Tau_arr))                ; % Expected error collectively
epsilon_sta_train = NaN(length(tau0_arr),length(Tau_arr),nDat)           ; % Expected error per station
epsilon_train     = NaN(length(tau0_arr),length(Tau_arr))                ; % Expected error collectively
epsilon_raw       = NaN(length(tau0_arr),length(Tau_arr))                ; % Expected error collectively, not normalized
compL             = zeros(size(tau0_arr))                                ; % Binary, 1 = L was complex, else 0 
nyquists          = zeros(length(tau0_arr),length(Tau_arr))              ; % Binary, 1 = Nyquist frequency found
gamma1max         = zeros(length(tau0_arr),length(Tau_arr))              ; % 
firstnyq          = zeros(size(tau0_arr))                                ; % Smallest tau0 where Nyquist was found

% Run hindcasts for a series of tau0s
for iT0 = 1:length(tau0_arr)
    
    % Do LIM for this tau0
    [L,Q,B,Gtau,C0,Ctau,ualpha,valpha,galpha,tau_decay_alpha,T_mode_oscil] = dolim(x(train_inds,:),tau0_arr(iT0));
    if ~isreal(L); compL(iT0) = 1; end
   
    % Loop over lead time array (Tau_arr)
    for iTau = 1:length(Tau_arr) 

        % Check for Nyquists
        nyqtmp = abs(imag(diag(B))*(Tau_arr(iTau)));
        if any((nyqtmp>=3.138) & (nyqtmp<=3.142))
            nyquists(iT0,iTau) = 1;
        end
 
        % Define Green function at this lag 
        G_tau = real(ualpha * diag(galpha).^(Tau_arr(iTau)/tau0_arr(iT0)) * transpose(valpha));
        eps = C0-(G_tau*(C0*transpose(G_tau)));
        for k = 1:nDat
            epsilon_sta_train(iT0,iTau,k) =  eps(k,k) / var(x(train_inds,k),'omitnan');
            epsilon_sta_val(iT0,iTau,k) =  eps(k,k) / var(x(validation_inds,k),'omitnan');
        end
        epsilon_train(iT0,iTau) = trace(eps) / sum(var(x(train_inds,:),'omitnan'));
        epsilon_val(iT0,iTau) = trace(eps) / sum(var(x(validation_inds,:),'omitnan'));
        epsilon_raw(iT0,iTau) = trace(eps) / nDat;
        
        % Calculate gamma
        [Psi_tmp,gamma_tmp] = eig(transpose(G_tau)*G_tau);
        gamma1max(iT0,iTau) = max(diag(gamma_tmp));
        
        fx = NaN(nDat,length(x));

        for iT = 1:length(x)-Tau_arr(iTau)    % Loop over timesteps 
        
            tmp = x(iT,:);
            ii = find(isnan(tmp));
            tmp(ii) = 0;
            thisfx = G_tau*transpose(tmp);
            thisfx(ii) = NaN;
            % Fx: compute forecast 
            fx(:,iT+Tau_arr(iTau)) = thisfx;
             
        end    

        fx_all(iT0,iTau,:,:) = fx;

    end

    tmp = find(nyquists(iT0,:) == 1,1,'first'); if isempty(tmp); tmp = max(Tau_arr); end
    firstnyq(iT0) = tmp;
        
end

% Pick your tau_0. We will back off two steps from the earliest Nyquist
tau_0 = find(sum(nyquists) == 1,1,'first')-2;
if isempty(tau_0); tau_0 = 4; end
tau_0 = 4; disp('You selected your own tau.');
disp(['Selected tau0 is ',num2str(tau_0)]);
% Run the code one more time using our tau_0 pick. Run on the training set.
[L,Q,B,Gtau,C0,Ctau,ualpha,valpha,galpha,tau_decay_alpha,T_mode_oscil] = dolim(x(train_inds,:),tau_0);
% Save the weights to hdf5
if save_weights == 1

    % query some info on this code release
    cd(code_path)
    [status, cmdout] = system('git describe --tags --abbrev=0');
    if status == 0
        release_version = strtrim(cmdout);
    else
        release_version = 'unknown';
    end

    % filename
    filename = [code_path,'bering_lim_model.h5'];

    % start clean
    if exist(filename,'file'); delete(filename); end

    % store output from dolim, including tau0 reference input
    h5create(filename,'/model/B',[size(B),2],'Datatype','double');
    h5write(filename,'/model/B',cat(3,real(B),imag(B)));
    h5writeatt(filename,'/model/B','description','eigenvalue of L, Beta_alpha; dimension three is real and imag parts');

    h5create(filename,'/model/C0',size(C0),'Datatype','double');
    h5write(filename,'/model/C0',C0);
    h5writeatt(filename,'/model/C0','description','contemporaneous covariance matrix');

    h5create(filename,'/model/Ctau',size(Ctau),'Datatype','double');
    h5write(filename,'/model/Ctau',Ctau);
    h5writeatt(filename,'/model/Ctau','description','lagged (by tau0) covariance matrix');

    h5create(filename,'/model/L',[size(L),2],'Datatype','double');
    h5write(filename,'/model/L',cat(3,real(L),imag(L)));
    h5writeatt(filename,'/model/L','description','linear operator, L; dimension three is real and imag parts');

    h5create(filename,'/model/Q',size(Q),'Datatype','double');
    h5write(filename,'/model/Q',Q);
    h5writeatt(filename,'/model/Q','description','noise covariance matrix, Q');
 
    h5create(filename,'/model/tau_decay_alpha',size(tau_decay_alpha),'Datatype','double');
    h5write(filename,'/model/tau_decay_alpha',tau_decay_alpha);
    h5writeatt(filename,'/model/tau_decay_alpha','description','decay time');

    h5create(filename,'/model/T_mode_oscil',size(T_mode_oscil),'Datatype','double');
    h5write(filename,'/model/T_mode_oscil',T_mode_oscil);
    h5writeatt(filename,'/model/T_mode_oscil','description','oscillation period');

    h5create(filename,'/model/galpha',[size(galpha),2],'Datatype','double');
    h5write(filename,'/model/galpha',cat(3,real(galpha),imag(galpha)));
    h5writeatt(filename,'/model/galpha','description','eigenvalues of Gtau; dimension three is real and imag parts');

    h5create(filename,'/model/ualpha',[size(ualpha),2],'Datatype','double');
    h5write(filename,'/model/ualpha',cat(3,real(ualpha),imag(ualpha)));
    h5writeatt(filename,'/model/ualpha','description','eigenvectors of Gtau; dimension three is real and imag parts');

    h5create(filename,'/model/valpha',[size(valpha),2],'Datatype','double');
    h5write(filename,'/model/valpha',cat(3,real(valpha),imag(valpha)));
    h5writeatt(filename,'/model/valpha','description','adjoints of Gtau; dimension three is real and imag parts');

    h5create(filename,'/model/tau0',size(tau_0),'Datatype','double');
    h5write(filename,'/model/tau0',tau_0);
    h5writeatt(filename,'/model/tau0','description','lag (tau0) that was used too develop the model');

    % globals
    h5writeatt(filename,'/model','description','Bering Sea LIM model weights based on training set from Cox & Penland (2026) JGR-Oceans');
    h5writeatt(filename,'/model','git_release',release_version);
    h5writeatt(filename,'/model','morder','F');
    h5writeatt(filename,'/model','sorting','modes are presorted in descending order');
    h5writeatt(filename,'/model','training_range',[num2str(train_years(1)),'-',num2str(train_years(end))]);

    % send us back where we were
    cd(fig_path)

end
%%
%[text] ## VII. Analysis of LIM diagnostics
%[text] - (1) MAC
%[text] - (2) Optimal Structure
%[text] - (3) Tau Test \
%[text] **(1) Calculate the Maximum Amplification Curve (MAC)**
%[text] The MAC is calculated corresponding to each (ranked) eigenvalue (n = nDat). We are looking for peaks some lead time 𝜏 \> 1 and amplitude \> 1. Else, the solution to the linear system does not add information over damped persistence. The figure shows the leading eigenvalue corresponds to a significant MAC, peaking at 𝜏 = 5. 
gamma1max = [];
for j = 1:max(Tau_arr)
    Gt_taum = real(ualpha * diag(galpha).^(j/tau_0) * transpose(valpha));
    [Psi1,gamma1] = eig(transpose(Gt_taum)*Gt_taum); %[EIGVEC,EIGVAL]
    gamma1 = diag(gamma1);
    if any(imag(gamma1)); break; end
    [ai,ab] = sort(gamma1,'descend');
    gamma1max = [gamma1max gamma1(ab)];
end

tau_peak = [];
peak_amp = [];
legs = [];
clf;
cmap=colormap(jet(nDat));
for k = 1:nDat
    ii = find(gamma1max(k,:) == max(gamma1max(k,:)));
    tau_peak = [tau_peak; ii];
    peak_amp = [peak_amp; gamma1max(k,ii)];
    legs = [legs; {['\gamma',num2str(k)]}];
    plot(gamma1max(k,:),'color',cmap(k,:),'linewidth',2); hold on;
end
tau_m = tau_peak(1);
plot(tau_peak,peak_amp,'ko-','linewidth',2,'markersize',8)
for ii = [1 2 3:3:nDat]% 1:length(tau_peak)
    if ii == 1
        text(tau_peak(ii),peak_amp(ii)-0.2,['\gamma',num2str(ii),' at \tau_m = ',num2str(tau_m)],'Color','k','fontsize',14)
    elseif ii == 2
        text(tau_peak(ii)+0.4,peak_amp(ii)+0.1,['\gamma',num2str(ii)],'Color','k','fontsize',14)
    else
        text(tau_peak(ii)-1.25,peak_amp(ii),['\gamma',num2str(ii)],'Color','k','fontsize',14)
    end
end
xlabel('\tau'); ylabel('amplification'); hline(1);
xlim([-0.5 20]); grid on;
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 24 18]);
set(gca,'fontsize',14); if sav_fig == 1 & deseas == 2; print('-dpng','-r300',[fig_path,'fig5.png']); end 
%[text] **(2) Analyze the optimal structure**
%[text] The leading eigenvalue has a peak \> 𝜏 = 1 at and amplitude \> 1. All other eigenvalues produce MACs that ~ \<= 1 peak ampltide and generally fall off sharply after 𝜏 =1. Thus, there is only one optimal structure. This should be like saying there is a particular orientation or shape of the ice edge that is associated with high predictability due to stochasitc growth. Initial indications were that this appears when the profile of the ice edge has an anomalously meridional component to it, specifically that the ice is farther north in the west and farther south in the east. Subsequent analysis indicated that the first conclusion, i.e., that only the leading eigenvalue is significant, is robust, but that the second, i.e., a more meridional initial condition produces a more confident forecast, was not. Specifically, we determined that the source of the orientation-based skill was due to the 2nd and/or 3rd components of the seasonal cycle, possibly infuenced by the treatment of the time series to manage these sources of variability.
%[text] Here we plot the loading patterns. We will display values for 𝜏 = 3 to 17 to be able to assess the consistency of the structure with 𝜏, but recall that the optimal structure was determined to be at 𝜏 = 6 from the previous plot. 
cmap2 = flipud(cbrewer('seq', 'OrRd',17, 'linear'));
cmap = [0.000 0.447 0.741; ... % Blue
        0.850 0.325 0.098; ... % Orange
        0.929 0.694 0.125; ... % Yellow
        0.494 0.184 0.556; ... % Purple
        0.466 0.674 0.188; ... % Green
        0.301 0.745 0.933; ... % Light Blue
        0.635 0.078 0.184];   % Dark Red
h = [];
legs = [];
clf;
for k = 3:17
    tau_m_tmp = k;
    Gt_taum_tmp = real(ualpha * diag(galpha).^(tau_m_tmp/tau_0) * transpose(valpha));
    [Psi1,gamma1] = eig(transpose(Gt_taum_tmp)*Gt_taum_tmp); %[EIGVEC,EIGVAL]
    ind = find(diag(gamma1) == max(diag(gamma1)));

    optimStructtmp = real(Psi1(:,ind));
    if sign(optimStructtmp(1)) == -1
        %disp(['Optimal structure multiplied by -1 for ', 964,'m = ', num2str(k)]); % the sign is arbitrary. this step is to facilitate a cleaner figure.
        optimStructtmp = optimStructtmp.*-1;
    end

    htmp = plot(optimStructtmp,'-','color',cmap2(k-2,:),'linewidth',0.5,'markersize',10); hold on;
    h = [h; htmp];
    legs = [legs; {['\tau=',num2str(k)]}];

    if k == tau_m
        optimStruct = optimStructtmp;
        Gt_taum = Gt_taum_tmp;
    end
    
end

htmp = plot(optimStruct,'-*','color',cmap2(tau_m,:),'linewidth',2,'markersize',10); hold on; 
h = [h; htmp]; legs = [legs; {['Optimal Structure \tau_m=',num2str(tau_m)]}];
htmp = plot(Gt_taum*optimStruct,'-^','linewidth',2,'markersize',10,'color',cmap(1,:));
h = [h; htmp]; legs = [legs; {['Final Structure']}];
xlim([0.5 length(stas(stind))+0.5]); ylim([-0.65 0.65]); grid on;
xl = stas(stind);
ylabel('Loading'); xlabel('station [\circE]'); set(gca,'xtick',1:length(stas(stind)),'xticklabel',xl(1:length(stas(stind))));
h1=legend(h,legs,'location','southeast');
h1.FontSize = 14;
h1.NumColumns = 4;
xlim([0.5 length(stas(stind))+0.5]); ylabel('Maximum Amplitude'); set(gca,'xtick',1:length(stas(stind))); %grid on;
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 24 18]);
set(gca,'fontsize',14); if sav_fig == 1 & deseas == 2; print('-dpng','-r300',[fig_path,'fig6.png']); end
%[text] **(3) Conduct the "tau test"**
%[text] Figure plots the theoretical (expected) error, ε, for a series of 𝜏0 values. We will color code based on whether the 𝜏0 includes Nyquists at 𝜏 \< 𝜏0 (which invalidates it) in red (if occurs) and we will plot the 𝜏0 used for the model in green. Dashed lines refer to 𝜏0 \<= 3. The test is subjective: expected behavior is that ε is insensitive to the choice of 𝜏0, but the tolerance is flexible.
%[text] The test shows agreement for 𝜏0 \> 2, but we consistently find some deviation for 𝜏0 \< 3. I speculate that the lower errors for 𝜏0 \<= 3 days is associated with artificial smoothing of short-term variability by a combination of the coarseness of the gridded product, which is exacerbated by the regriddiing procedure, as well as the nature of the underlying data, which are temporal aggregates of a limted number of overpasses. It is probably best to avoid 𝜏0 \<= 3, though the errors incurred are not expected to be large.
clf;
plot([0 Tau_arr],[transpose(zeros(size(epsilon_train(compL==0 & find(compL>-1) <= 3,1))));transpose(epsilon_train(compL==0 & find(compL>-1) <= 3,:))],'b--','linewidth',1); hold on;
h1=plot([0 Tau_arr],[transpose(zeros(size(epsilon_train(compL==0 & find(compL>-1) > 3,1))));transpose(epsilon_train(compL==0 & find(compL>-1) > 3,:))],'b','linewidth',1);
try
    plot([0 Tau_arr],[transpose(zeros(size(epsilon_train(compL==1 & find(compL>-1) <= 3,1))));transpose(epsilon_train(compL==1 & find(compL>-1) <= 3,:))],'r--','linewidth',1); 
end
try
    h2=plot([0 Tau_arr],[transpose(zeros(size(epsilon_train(compL==1 & find(compL>-1) > 3,1))));transpose(epsilon_train(compL==1 & find(compL>-1) > 3,:))],'r','linewidth',1);
end
h3=plot([0 Tau_arr],[transpose(zeros(size(epsilon_train(tau_0,1))));transpose(epsilon_train(tau_0,:))],'g','linewidth',3);
ylim([0 1]);
xlabel('\tau'); ylabel('\epsilon^2 (\tau,\tau_0)'); title(['Tau test for \tau ',num2str(min(tau0_arr)),'-',num2str(max(tau0_arr))]); grid on;
try
h=legend([h1(1) h2(1) h3],'\tau_0 (L is real; no Nyquists)','\tau_0 (L is complex; has Nyquists)',['Chosen \tau_0 = ',num2str(tau_0)],'location','southeast');
catch
h=legend([h1(1) h3],'\tau_0 (L εs real; no Nyquists)',['Chosen \tau_0 = ',num2str(tau_0)],'location','southeast');
end
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 24 18]);
set(gca,'fontsize',16); grid on; if sav_fig == 1 & deseas == 2; print('-dpng','-r300',[fig_path,'figS1.png']); end
%%
%[text] ## VIII. Forecasts of opportunity
%[text] Forecasts of opportunity comprise initial conditions (IC) that are associated with particularly skillful forecasts. It's a distribution so a threshold is needed to classify skillful ICs. Penland et al. (2021) recommends the top tercile. There are two ways to idenfity the values using the LIM data:
%[text] Optimal Projection: Recall that the optimal structure is the loading pattern representing the optimal IC profile for the iceedge w.r.t. predictablity. While we found that the most salient features were linked to the third component of the seasonal cycle, the optimal loading pattern can still be projected onto the ICs to isolate optimal times.
%[text] - Use upper tercile for defining. \
% Calculate the optimal projection from the optimal stucture in the previous section
optProj = NaN(length(x),1);
optProjF = NaN(length(x),1); 
for k = 1:length(x)-tau_m
    optProj(k) = sum(optimStruct'.*x(k,:),'omitnan'); 
    optProjF(k) = sum((Gt_taum*optimStruct)'.*x(k+tau_m,:),'omitnan'); 
end
optProj(~ismember(mm,mo_inc)) = NaN;
optProjF(~ismember(mm,mo_inc)) = NaN;

clf;
histogram(abs(optProj),0:0.1:10,'normalization','percentage');
vline(prctile(abs(optProj),66.66));
xlim([0 8]); ylim([0 1.2]); grid on;
xlabel('|\Psi^T_1{\itx}(0)|');
ylabel('[%]'); grid on;
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 24 18]);
set(gca,'fontsize',16);  if sav_fig == 1 & deseas == 2; print('-dpng','-r300',[fig_path,'fig_extra.png']); end
corOS  = NaN(length(x),1);
corFS  = NaN(length(x),1);
for k = 1:length(x)-(tau_m+1)

    tmp = nancorrcoef(x(k,:),optimStruct);
    corOS(k) = tmp(2);
    tmp = nancorrcoef(x(k+tau_m,:),Gt_taum*optimStruct);
    corFS(k) = tmp(2); 

end

clf;
cmap = (cbrewer('div', 'RdBu',50, 'linear'));
scatter(corOS,corFS,20,[0;movmean(diff(mean(x,2,'omitnan')),14)],'filled','markeredgecolor',[0.75 0.75 0.75]); 
colormap(cmap);     
xlabel([{'Optimal Structure'};{'Cor(|\Psi_1^Tx(0)|,x(0))'}]);
ylabel([{'Final Structure'};{'Cor(|G(\tau=5)\Psi_1^Tx(\tau=5)|,x(\tau=5))'}]);
grid on; xlim([-0.75 0.75]); set(gca,'xtick',-0.8:0.2:0.8);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 24 14]);
box on; clim([-0.25 0.25]);
[cr,cp] = corrcoef(corOS,corFS,'rows','pairwise');
set(gca,'fontsize',14); if sav_fig == 1 & deseas == 2; print('-dpng','-r300',[fig_path,'fig8_raw.png']); end; disp(['Correlation, r = ',num2str(cr(2))]); 
% These vectors were created manually by setting tau_0 at the end of Section V in turn as in tau0_arr then running Sections VI and VII to find tau_m, peak_amp(1), and corrcoef(corOS,corFS)
% Calculated with removal of 1st, 2nd cycles (deseas = 2)

% % Version all data (x=x)
% tau0_arr = 3:12; % specified tau_0 
% tau_m_arr = [5,6,6,7,7,5,4,5,6,4]; % tau_m
% peak_amp_arr = [1.5517,2.01782,2.4706,2.7138,3.1344,3.9861,5.2457,5.85,6.6290,7.9550]; %  % maximum amplitude, first eigenval (ie @ tau_m), peak_amp(1)
% corcor_arr = [0.6568,0.5858,0.5732,0.5070,0.4657,0.4404,0.4096,0.3306,0.3082,0.2597]; % the correlation from the plot above
% first_nyq_ind = 8; % index (in tau0_arr) where the first Nyquist appears
% 
% % Version training set only, 1988-2017
% tau0_arr = 3:12; % specified tau_0 
% tau_m_arr = [6,5,5,4,5,3,3,7,6,5]; % tau_m
% peak_amp_arr = [1.6049,1.8340,2.0413,2.1411,2.3523,4.4181,5.9363,4.1702,5.1582,5.2576]; %  % maximum amplitude, first eigenval (ie @ tau_m), peak_amp(1)
% corcor_arr = [0.61758,0.64725,0.65027,0.66981,0.58012,0.55598,0.52779,0.36933,0.32847,0.34729]; % the correlation from the plot above
% first_nyq_ind = 9; % index (in tau0_arr) where the first Nyquist appears

% Version training set only, 1988-2019
tau0_arr = 3:12; % specified tau_0 
tau_m_arr = [5,5,5,4,6,4,5,3,3,3]; % tau_m
peak_amp_arr = [1.5664,1.9700,2.3135,2.5262,2.8500,4.0772,5.0162,5.4257,13.5852,10.7432]; %  % maximum amplitude, first eigenval (ie @ tau_m), peak_amp(1)
corcor_arr = [0.63451,0.62067,0.63453,0.64584,0.51731,0.44686,0.42135,0.26952,0.1324,0.10188]; % the correlation from the plot above
first_nyq_ind = 9; % index (in tau0_arr) where the first Nyquist appears

clf;
scatter(peak_amp_arr(1:first_nyq_ind)'    ,corcor_arr(1:first_nyq_ind)'    ,200,tau_m_arr(1:first_nyq_ind)'    ,'o','filled','MarkerEdgeColor','k'); hold on; % before first Nyquist
scatter(peak_amp_arr(first_nyq_ind+1:end)',corcor_arr(first_nyq_ind+1:end)',200,tau_m_arr(first_nyq_ind+1:end)','^','filled','MarkerEdgeColor','k'); % after first Nyquist
box on; grid on;
colormap(sky());  h = colorbar; 
h.Ticks = [3.4 4.2 5 5.8 6.6];
h.TickLabels = {'t3','t4','t5','t6','t7'};
h.Label.String = '\tau_m';
xlabel('Maximum Amplitude');
ylabel('Optimal vs Final Proj. Correlation (Fig. 5) [r]');
labelpoints(peak_amp_arr+0.1,corcor_arr+0.01,[{'\tau_0=3'} {'\tau_0=4'} {'\tau_0=5'} {'\tau_0=6'} {'\tau_0=7'} {'\tau_0=8'} ...
                   {'\tau_0=9'} {'\tau_0=10'} {'\tau_0=11'} {'\tau_0=12'}],'N',0.15,'FontSize',12)
ylim([0 1]); xlim([0 15]);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 24 16]);
set(gca,'fontsize',12); if sav_fig == 1 & deseas == 2; print('-dpng','-r300',[fig_path,'figS3.png']); end
%%
%[text] ## IX. Error Assessment
%[text] LIM skill is calculated using error variance, but we will also plot the unstandardized error in kilometers ("raw"). Several other metrics will be shown, including calculating the error assuming persistence, a first-order autoregressive forecast, and the LIM's theoretical (expected) error, ε. Skill will also be shown for top-tercile optimal projections, as well as both positive and negative anomalies in ALBSA. The ALBSA anomalies are calculated usiing 1±σ thresholds. Note that the mean of ALBSA in winter is only marginally different than 0 and its distribution is approximately symmetrical with a skewness of 0.33.

% Choose one:
%selected_inds = validation_inds; suffix = 'val';
selected_inds = train_inds; suffix = 'train';

if strcmp(suffix,'val')
    epsilon_tmp = epsilon_val;
elseif strcmp(suffix,'train')
    epsilon_tmp = epsilon_train;
end    

% ALBSA subet
sigth = 1; % standard deviation threshold for anomaly 
indalbsap = find(albsalim(selected_inds)>(mean(albsalim,'omitnan')+sigth*std(albsalim,'omitnan')));
indalbsan = find(albsalim(selected_inds)<(mean(albsalim,'omitnan')-sigth*std(albsalim,'omitnan')));

% top tercile, optProj 
indgoodfx_optprj = find(abs(optProj(selected_inds)) > prctile(abs(optProj),66.66));

xhat                     = squeeze(fx_all(tau_0,:,:,selected_inds))        ; % LIM forecast, xhat, using chosen tau0

% standardize errors
normErrLIM               = NaN(nDat,length(Tau_arr))           ; % standardized LIM error
normErrPers              = NaN(nDat,length(Tau_arr))           ; % standardized persistence error
ar1                      = NaN(nDat,length(Tau_arr))           ; % standardized 1st order autoregressive error
normErrLIMgoodfx_optprj  = NaN(nDat,length(Tau_arr))           ; % standardized LIM error based on optimal projection
normErrLIMgoodfx_rhoinf  = NaN(nDat,length(Tau_arr))           ; % standardized LIM error based on rho infinite
normErrLIMalbsapfx       = NaN(nDat,length(Tau_arr))           ; % standardized LIM error based on positive ALBSA
normErrLIMalbsanfx       = NaN(nDat,length(Tau_arr))           ; % standardized LIM error based on negative ALBSA

% as before, unstandardized
normErrLIMr              = NaN(nDat,length(Tau_arr))           ; % raw LIM error
normErrPersr             = NaN(nDat,length(Tau_arr))           ; % raw persistence error
ar1r                     = NaN(nDat,length(Tau_arr))           ; % raw 1st order autoregressive error
normErrLIMgoodfxr_optprj = NaN(nDat,length(Tau_arr))           ; % raw LIM error based on optimal projection
normErrLIMgoodfxr_rhoinf = NaN(nDat,length(Tau_arr))           ; % raw LIM error based on rho infinite
normErrLIMalbsapfxr      = NaN(nDat,length(Tau_arr))           ; % raw LIM error based on positive ALBSA
normErrLIMalbsanfxr      = NaN(nDat,length(Tau_arr))           ; % raw LIM error based on negative ALBSA

for k = 1:nDat

    for h = 1:length(Tau_arr)

        % Some values will be used more than once so calculate up front
        xobs = x(selected_inds,k);
        xvar = var(x(selected_inds,k),'omitnan');

        xobs_optprj = x(selected_inds(indgoodfx_optprj),k);
        xvar_optprj = var(x(selected_inds(indgoodfx_optprj),k),'omitnan');

        xobs_ap = x(selected_inds(indalbsap),k);
        xvar_ap = var(x(selected_inds(indalbsap),k),'omitnan');

        xobs_an = x(selected_inds(indalbsan),k);
        xvar_an = var(x(selected_inds(indalbsan),k),'omitnan');

        % Calculate standardized errors         
        normErrLIM(k,h) = var(squeeze(xhat(h,k,:))-xobs,'omitnan') / xvar; 
        normErrPers(k,h) = var(xobs(Tau_arr(h)+1:end)-xobs(1:end-Tau_arr(h)),'omitnan') / xvar; 
        aproc = exp((-1/sta_eF(stind(k)))*Tau_arr(h)) * xobs;
        ar1(k,h) = var(xobs(Tau_arr(h)+1:end) - aproc(1:end-Tau_arr(h)),'omitnan') / xvar; 
        normErrLIMgoodfx_optprj(k,h) = var(squeeze(xhat(h,k,indgoodfx_optprj))-xobs_optprj,'omitnan') / xvar_optprj; 
        normErrLIMalbsapfx(k,h) = var(squeeze(xhat(h,k,indalbsap))-xobs_ap,'omitnan') / xvar_ap; 
        normErrLIMalbsanfx(k,h) = var(squeeze(xhat(h,k,indalbsan))-xobs_an,'omitnan') / xvar_an; 
        
        % Calculate unstandardized (raw) errors
        normErrLIMr(k,h) = var(squeeze(xhat(h,k,:))-xobs,'omitnan'); 
        normErrPersr(k,h) = var(xobs(Tau_arr(h)+1:end)-xobs(1:end-Tau_arr(h)),'omitnan'); 
        aprocr = exp((-1/sta_eF(stind(k)))*Tau_arr(h)) * xobs;
        ar1r(k,h) = var(xobs(Tau_arr(h)+1:end) - aproc(1:end-Tau_arr(h)),'omitnan'); 
        normErrLIMgoodfxr_optprj(k,h) = var(squeeze(xhat(h,k,indgoodfx_optprj))-xobs_optprj,'omitnan'); 
        normErrLIMalbsapfxr(k,h) = var(squeeze(xhat(h,k,indalbsap))-xobs_ap,'omitnan'); 
        normErrLIMalbsanfxr(k,h) = var(squeeze(xhat(h,k,indalbsan))-xobs_an,'omitnan'); 

    end

end

ta1 = 0:0.01:20;
nl2 = interp1(Tau_arr,mean(normErrLIM),ta1);
ta1(find(nl2 > 0.4,1,'first'))

cmap = [0.000 0.447 0.741; ... % Blue
        0.850 0.325 0.098; ... % Orange
        0.929 0.694 0.125; ... % Yellow
        0.494 0.184 0.556; ... % Purple
        0.466 0.674 0.188; ... % Green
        0.301 0.745 0.933; ... % Light Blue
        0.635 0.078 0.184];   % Dark Red
clf;
plot([0 Tau_arr],[0 mean(normErrLIM)],'ko-','color',cmap(1,:),'linewidth',2); hold on;
plot([0 Tau_arr],[0 epsilon_tmp(tau_0,:)],'ko-','color',cmap(2,:),'linewidth',1);
plot([0 Tau_arr],[0 mean(ar1)],'ko-','color',cmap(3,:),'linewidth',2);
plot([0 Tau_arr],[0 mean(normErrPers)],'ko-','color',cmap(4,:),'linewidth',2);
plot([0 Tau_arr],[0 mean(normErrLIMgoodfx_optprj)],'ko-','color',cmap(5,:),'linewidth',1);
plot([0 Tau_arr],[0 mean(normErrLIMalbsapfx)],'ko-','color',cmap(6,:),'linewidth',2);
plot([0 Tau_arr],[0 mean(normErrLIMalbsanfx)],'ko-','color',cmap(7,:),'linewidth',2);
grid on; hline(0.4); ylim([0 1.1]); xlim([0 20]);
set(gca,'ytick',0:0.1:1.1,'xtick',0:21);
xlabel('\tau [days]'); ylabel('Norm. Error variance'); 
grid on;
legend('LIM','Expected Error','AR1','Persistence','Optimal Fx','+ALBSA','-ALBSA','location','northwest');
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 30 18]);
set(gca,'fontsize',18); if sav_fig == 1; print('-dpng','-r300',[fig_path,'fig4_relative_',suffix,'_seas',num2str(deseas),'.png']); end

% East vs West
% cmap = lines(256);
% clf;
% plot([0 Tau_arr],[0 mean(normErrLIM)],'ko-','color',cmap(1,:),'linewidth',2); hold on;
% 
% plot([0 Tau_arr],[0 mean(normErrLIM(1:9,:))],'ro-','color',cmap(2,:),'linewidth',2); hold on;
% 
% plot([0 Tau_arr],[0 mean(normErrLIM(1:4,:))],'ro--','color',cmap(2,:),'linewidth',2); hold on;
% 
% plot([0 Tau_arr],[0 mean(normErrLIM(10:18,:))],'go-','color',cmap(3,:),'linewidth',2); hold on;
% 
% plot([0 Tau_arr],[0 mean(normErrLIM(15:18,:))],'go--','color',cmap(3,:),'linewidth',2); hold on;
% 
% grid on; hline(0.4); ylim([0 1.1]); xlim([0 20]);
% set(gca,'ytick',0:0.1:1.1,'xtick',0:21);
% xlabel('\tau [days]'); ylabel('Error variance'); 
% grid on;
% legend('LIM','West','Far West','East','Far East','location','northwest');
% set(gcf,'PaperPositionMode','manual');
% set(gcf,'PaperUnits','centimeters');
% set(gcf,'PaperPosition',[0 0 30 18]);




clf;
plot([0 Tau_arr],[0 mean(normErrLIMr)],'ko-','color',cmap(1,:),'linewidth',2); hold on;
plot([0 Tau_arr],[0 epsilon_raw(tau_0,:)],'ko-','color',cmap(2,:),'linewidth',1);
plot([0 Tau_arr],[0 mean(ar1r)],'ko-','color',cmap(3,:),'linewidth',2);
plot([0 Tau_arr],[0 mean(normErrPersr)],'ko-','color',cmap(4,:),'linewidth',2);
plot([0 Tau_arr],[0 mean(normErrLIMgoodfxr_optprj)],'ko-','color',cmap(5,:),'linewidth',1);
plot([0 Tau_arr],[0 mean(normErrLIMalbsapfxr)],'ko-','color',cmap(6,:),'linewidth',2);
plot([0 Tau_arr],[0 mean(normErrLIMalbsanfxr)],'ko-','color',cmap(7,:),'linewidth',2);
tmp = 0:0.1:1.1;
ylabs = [];
for k = 1:length(tmp)
    ylabs = [ylabs; {[num2str(tmp(k)),' (',num2str(tmp(k)*111),')']}];
end
grid on; hline(0.4); ylim([0 1.1]); xlim([0 20]); 
set(gca,'ytick',0:0.1:1.1,'xtick',0:21,'yticklabel',ylabs);
xlabel('\tau [days]'); ylabel('Mean Error [\circ, km]'); 
%title('Fx Errors'); 
grid on;
legend('LIM','Expected Error','AR1','Persistence','Optimal Fx','+ALBSA','-ALBSA','location','northwest');
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 30 18]);
set(gca,'fontsize',18); if sav_fig == 1; print('-dpng','-r300',[fig_path,'fig4_absolute_',suffix,'_seas',num2str(deseas),'.png']); end
fx_skill_all = [];
fx_skill_albsan = [];
fx_skill_albsap = [];
for k = 1:18

    fx_skill_all = [fx_skill_all; find(normErrLIM(k,:) > 0.4,1,'first')];
    fx_skill_albsan = [fx_skill_albsan; find(normErrLIMalbsanfx(k,:) > 0.4,1,'first')];
    fx_skill_albsap = [fx_skill_albsap; find(normErrLIMalbsapfx(k,:) > 0.4,1,'first')];

end

if strcmp(suffix,'train')
    clf
    bar(181:198,[fx_skill_albsan fx_skill_albsap fx_skill_all])
    legend('Negative ALBSA','Positive ALBSA','All','location','northwest');
    xlabel('Longitude [\circE]');
    ylabel('Skillfuul Forecast Lead Time [days]');
    xlim([180.5 198.5]);
    hline(mean(fx_skill_all));
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperPosition',[0 0 24 16]);
    set(gca,'fontsize',14); if sav_fig == 1 & strcmp(suffix,'train'); print('-dpng','-r300',[fig_path,'fig7_spatialskill.png']); end
end

xx = 0:1/24:max(Tau_arr);
YY = NaN(nDat,5);
for k = 1:nDat
    yy = interp1([0 Tau_arr],[0 normErrPers(k,:)],xx);
    tmp = xx(find(yy > 0.4,1,'first')); if isempty(tmp); tmp = max(Tau_arr); end
    YY(k,1) = tmp;
    yy = interp1([0 Tau_arr],[0 ar1(k,:)],xx);
    tmp = xx(find(yy > 0.4,1,'first')); if isempty(tmp); tmp = max(Tau_arr); end
    YY(k,2) = tmp;
    yy = interp1([0 Tau_arr],[0 normErrLIM(k,:)],xx);
    tmp = xx(find(yy > 0.4,1,'first')); if isempty(tmp); tmp = max(Tau_arr); end
    YY(k,3) = tmp;  
    yy = interp1([0 Tau_arr],[0 normErrLIMgoodfx_optprj(k,:)],xx);
    tmp = xx(find(yy > 0.4,1,'first')); if isempty(tmp); tmp = max(Tau_arr); end
    YY(k,4) = tmp;  
    yy = interp1([0 Tau_arr],[0 normErrLIMalbsanfx(k,:)],xx);
    tmp = xx(find(yy > 0.4,1,'first')); if isempty(tmp); tmp = max(Tau_arr); end
    YY(k,5) = tmp;  
end
eval(['YYrem',num2str(deseas),' = YY;']); % eval...i know. my apologies.
%if strcmp(suffix,'train'); save([data_path,'rem',num2str(deseas),'.mat'], ['YYrem',num2str(deseas)]); end

yy = interp1([0 Tau_arr],[0 mean(normErrPers)],xx);
disp(['Persistence skill at 0.4 ev = ',num2str(xx(find(yy > 0.4,1,'first'))),' days']); 
yy = interp1([0 Tau_arr],[0 mean(ar1)],xx);
disp(['AR1 skill at 0.4 ev = ',num2str(xx(find(yy > 0.4,1,'first'))),' days']); 
yy = interp1([0 Tau_arr],[0 mean(normErrLIM)],xx);
disp(['LIM skill at 0.4 ev = ',num2str(xx(find(yy > 0.4,1,'first'))),' days']); 
yy = interp1([0 Tau_arr],[0 mean(normErrLIMgoodfx_optprj)],xx);
disp(['LIM optimal skill (projection) at 0.4 ev = ',num2str(xx(find(yy > 0.4,1,'first'))),' days']);
yy = interp1([0 Tau_arr],[0 mean(normErrLIMgoodfx_rhoinf)],xx);
disp(['LIM optimal skill (\rho_i_n_f) at 0.4 ev = ',num2str(xx(find(yy > 0.4,1,'first'))),' days']);
%%
%[text] ## X. Case Study
%[text] Plot an example of the time series from the validation set in winter 2023. The top panel is the western side of the domain and the lower panel is the eastern side.
ind = find(yyyy == 2023);
indgoodfx_optprj = find(abs(optProj) > prctile(abs(optProj),66.66));

clf;
subplot(2,1,1)
lonind = find(stas(stind) <= 189);
thedata = mean(x(ind,lonind),2,'omitnan');
fx_all_case = squeeze(fx_all(tau_0,:,lonind,ind));
% rearrange. currently the fx is matched for errors, not for
% forward-looking fx. e.g. lead array at Jan 10 is 1, 2, 3 day fx of Jan 10
% initializing Jan 9, 8, 7 rather than being 1, 2, 3, day fx for Jan 10, 11, 12
fx_all_case_r = NaN.*fx_all_case;
for h = 1:365
    for t = 1:10
        if h+t > 365; continue; end
        fx_all_case_r(t,:,h) = fx_all_case(t,:,h+t);
    end
end

delta_west = [];

h1=plot(iceedge_dn(ind),thedata,'linewidth',2,'color',cmap(2,:)); hold on;

for initday = 1:3:length(ind)
    ds = [];
    val = [];
    for fxlead = 0:7      
        if fxlead == 0
            subi = ind(initday);
            ds = [ds; iceedge_dn(subi)];
            val = [val; thedata(initday)];
            if ismember(ind(initday),indgoodfx_optprj)
                h2b=plot(ds,val,'k.','markersize',12);
                h2a=plot(ds,val,'ko','markersize',12);
            else
                h2b=plot(ds,val,'k.','markersize',12);
            end
        else
            subi = ind(initday)+fxlead-1;
            if subi > length(iceedge_dn); continue; end
            ds = [ds; iceedge_dn(subi)];
            val = [val; mean(fx_all_case_r(fxlead,:,initday),'omitnan')];
        end
    end
    delta_west = [delta_west; val(end)-val(1)];
    h3=plot(ds,val,'k-');
end
set(gca,'xtick',iceedge_dn(ind(1)-1):7:iceedge_dn(ind(1)+92));
xlim([iceedge_dn(ind(1)-1) iceedge_dn(ind(1)+92)]);
datetick('x','mmmdd','keepticks','keeplimits')
grid on;
ylabel('Rel. Latitude [\circ]');
ylim([-3.1 3.1]);
set(gca,'fontsize',14);

yyaxis right

h4=plot(iceedge_dn(ind),albsalim(ind),'linewidth',1,'color',cmap(1,:)); 
ylim([-800 800]);
ylabel('ALBSA [m]');

ax = gca;
ax.YAxis(1).Color = cmap(2,:);
ax.YAxis(2).Color = cmap(1,:);

legend([h1 h2b h2a h3 h4],'Observation','Initialization','Initialization Optimal', ...
    'Forecast','ALBSA','location','northwest');

subplot(2,1,2)
lonind = find(stas(stind) >= 190);
thedata = mean(x(ind,lonind),2,'omitnan');
fx_all_case = squeeze(fx_all(tau_0,:,lonind,ind));
% rearrange. currently the fx is matched for errors, not for
% forward-looking fx. e.g. lead array at Jan 10 is 1, 2, 3 day fx of Jan 10
% initializing Jan 9, 8, 7 rather than being 1, 2, 3, day fx for Jan 10, 11, 12
fx_all_case_r = NaN.*fx_all_case;
for h = 1:365
    for t = 1:10
        if h+t > 365; continue; end
        fx_all_case_r(t,:,h) = fx_all_case(t,:,h+t);
    end
end

delta_east = [];

h1=plot(iceedge_dn(ind),thedata,'linewidth',2,'color',cmap(2,:)); hold on;

for initday = 1:3:length(ind)
    ds = [];
    val = [];
    for fxlead = 0:7      
        if fxlead == 0
            subi = ind(initday);
            ds = [ds; iceedge_dn(subi)];
            val = [val; thedata(initday)];
            if ismember(ind(initday),indgoodfx_optprj)
                h2b=plot(ds,val,'k.','markersize',12);
                h2a=plot(ds,val,'ko','markersize',12);
            else
                h2b=plot(ds,val,'k.','markersize',12);
            end
        else
            subi = ind(initday)+fxlead-1;
            if subi > length(iceedge_dn); continue; end
            ds = [ds; iceedge_dn(subi)];
            val = [val; mean(fx_all_case_r(fxlead,:,initday),'omitnan')];
        end
        delta_east = [delta_east; val(end)-val(1)];
    end
    h3=plot(ds,val,'k-');
end
set(gca,'xtick',iceedge_dn(ind(1)-1):7:iceedge_dn(ind(1)+92));
xlim([iceedge_dn(ind(1)-1) iceedge_dn(ind(1)+92)]);
datetick('x','mmmdd','keepticks','keeplimits')
grid on;
ylabel('Rel. Latitude [\circ]');
ylim([-3.1 3.1]);
set(gca,'fontsize',14);

yyaxis right

h4=plot(iceedge_dn(ind),albsalim(ind),'linewidth',1,'color',cmap(1,:)); 
ylim([-800 800]);
ylabel('ALBSA [m]');

ax = gca;
ax.YAxis(1).Color = cmap(2,:);
ax.YAxis(2).Color = cmap(1,:);

legend([h1 h2b h2a h3 h4],'Observation','Initialization','Initialization Optimal', ...
    'Forecast','ALBSA','location','northwest');

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 28 26]);
set(gca,'fontsize',14); if sav_fig == 1 & deseas == 2; print('-dpng','-r300',[fig_path,'fig3raw']); end
%%
%[text] ## XI. Figure 8, ALBSA maps
%[text] Plot the GPH for +ALBSA and -ALBSA scenarios that were tested earlier. Climatological mean geopotential height at 850 hPa during days from January-March when ALBSA \< -1σ (a) and when ALBSA \> +1σ, corresponding to the cyan and maroon lines in Figure 4.  
indalbsap = find(albsalim>(mean(albsalim,'omitnan')+sigth*std(albsalim,'omitnan')));
indalbsan = find(albsalim<(mean(albsalim,'omitnan')-sigth*std(albsalim,'omitnan')));

indalbsap_val = find(albsalim(validation_inds)>(mean(albsalim,'omitnan')+sigth*std(albsalim,'omitnan')));
indalbsan_val = find(albsalim(validation_inds)<(mean(albsalim,'omitnan')-sigth*std(albsalim,'omitnan')));


u_n = mean(era_uv.u(:,:,validation_inds(indalbsan_val)),3,'omitnan');
v_n = mean(era_uv.v(:,:,validation_inds(indalbsan_val)),3,'omitnan');
u_p = mean(era_uv.u(:,:,validation_inds(indalbsap_val)),3,'omitnan');
v_p = mean(era_uv.v(:,:,validation_inds(indalbsap_val)),3,'omitnan');

clf;
subplot(1,2,1)
cmap = flipud(cbrewer('div', 'BrBG',nDato, 'linear'));
m_proj('lambert','long',[150 240],'lat',[40 85]);
    [~,~] = m_contourf(repmat(albsa.longitude,1,length(albsa.latitude)), ...
        repmat(albsa.latitude,1,length(albsa.longitude))',squeeze(mean(amap(:,:,validation_inds(indalbsan_val)),3,'omitnan')),1150:10:1550,'edgecolor','none'); hold on;
m_quiver(repmat(albsa.longitude,1,length(albsa.latitude)),repmat(albsa.latitude,1,length(albsa.longitude))',u_n,v_n,1,'color','k'); 
m_line([175 175 205 205 175],[50 72 72 50 50],'linewi',2,'color','k');
m_coast('linewidth',1,'color','k');
m_grid('linestyle','none','tickdir','in','linewidth',1);  
clim([1150 1550]);
colormap(cmap)
h = colorbar; h.Label.String = '[m]';
set(gca,'fontsize',14); 
subplot(1,2,2)
cmap = flipud(cbrewer('div', 'BrBG',nDato, 'linear'));
m_proj('lambert','long',[150 240],'lat',[40 85]);
    [~,~] = m_contourf(repmat(albsa.longitude,1,length(albsa.latitude)), ...
        repmat(albsa.latitude,1,length(albsa.longitude))',squeeze(mean(amap(:,:,validation_inds(indalbsap_val)),3,'omitnan')),1150:10:1550,'edgecolor','none'); hold on; 
m_quiver(repmat(albsa.longitude,1,length(albsa.latitude)),repmat(albsa.latitude,1,length(albsa.longitude))',u_p,v_p,1,'color','k'); 
m_line([175 175 205 205 175],[50 72 72 50 50],'linewi',2,'color','k');
m_coast('linewidth',1,'color','k');
m_grid('linestyle','none','tickdir','in','linewidth',1);  
clim([1150 1550]);
colormap(cmap)
h = colorbar; h.Label.String = '[m]';
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 34 20]);
set(gca,'fontsize',14); if sav_fig == 1 & deseas == 2; print('-dpng','-r300',[fig_path,'fig9_raw.png']); end
%%
%[text] ## XII. Figure 9, boxplots
%[text] Boxplots of skillful lead time for four different models (x-axis) and three differences treatments of the seasonal cycle (colors) where blues are time series where only the first component of the seasonal cycle was removed, reds have both the first and second components removed, and grays have the first, second, and third components removed. “Skillful” is defined as the maximum lead time when the error variance \< 0.4 (see Figure 4).
% These data were saved in Section IX
load([data_path,'rem1.mat'])
load([data_path,'rem2.mat'])
load([data_path,'rem3.mat'])

legs = [{'Persistence'} {'AR1'} {'LIM'} {'ALBSA < -1\sigma'} {'LIM (optProj)'}];
xlabs = [{'1st,2nd,3rd'},{'1st,2nd'},{'1st'}];

cmap1 = (cbrewer('seq', 'Blues',nDat+1, 'linear'));
cmap2 = (cbrewer('seq', 'Greys',nDat+1, 'linear'));
cmap3 = (cbrewer('seq', 'Oranges',nDat+1, 'linear'));

data2 = [YYrem3(:,[1 2 3 5 4]); YYrem2(:,[1 2 3 5 4]); YYrem1(:,[1 2 3 5 4])];
group_inx = [ones(1,nDat), 2.*ones(1,nDat), 3.*ones(1,nDat)];

clf;
h = daboxplot(data2,'groups',group_inx,'xtlabels', legs,...
    'fill',0,'whiskers',0,'scatter',2,'outsymbol','k*',...
    'outliers',1,'scattersize',10^-6,'flipcolors',1,'boxspacing',1.2,'mean',1); 
hold on;
for k = 1:nDat
    h11=plot(1.23,YYrem1(k,1),'ko','markerfacecolor',cmap1(k+1,:),'markersize',10,'markeredgecolor','k');
    plot(2.23,YYrem1(k,2),'ko','markerfacecolor',cmap1(k+1,:),'markersize',10,'markeredgecolor','k');
    plot(3.23,YYrem1(k,3),'ko','markerfacecolor',cmap1(k+1,:),'markersize',10,'markeredgecolor','k');
    plot(5.23,YYrem1(k,4),'ko','markerfacecolor',cmap1(k+1,:),'markersize',10,'markeredgecolor','k');
    plot(4.23,YYrem1(k,5),'ko','markerfacecolor',cmap1(k+1,:),'markersize',10,'markeredgecolor','k');

    h33=plot(0.77,YYrem3(k,1),'ko','markerfacecolor',cmap2(k+1,:),'markersize',10,'markeredgecolor','k');
    plot(1.77,YYrem3(k,2),'ko','markerfacecolor',cmap2(k+1,:),'markersize',10,'markeredgecolor','k');
    plot(2.77,YYrem3(k,3),'ko','markerfacecolor',cmap2(k+1,:),'markersize',10,'markeredgecolor','k');
    plot(4.77,YYrem3(k,4),'ko','markerfacecolor',cmap2(k+1,:),'markersize',10,'markeredgecolor','k');
    plot(3.77,YYrem3(k,5),'ko','markerfacecolor',cmap2(k+1,:),'markersize',10,'markeredgecolor','k');

    h22=plot(1,YYrem2(k,1),'ko','markerfacecolor',cmap3(k+1,:),'markersize',10,'markeredgecolor','k');
    plot(2,YYrem2(k,2),'ko','markerfacecolor',cmap3(k+1,:),'markersize',10,'markeredgecolor','k');
    plot(3,YYrem2(k,3),'ko','markerfacecolor',cmap3(k+1,:),'markersize',10,'markeredgecolor','k');
    plot(5,YYrem2(k,4),'ko','markerfacecolor',cmap3(k+1,:),'markersize',10,'markeredgecolor','k');
    plot(4,YYrem2(k,5),'ko','markerfacecolor',cmap3(k+1,:),'markersize',10,'markeredgecolor','k');

    if k == 12
        h1 = h11;
        h2 = h22;
        h3 = h33;
    end
end
plot(0.5,1.1,'ks','markerfacecolor',cmap3(1,:),'markersize',10,'markeredgecolor','k');
plot(0.75,1.1,'ks','markerfacecolor',cmap3(9,:),'markersize',10,'markeredgecolor','k');
text(0.58,1.1,'\rightarrow');
plot(1,1.1,'ks','markerfacecolor',cmap3(18,:),'markersize',10,'markeredgecolor','k');
text(0.82,1.1,'\rightarrow');
text(1.05,1.1,' = 181\circ W \rightarrow 189\circ W \rightarrow 198\circ W');
hh=legend([h1 h2 h3],[{'Rem. 1st'},{'Rem. 1st & 2nd'},{'Rem. 1st, 2nd, & 3rd'}],'location','northwest');
ylabel('Skillful Lead Time [days]'); set(gca,'yscale','log')
grid on; box on; ylim([1 90]);
set(gca,'YTick',[1:10 20:10:100]', 'YTickLabel', cellstr(num2str([1:10 20:10:100]'))); 
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 24 16]);
set(gca,'fontsize',12); if sav_fig == 1 & deseas == 2; print('-dpng','-r300',[fig_path,'fig10.png']); end
%%
%[text] ## XIII. Uncertainty
%[text] 
% I randomly shuffled training/validation, 32/6 years, 50 times and recorded the tau where the mean LIM error variance crosses 0.4 for both sets. The code for this is commened at the top of VI. Do LIM
train_skills = [6.19, 5.55, 5.92, 6.21, 5.81, 6.28, 5.61, 5.82, 5.52, 6.13, 6.02, 6.18, 6.18, 5.71, 5.78, 6.13, 5.98, 5.97, 6.02, 5.94, 5.97, 6.39, 5.77, 5.65, 5.71, ...
                5.92, 6.02, 5.75, 6.32, 5.71, 5.86, 5.65, 5.66, 5.56, 5.69, 6.22, 5.99, 5.35, 5.80, 5.56, 5.98, 5.57, 6.12, 6.24, 5.58, 5.82, 5.78, 5.67, 6.05, 5.74];

val_skills   = [4.75, 6.93, 4.67, 4.08, 5.61, 3.77, 6.97, 5.61, 7.02, 3.93, 4.56, 4.41, 3.96, 6.67, 5.69, 3.89, 4.75, 4.29, 3.84, 4.81, 4.67, 3.62, 5.88, 6.45, 6.32, ...
                4.74, 4.63, 5.85, 3.78, 6.42, 5.48, 6.34, 7.16, 7.91, 6.18, 4.13, 4.63, 8.29, 5.68, 7.91, 4.59, 7.15, 4.21, 4.29, 6.97, 5.64, 5.98, 6.63, 4.74, 5.71];
                
clf;
histogram(train_skills,3:0.1:9); hold on; 
histogram(val_skills,3:0.1:9);

mean(train_skills)
std(train_skills)

mean(val_skills)
std(val_skills)
histogram(val_skills-train_skills,-3:0.1:3)

% A more robust approach to respect the chronology of the time series is a fixed-origin rolling cross-validation, here implemented in 19 steps beginning 1988-2001 (train) / 2002-2007 (val) and going to 1988-2019 (train) / 2020-2025 (val)
% results are very similar
train_skills = [5.52 5.53 5.65 5.83 5.84 5.98 6.01 6.00 6.00 5.98 6.01 6.06 5.66 5.16 5.05 5.12 5.14 4.96 5.05];
val_skills = [6.97 6.56 5.64 4.73 4.73 4.28 3.77 3.62 3.89 4.66 4.38 5.02 6.63 7.95 7.94 7.85 7.81 8.27 6.47];
%[text]  
%[text] 
%[text] 
%[text] 
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":29.8}
%---
