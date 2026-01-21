% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHORS:
%
%   Christopher Cox (NOAA) christopher.j.cox@noaa.gov
%
% REFERENCES:
%
%   Fetterer, F. Knowles, K., Meier, W. N., Savoie, M., and Windnagel, 
%       A. K. (2017), Sea Ice Index (G02135, Version 3). [Dataset]. 
%       National Snow and Ice Data Center. Doi:10.7265/N5K072F8
% 
%
% PURPOSE:
%
%   Loads sea ice index data and does some processing for Cox and Penland
%   (in prep), as described in the paper.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -- OPTIONS --------------------------------------------------------------

% Starting year
yr_start = 1987; % because that is when the switched to daily...

% Stations (meridions)
stas = 175:200; % degrees longitude

% Output to geographic grid with the following specifications
lon = 0:1:360; % degrees longitude
lat = 50:0.1:90; % degrees latitude

% -------------------------------------------------------------------------


% -- MAIN -----------------------------------------------------------------

% Data directory
data_dir = '/path';
d = [dir([data_dir,'1*']);dir([data_dir,'2*'])];
for k = 1:length(d); if strcmp(d(k).name,num2str(yr_start)); break; end; end
d = d(k:end); % begin at yr_start 

% Load the grid, NSIDC's polar sterographic projection
fid = fopen([data_dir,'psn25lats_v3.dat']);
grlats = fread(fid,[304 448],'int',0,'n')./100000;
fclose(fid);
fid = fopen([data_dir,'psn25lons_v3.dat']);
grlons = fread(fid,[304 448],'int',0,'n')./100000;
grlons = mod(grlons,360);
fclose(fid);

% Regrid to a geographic grid, 0.1 degrees
latflip = fliplr(lat);
[xq,yq] = meshgrid(lon,lat);

% an array of (STAS x DAYS) and another for time in matlab serial dates
iceedge = NaN(length(stas),length(d)*366);
iceedge_dn = NaN(1,length(d)*366);

% Load the data. This takes about an hour on my laptop.
count = 1;
for y = 1:length(d)

    cd([d(y).folder,'/',d(y).name]);
    dd = [dir('0*_*');dir('1*_*')]; 

    for m = 1:length(dd)

        disp(['Month: ',num2str(m),' Year: ',d(y).name])
        cd([dd(m).folder,'/',dd(m).name]);

        ddd = dir('*_concentration_*');

        for i = 1:length(ddd)
            file = importdata(ddd(i).name);
            sic = griddata(grlons,grlats,double(file.cdata'),double(xq),double(yq))./10;
            
            iceedge_dn(count) = datenum([str2num(ddd(i).name(3:6)) str2num(ddd(i).name(7:8)) str2num(ddd(i).name(9:10)) 0 0 0]);
            for l = stas

                tmp = flipud(sic(:,l));

                % this buffers the land, which due to the regrid has become fuzzy
                tmp2 = tmp;
                tmp(tmp > 100) = 254;
                ind = find(tmp == 254);
                ind2 = [];
                for k = 1:length(ind)
                    ind2 = [ind2; ind(k)-3 ind(k)-2 ind(k)-1 ind(k) ind(k)+1 ind(k)+2 ind(k)+3];
                end
                ind2 = unique(ind2);
                ind2(ind2 > length(tmp) | ind2 < 1) = [];
                tmp2(ind2) = 254;
                tmp = tmp2;
                    
                indice = find(tmp <= 15,1,'first');
                if indice > 1 && indice < length(tmp)
                    if (tmp(indice+1) == 254) || (tmp(indice-1) == 254)
                        % then ignore, we are up against land. shall be nan
                    else
                        iceedge(l-(stas(1)-1),count) = latflip(indice);
                    end
                end
            end
            count = count+1;
        end

    end
    
end

iceedge(:,count:end) = [];
iceedge_dn(count:end) = [];

clearallbut iceedge iceedge_dn
save nsidc_data

% -------------------------------------------------------------------------