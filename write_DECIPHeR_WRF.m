% Write DECIPHeR input
%  This script writes DECIPHeR input files into the NEW HOURLY format (separate
%  files for discharge, rainfall and PET instead of one confusing file
%  containing everything). 
%
%  Rosie Lane - 10/11/2017
%  Added to by Gemma Coxon 29th Jan 2018 to work with the new HRU Meta
%  files instead of the functional units file
%  Iskra 01/11/2018 Edited for the format for hourly discharge, PET and
%  1 km gridded rainfall.

clear;

cd('/home/pm15021/Matlab_hourly')

start_date = '27/06/2012';
end_date = '29/06/2012';
stunde=12; %time in which WRF simulation started. If it's whole day, write 0
filename_output = 'Tyne_23016_';
header_output = 'Extracted for Ouseburn at Crag Hall 23016';
folder_output = '/home/pm15021/DECIPHeR_v2/GEAR_valid_Newc_27b/INPUTS';
gauge_for_analysis = 23016;
hrumeta_file = '/home/pm15021/DECIPHeR_v2/HRU_23016_WRF/23016000_hru_meta.dat';
catcharea_file = '/home/pm15021/DECIPHeR_v2/DTA_23016/masks/mask_info.txt';
folder_Q = '/home/pm15021/Hourly_Discharge';
folder_rain = '/home/pm15021/Hourly_Rainfall/WRF_rainfall_Newcastle/NewcastleWRF27b';    
folder_PET = '/home/pm15021/Hourly_PET';


%defaults:
nodata_value = -9999;

% NB. A single PET is used as input and based on the largest catchment
% area (i.e. the downstream gauge)
% Currently based on 10km rainfall

dates = (datenum(start_date, 'dd/mm/yyyy'):1:datenum(end_date, 'dd/mm/yyyy'))';
dates24 = repmat(dates,1,24)';
dates24 = dates24(:);
spinup = size(dates24,1);
dates24 = dates24(stunde+1:spinup);
no_timesteps = size(dates24,1);
% dates24 = repmat(dates,[24 1]);
% dates24 = dates24(:);

%Read in HRU file - looking for gauge numbers, rainfall and PET grid
%numbers so first just read in the header line
fid = fopen(hrumeta_file);
for i = 1 : 5
tmp = fgetl(fid);
end
%hrumeta_headerline = strsplit(tmp, ' ');
hrumeta_headerline = regexp(tmp,'  ','split');
fclose(fid);
%Then read in the rest of the data
hrumeta = dlmread(hrumeta_file, '\t', 5, 0);

% Number of gauges is always in column 5
n_cat = max(hrumeta(:, 5)); 

% Start a log file to record discharge data processing
if exist(folder_output,'dir') ~= 7
    fprintf('\nERROR: Output folder given does not exist.')
    fprintf('\ncannot find folder %s',folder_output)
    cd ~
    mkdir('DYMOND_INPUT_FILES')
    folder_output = '~/DYMOND_INPUT_FILES';
end

fid_log = fopen([folder_output,'/',filename_output,'_WriteInput_log.txt'],'w');
fprintf(fid_log,'Log file: write_DYMOND_input_Rosie');
fprintf(fid_log,'\nheader_output = %s',header_output);
fprintf(fid_log,'\nfolder_rain = %s',folder_rain);
fprintf(fid_log,'\nfolder_PET = %s',folder_PET);
fprintf(fid_log,'\nfolder_Q = %s',folder_Q);
fprintf(fid_log,'\nExtracting data from %s to %s',start_date, end_date);
fprintf(fid_log,'\nHRU Meta Data file selected: %s',hrumeta_file);
fprintf(fid_log,'\nCatchment area file selected: %s',catcharea_file);

% Process the discharge data
fprintf(fid_log,'\n\n===================================================');
fprintf(fid_log,'\n\tProcessing the discharge data');
fprintf(fid_log,'\n===================================================');
fprintf('\n\nProcessing the dicharge data:')
%check given filename exists
if exist(folder_Q,'dir') == 7
    cd(folder_Q)
else
    fprintf(fid_log,'\nERROR: Discharge input folder given does not exist.');
    fprintf('\nERROR: Discharge input folder given does not exist.');
    fprintf(fid_log,'\ncannot find folder %s',folder_Q);
    fprintf('\ncannot find folder %s',folder_Q);
    return
end

%1a. Get catchment areas 
if exist(catcharea_file)==2
fid = fopen(catcharea_file);
fgetl(fid);
area_all = textscan(fid,'%f %f %f %f %f %*[^\n]');
area_all = cell2mat(area_all);      %catch_id, x, y, area
fclose(fid);
else
    fprintf('\nERROR: mask_info.txt file missing.')
    fprintf('\nCannot find file: %s',catcharea_file)
    return
end

%1b. create empty discharge array for each catchment
obsq = nan(length(dates24),n_cat); %row for each day, col for each catch_id

%1c. grab discharge data for all catchments!

% Find the column from the HRU meta file that contains all the gauge
% numbers
id_cat = ismember(hrumeta_headerline, 'GAUGE_NUM');

for c = 1:n_cat
    %get index of catchment
    ix = find(hrumeta(:, 5) == c, 1, 'first');
    catchid(c) = hrumeta(ix, id_cat)/1000; %discharge_file
    catcharea(c) = area_all(catchid(c)*1000 == area_all(:, 1), 5) / 1000000;
    
    if catchid(c)>0 % If catchid(c) is an actual gauge number 
        area_avail=true;
        fprintf('\nValid point c= %d, catchid= %d, catcharea= %d',c,catchid(c),catcharea(c));
        
        % If the file exists, then read it.
        if exist([num2str(catchid(c)) '_to_HOURLY.csv'], 'file') == 2 && area_avail
            obsqall = load([num2str(catchid(c)) '_to_HOURLY.csv']);
            %convert obsqall from m3/seconds to mm/day
            obsqall(:, 5)= (obsqall(:, 5)*60*60*24)/(catcharea(c)*1000);

            qdates = datenum(obsqall(:, 3), obsqall(:, 2), obsqall(:, 1));
            ib = find(qdates(:,1) == dates(1,1), 1, 'first') + 1; %first time date appears is at 24 am, but days actually start at 1 am.
            [~, ia, dummy] = intersect(dates, qdates);
            %print warning if discharge data not available for selected dataes
            if isempty(ia) == 1
                fprintf(fid_log,'\ncatchment %d / %d. No discharge data for %d for chosen time period',c,n_cat,catchid(c));
                fprintf('\ncatchment %d / %d. No discharge data for %d for chosen time period',c,n_cat,catchid(c));
            else
                erste = ib+stunde;
                letzte = erste+no_timesteps-1;
                ib24 = (erste:letzte)'; % Row indexes of the matrix with all data, for the selcted period 
                obsq(1:no_timesteps, c) = obsqall(ib24, 5);
            end
        % If the file doesn't exist then it is usually ungauged
        else
            fprintf('\ncatchment %d / %d. %d_to_HOURLY.csv does not exist - assuming ungauged point and filling blanks with NaNs',c,n_cat,catchid(c));
            fprintf(fid_log,'\ncatchment %d / %d. %d_to_HOURLY.csv does not exist - assuming ungauged point and filling blanks with NaNs',c,n_cat,catchid(c));
            obsq(:, c) = NaN;
        end
    else % If catchid(c) is NOT a gauge number
        fprintf('\nIgnoring point c= %d, catchid= %d, catcharea= %d',c,catchid(c),catcharea(c));
        area_avail=false;
    end 
end
header_second_discharge = 'Extracted from the National River Flow Archive, aggregated from 15 minute data.';

%===============================================================================================================

fprintf(fid_log,'\n\n===================================================');
fprintf(fid_log,'\n\tProcessing the precipitation data');
fprintf(fid_log,'\n===================================================');
fprintf('\n\nProcessing the precipitation data:')
if exist(folder_rain,'dir') == 7
    cd(folder_rain)
else
    fprintf('\nERROR: Rainfall input folder given does not exist.')
    fprintf('\ncannot find folder %s',folder_rain)
        fprintf(fid_log,'\nERROR: Rainfall input folder given does not exist.');
    fprintf(fid_log,'\ncannot find folder %s',folder_rain);
    return
end

%grab the data - format depending on source
if strfind(folder_rain,'GEAR_1km_grid') > 0    %rain input is GEAR 10km
    
    % Looking for rainfall grid IDs
    id_rain = ismember(hrumeta_headerline, 'RAIN_NUM');
    id2_rain = ismember(hrumeta_headerline, 'RAIN_ID');
    
    n_rain = max(hrumeta(:, id2_rain));
    
    fprintf('\nRainfall input is GEAR 1km grid data.')
    fprintf(fid_log,'\nRainfall input is GEAR 1km grid data.');
    for i = 1 : n_rain
        ix = find(hrumeta(:, id2_rain) == i, 1, 'first');
        raingridid(i) = hrumeta(ix, id_rain); %10
        rain_gridcell = load(['rainfall_1km_gridID' num2str(raingridid(i)) '.txt']);
        rain_dates = datenum(rain_gridcell(:, 1), rain_gridcell(:, 2), rain_gridcell(:, 3));
        ib = find(rain_dates(:,1) == dates(1,1), 1, 'first') + 1; %first time date appears is at 24 hrs, but days actually start at 1 am.
        [~, ia, dummy] = intersect(dates, qdates);
        ib24 = (ib(1,1):(ib(1,1)+no_timesteps-1))'; % Row indexes from rainfall file with all dates, for selected period.
        daily_rain(:, i) = rain_gridcell(ib24,5);
    end
    header_second_rain = 'Extracted from CEH GEAR dataset, 1km cells.';

elseif strfind(folder_rain,'WRF') > 0    %rain input is from WRF 
    
    % Looking for rainfall grid IDs
    id_rain = ismember(hrumeta_headerline, 'RAIN_NUM');
    id2_rain = ismember(hrumeta_headerline, 'RAIN_ID');
    
    n_rain = max(hrumeta(:, id2_rain));
    
    fprintf('\nRainfall input is WRF 2 km grid data.')
    fprintf(fid_log,'\nRainfall input is WRF 2 km grid data.');
    for i = 1 : n_rain
        ix = find(hrumeta(:, id2_rain) == i, 1, 'first');
        raingridid(i) = hrumeta(ix, id_rain);
        rain_gridcell = load(['rainfall_WRF_gridID' num2str(raingridid(i)) '.txt']);
        rain_dates = datenum(rain_gridcell(:, 1), rain_gridcell(:, 2), rain_gridcell(:, 3));
        %daily_rain(:, i) = rain_gridcell(ismember(rain_dates, dates), 4);
        daily_rain(:, i) = rain_gridcell(:, 5);
    end
    header_second_rain = 'Extracted from WRF 2 km grid.';


elseif strfind(folder_rain,'GEAR_5km_grid') > 0    %rain input is GEAR 10km
    
    % Looking for rainfall grid IDs
    id_rain = ismember(hrumeta_headerline, 'RAIN_NUM');
    id2_rain = ismember(hrumeta_headerline, 'RAIN_ID');
    
    n_rain = max(hrumeta(:, id2_rain));
    
    fprintf('\nRainfall input is GEAR 5km grid data.')
    fprintf(fid_log,'\nRainfall input is GEAR 5km grid data.');
    for i = 1 : n_rain
        ix = find(hrumeta(:, id2_rain) == i, 1, 'first');
        raingridid(i) = hrumeta(ix, id_rain);
        rain_gridcell = load(['CEHGEAR_5km_gridID' num2str(raingridid(i)) '.txt']);
        rain_dates = datenum(rain_gridcell(:, 1), rain_gridcell(:, 2), rain_gridcell(:, 3));
        daily_rain(:, i) = rain_gridcell(ismember(rain_dates, dates), 4);
    end
    header_second_rain = 'Extracted from CEH GEAR dataset, 1km cells averaged to 5km.';
    
elseif strfind(folder_rain,'GEAR_lumped') > 0   %rain input is GEAR lumped data
    fprintf('\nRainfall input is GEAR lumped data.')
    fprintf(fid_log,'\nRainfall input is GEAR lumped data.');
    fname_rain = ['GEAR_daily_rainfall_',num2str(catchid(catcharea==max(catcharea))),'.txt'];
    rainall=load([folder_rain,'/',fname_rain]);
    fprintf('\nLumped data used: %s',fname_rain)
    raindates = datenum(rainall(:, 1), rainall(:, 2), rainall(:, 3));
    daily_rain = rainall(ismember(raindates, dates), 4);
    header_second_rain = 'Extracted from CEH GEAR dataset, using CEH catchment outlines.';
    n_rain=1;
    
elseif strfind(folder_rain, 'GEAR_WAH_grid') > 0 %rain input is GEAH WAH grid
    
    % Looking for rainfall grid IDs
    id_rain = ismember(hrumeta_headerline, 'RAIN_NUM');
    id2_rain = ismember(hrumeta_headerline, 'RAIN_ID');
    
    n_rain = max(hrumeta(:, id2_rain));
    
    fprintf('\nRainfall input is GEAR WAH grid data.')
    fprintf(fid_log,'\nRainfall input is GEAR WAH grid data.');
    for i = 1 : n_rain
        ix = find(hrumeta(:, id2_rain) == i, 1, 'first');
        raingridid(i) = hrumeta(ix, id_rain);
        rain_gridcell = load(['rainfall_WAH_gridID' num2str(raingridid(i)) '.txt']);
        rain_dates = datenum(rain_gridcell(:, 1), rain_gridcell(:, 2), rain_gridcell(:, 3));
        daily_rain(:, i) = rain_gridcell(ismember(rain_dates, dates), 4);
    end 
    header_second_rain = 'Extracted from the CEH GEAR dataset. WAH grid.';
    
else                                            %rain input is not a recognised source
    fprintf('\nERROR: Rainfall input is not a recognised format.')
    fprintf('\nOnly GEAR_10km_grid and GEAR_lumped data can automatically be processed.')
    fprintf(fid_log,'\nERROR: Rainfall input is not a recognised format.');
    fprintf(fid_log,'\nOnly GEAR_10km_grid and GEAR_lumped data can automatically be processed.');
    return
end
    
% Process the PET data
fprintf(fid_log,'\n\n===================================================');
fprintf(fid_log,'\n\tProcessing the PET data');
fprintf(fid_log,'\n===================================================');
fprintf('\n\nProcessing the PET data:')

%Create a matrix of gauges and areas without negative numbers
catchmatrix = [catchid',catcharea'];
catchmatrix(catchmatrix(:,1) < 0,:) = [];

%--------------IMPORTANT EDIT---------------------
%Keep only gauges in catchment (i.e. if analysed gauge is 23016, keep
%gauge numbers 23XXX)
catchmatrix = [catchmatrix,(1:size(catchmatrix,1))']; %This creates a ghost 3rd colum
catchmatrix(:,3) = floor(catchmatrix(:,1)/1000); %3rd column is now first two number of gauges
catchmatrix = catchmatrix(catchmatrix(:,3) == floor(gauge_for_analysis/1000),:);

if exist(folder_PET,'dir') == 7
    cd(folder_PET)
else
    fprintf('\nERROR: PET folder given does not exist.')
    fprintf('\ncannot find folder %s',folder_PET)
   fprintf(fid_log,'\nERROR: PET folder given does not exist.');
    fprintf(fid_log,'\ncannot find folder %s',folder_PET);
    return
end

%grab the data - format depending on source
if strfind(folder_PET,'CHESS_PET_10km_Grid') > 0    %PET input is CHESS PET 10km
    
    % Looking for PET grid IDs
    id_pet = ismember(hrumeta_headerline, 'PET_NUM');
    id2_pet = ismember(hrumeta_headerline, 'PET_ID');
    
    n_pet = max(hrumeta(:, id2_pet));
    
    fprintf('\nPET input is CHESS PET 10km grid data.')
    fprintf(fid_log,'\nPET input is CHESS PET 10km grid data.');
    for i = 1 : n_pet
        ix = find(hrumeta(:, id2_pet) == i, 1, 'first');
        petgridid(i) = hrumeta(ix, id_pet);
        pet_gridcell = load(['CHESSPET_10km_gridID' num2str(petgridid(i)) '.txt']);
        pet_dates = datenum(pet_gridcell(:, 1), pet_gridcell(:, 2), pet_gridcell(:, 3));
        evap(:, i) = pet_gridcell(ismember(pet_dates, dates), 4);
    end
    header_second_pet = 'Extracted from CHESS PET dataset, 1km cells averaged to 10km.';
    
elseif strfind(folder_PET, 'Hourly_PET') > 0 
    fname_PET = ['hourly_peti_1961_2015_' num2str(catchmatrix(catchmatrix(:,2) == max(catchmatrix(:,2)))) '.txt'];
    fprintf('\nLumped hourly PETi data used: %s',fname_PET)
    fprintf(fid_log,'\nLumped hourly PETi data used: %s',fname_PET);
    evapall=dlmread(fname_PET, ',', 1, 0);
    evapdates = datenum(evapall(:, 1), evapall(:, 2), evapall(:, 3));
       ib = find(evapdates(:,1) == dates(1,1), 1, 'first'); %first time date appears is at 1 am.
       erste = ib+stunde;
       letzte = erste+no_timesteps-1;
       ib24 = (erste:letzte)'; % Row indexes of the matrix with all data, for the selcted period 
    evap = evapall(ib24,5);
    header_second_pet = 'Extracted from the Hourly CHESS PETi dataset.';
    n_pet=1;

elseif strfind(folder_PET,'CHESS_PET_5km_Grid') > 0    %PET input is CHESS PET 10km
    
    % Looking for PET grid IDs
    id_pet = ismember(hrumeta_headerline, 'PET_NUM');
    id2_pet = ismember(hrumeta_headerline, 'PET_ID');
    
    n_pet = max(hrumeta(:, id2_pet));
    
    fprintf('\nPET input is CHESS PET 5km grid data.')
    fprintf(fid_log,'\nPET input is CHESS PET 5km grid data.');
    for i = 1 : n_pet
        ix = find(hrumeta(:, id2_pet) == i, 1, 'first');
        petgridid(i) = hrumeta(ix, id_pet);
        pet_gridcell = load(['CHESSPET_5km_gridID' num2str(petgridid(i)) '.txt']);
        pet_dates = datenum(pet_gridcell(:, 1), pet_gridcell(:, 2), pet_gridcell(:, 3));
        evap(:, i) = pet_gridcell(ismember(pet_dates, dates), 4);
    end
    header_second_pet = 'Extracted from CHESS PET dataset, 1km cells averaged to 5km.';
   
elseif strfind(folder_PET, 'CHESS_PET_lumped') > 0 
    fname_PET = ['daily_pet_1961_2015_' num2str(catchid(catcharea == max(catcharea))) '.txt'];
    fprintf('\nLumped data used: %s',fname_PET)
    fprintf(fid_log,'\nLumped data used: %s',fname_PET);
    evapall= load(fname_PET);
    evapdates = datenum(evapall(:, 1), evapall(:, 2), evapall(:, 3));
    evap= evapall(ismember(evapdates, dates), 4);
    header_second_pet = 'Extracted from the CHESS PET dataset.';
    n_pet=1;

elseif strfind(folder_PET, 'CHESS_PET_WAH') >0
    fname_PET = ['daily_chesspet_wah_' num2str(catchid(catcharea == max(catcharea))) '.txt'];
    fprintf('\nLumped data used: %s',fname_PET)
    fprintf(fid_log,'\nLumped data used: %s',fname_PET);
    evapall= load(fname_PET);
    evapdates = datenum(evapall(:, 1), evapall(:, 2), evapall(:, 3));
    evap= evapall(ismember(evapdates, dates), 4);
    header_second_pet = 'Extracted from the CHESS PET (WAH) dataset.';
   
elseif strfind(folder_PET,'CHESS_PETi_WAH') == 1
    fname_PET = ['daily_chesspeti_wah_' num2str(catchid(catcharea == max(catcharea))) '.txt'];
    fprintf('\nLumped data used: %s',fname_PET)
    fprintf(fid_log,'\nLumped data used: %s',fname_PET);
    evapall= load(fname_PET);
    evapdates = datenum(evapall(:, 1), evapall(:, 2), evapall(:, 3));
    evap= evapall(ismember(evapdates, dates), 4);  
    header_second_pet = 'Extracted from the CHESS PETi WAH dataset.';
    
else
    fprintf('\nERROR: evap input is not in a recognised format')
    fprintf('\nRecommended to use CHESS PET')  
        fprintf(fid_log,'\nERROR: evap input is not in a recognised format');
    fprintf(fid_log,'\nRecommended to use CHESS PET')  ;
end

% Convert data
%convert from mm/day to m/day
input_rain = daily_rain./1000;
input_evap = evap./1000;
input_obsq = obsq./1000;

%convert missing data to nodata_value
input_obsq(isnan(input_obsq)) = nodata_value;
input_rain(isnan(input_rain)) = nodata_value;
input_evap(isnan(input_evap)) = nodata_value;

% Get dates and headers
input_dates(:,1) = str2num(datestr(dates24,'yyyy'));
input_dates(:,2) = str2num(datestr(dates24,'mm'));
input_dates(:,3) = str2num(datestr(dates24,'dd'));
hour24 = repmat(1:24,1,size(dates,1))'; % Repeat 1:24 for each day of the daily timeseries.
hour24 = hour24(stunde+1:size(hour24));
input_dates(:,4) = hour24;

all_rain = [input_dates, input_rain];
all_evap = [input_dates, input_evap];
all_obsq = [input_dates, input_obsq];

%create headers
header_first{1} = ['UK hourly '];
header_first{2} = ['data (m/day) ',datestr(dates(1)),' --> ',datestr(dates(end)),'.'];
header_rain{1} = 'YYYY'; header_rain{2} = 'MM'; header_rain{3} = 'DD'; header_rain{4} = 'HH';
header_obsq{1} = 'YYYY'; header_obsq{2} = 'MM'; header_obsq{3} = 'DD'; header_obsq{4} = 'HH';
header_evap{1} = 'YYYY'; header_evap{2} = 'MM'; header_evap{3} = 'DD'; header_evap{4} = 'HH';

for i = 1:n_pet
    header_evap{i+4} = ['PETID_',num2str(i)];
end

for i = 1:n_rain
    header_rain{i+4} = ['RainID_',num2str(i)];
end

for i = 1:n_cat
    header_obsq{i+4} = ['CatID_',num2str(i)];
end

% Save DECIPHeR input files
fprintf(fid_log,'\n\n===================================================');
fprintf(fid_log,'\n\tSaving DECIPHeR input files');
fprintf(fid_log,'\n===================================================');
fprintf('\n\nSaving DECIPHeR input files.')
fprintf('\nSaving files in: %s', folder_output)
fprintf(fid_log,'\nSaving files in: %s', folder_output);
cd(folder_output)

% save pet file
fid = fopen([filename_output '_PET.dat'], 'w');
fprintf(fid,'%sPET %s',header_first{1},header_first{2});
fprintf(fid,'\n%s',header_output);
fprintf(fid,'\n%s',header_second_pet);
fprintf(fid,'\n1\t %d \t%d \t %d \t(timestep_hours, n_timesteps, n_PETIDs, nodata_value)\n\n',length(dates24),n_pet,nodata_value);
for i = 1:length(header_evap)
    if i<=4
        fprintf(fid,'%s       ',header_evap{i});
    else
        fprintf(fid,'%s            ',header_evap{i});
    end
end
formatString = repmat('%4.10f     ', 1, size(all_evap, 2)-4);
formatString = ['\n%4d     %4d     %4d     %4d     ',formatString];
for k = 1 : size(all_evap, 1)
    fprintf(fid, formatString, all_evap(k, :));
end
fclose(fid);

fprintf('\nPET file saved successfully: %s/%s_PET.dat',folder_output,filename_output)
fprintf(fid_log,'\nPET file saved successfully: %s/%s_PET.dat',folder_output,filename_output);

%formatted discharge
fid = fopen([filename_output '_obsq.dat'], 'w');
fprintf(fid,'%sdischarge %s',header_first{1},header_first{2});
fprintf(fid,'\n%s',header_output);
fprintf(fid,'\n%s',header_second_discharge);
fprintf(fid,'\n1\t %d \t%d \t %d \t(timestep_hours, n_timesteps, n_catchments, nodata_value)\n\n',length(dates24),n_cat,nodata_value);
for i = 1:length(header_obsq)
    if i<=4
        fprintf(fid,'%s       ',header_obsq{i});
    else
        fprintf(fid,'%s            ',header_obsq{i});
    end
end
formatString = repmat('%4.10f     ', 1, size(all_obsq, 2)-4);
formatString = ['\n%4d     %4d     %4d     %4d     ',formatString];
for k = 1 : size(all_obsq, 1)
    fprintf(fid, formatString, all_obsq(k, :));
end
fclose(fid);
fprintf('\nDischarge file saved successfully: %s/%s_obsq.dat',folder_output,filename_output)
fprintf(fid_log,'\nDischarge file saved successfully: %s/%s_obsq.dat',folder_output,filename_output);

% save rain file
fid = fopen([filename_output '_rain.dat'], 'w');
fprintf(fid,'%srainfall %s',header_first{1},header_first{2});
fprintf(fid,'\n%s',header_output);
fprintf(fid,'\n%s',header_second_rain);
fprintf(fid,'\n1\t %d \t%d \t %d \t(timestep_hours, n_timesteps, n_RainIDs, nodata_value)\n\n',length(dates24),n_rain,nodata_value);
for i = 1:length(header_rain)
    if i<=4
        fprintf(fid,'%s       ',header_rain{i});
    else
        fprintf(fid,'%s        ',header_rain{i});
    end
end
formatString = repmat('%4.10f     ', 1, size(all_rain, 2)-4);
formatString = ['\n%4d     %4d     %4d     %4d     ',formatString];
for k = 1 : size(all_rain, 1)
    fprintf(fid, formatString, all_rain(k, :));
end
fclose(fid);
%dlmwrite([filename_output '_rain.dat'],all_rain,'-append','delimiter','\t',...
%    'roffset',1)
fprintf('\nRainfall file saved successfully: %s/%s_rain.dat\n',folder_output,filename_output)
fprintf(fid_log,'\nRainfall file saved successfully: %s/%s_rain.dat',folder_output,filename_output);
fclose(fid_log);

%end
