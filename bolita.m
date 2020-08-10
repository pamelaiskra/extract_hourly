%make_raingrid_hourly()

%this produces hourly timeseries for each 10km cell
%data taken from the dissag CEH GEAR layers

tic;
fprintf('Starting\n')

%load GC 10km grid
cd '/projects/The_Env_Virtual_observatory/DynaTOP_data/Input_Grids'

rain_grid = dlmread('UK_input_grid_1km.asc',' ', 6,0);

GC_10km = rain_grid; clear rain_grid;

fid = fopen('UK_input_grid_1km.asc');

fprintf('Checkpoint_1\n')

tmp = textscan(fid, '%s%f', 1);
ncols_RG = tmp{1,2};
tmp = textscan(fid, '%s%f', 1);
nrows_RG = tmp{1,2};
tmp = textscan(fid, '%s%f', 1);
xll_RG = tmp{1,2};
tmp = textscan(fid, '%s%f', 1);
yll_RG = tmp{1,2};
tmp = textscan(fid, '%s%f', 1);
cellsize_RG = tmp{1,2};
tmp = textscan(fid, '%s%f', 1);
exclude_RG = tmp{1,2};

fprintf('Checkpoint_2\n')

%create a refernce for each 10km grid and the associated cells we will take
%from the CEH sub daily grids
Max_ID = max(GC_10km(:));
IDs = 1:1:Max_ID; IDs = IDs';
display(length(IDs(:,1)))
for ID = 1:1:length(IDs(:,1));
    
    display(ID)

    [X,Y] = find (GC_10km == ID); %find cell coordinates we care about
    Min_X = min(Y); Min_Y = min(X);
    Max_X = max(Y); Max_Y = max(X); %corner CELLS
    %get the cells of interest from 1km layer
    rain_LLx = floor((Min_X / 20)+1); %20 of gc 50m cells in my 1km grid
    rain_LRx = floor((Max_X/20)+1); 
   rain_LLy = 1251 - (floor((length(GC_10km(:,1)) - Min_Y)/20)+1);
   rain_ULy = 1251 - (floor((length(GC_10km(:,1)) - Max_Y)/20)+1); 
  
    %assign to linkage table
    IDs(ID,2) = rain_LLx;
    IDs(ID,3) = rain_LRx;
    IDs(ID,4) = rain_LLy;
    IDs(ID,5) = rain_ULy;
end;

fprintf('Checkpoint_3\n')

%create timeseries to fill (1st jan 1990 10:00 - 1st jan 2013 09:00)
DS = '2003/01/01 10:00'; DE = '2008/01/01 09:00';
DS = datenum(DS); DE=datenum(DE);
dates = DS: (1/24) : DE; dates = dates';

% Pre-Allocate huge matrix to hold each ID data.
% size = dates, Max_ID;
VALUES = dates(:,1);
VALUES (:,2827) = NaN; VALUES = VALUES .*NaN;


%manipulate each rain day - add to value file
YEAR = 2003;
Row = 1; %row in the data we want to add to
Folder = ([num2str(YEAR), '_TBR_Disag']);
Folder = (['/projects/The_Env_Virtual_observatory/Sub_Daily_Hydrological_Data/Rainfall_Data/SubDaily_CEH/Grid_data/', Folder]);
cd (Folder)

fprintf('Checkpoint_4\n')

display(length(dates))
difff = 81168-76021+1;

tic;
for i = 1:24:length(dates); %10am on each day
    
    display(i)
   
    YEARi = datestr(dates(i), 'yyyy');
    
    if str2num(YEARi) > YEAR;
    YEAR = YEAR + 1;
    Folder = ([num2str(YEAR), '_TBR_Disag']);
    Folder = (['/projects/The_Env_Virtual_observatory/Sub_Daily_Hydrological_Data/Rainfall_Data/SubDaily_CEH/Grid_data/', Folder]);
    cd (Folder)
    end;
    
    %get the day in
    t = datestr(dates(i,1), 'dd-mmm-yyyy');
    DAY = (['hr_rain_',t,'.mat']);
    load(DAY); %we now have the DAY of data assocaited with the 10 am star time ofinterest
    
    %load each HOUR and save to file after assigning to correct IDs
    
    for HR = 1:1:24;
        rain = Grid_Rain(:,:,HR); rain = rot90(rain);
        
        %for ID = 1:1:length(IDs);
        %for ID = 76021:1:81168; % This was used for VALUES2
        for ID = 1:1:difff;
            IDz = 76021+ID-1;
            grid = rain(IDs(IDz,4):IDs(IDz,5), IDs(IDz,2):IDs(IDz,3));
            VALUES(Row,ID) = nanmean(grid(:));  
        end
        Row = Row +1;
    end;
       
       
end %should by this point have all the data in the big matrix

fprintf('Checkpoint_5\n')

toc

%Next we need to print out this into the text files 
OUT_DATES = [];
     for R = 1:1:length(VALUES(:,1));
        t = datestr(dates(R),'yyyy');
        t1 = datestr(dates(R),'mm');
        t2 =  datestr(dates(R),'dd');  
        t3 = datestr(dates(R), 'HH');
        OUT_DATES(R,1) = str2num(t);
        OUT_DATES(R,2) = str2num(t1);
        OUT_DATES(R,3) = str2num(t2);
        OUT_DATES(R,4) = str2num(t3);
     end

fprintf('Checkpoint_6\n')

cd '/home/pm15021/Hourly_Rainfall/bolita';

for i = 1:1:difff;
    archivo = 76021+i-1;
    Name = (['rainfall_1km_gridID',num2str(archivo),'.txt']);
    fout = fopen (Name,'w');
    for R = 1:1:length(VALUES(:,1));
        t = ([num2str(OUT_DATES(R,1)),' ',  num2str(OUT_DATES(R,2)), ' ', num2str(OUT_DATES(R,3)), ' ', num2str(OUT_DATES(R,4)), ' ', num2str(VALUES(R,i)) ]);
        fprintf(fout, t);
        fprintf(fout, '\n');
    end;
    fclose('all');
    display(Name)
end

%cd ..
dlmwrite('OUT_DATES_Birm.txt',OUT_DATES,'delimiter','\t','newline','pc')
dlmwrite('VALUES_Birm.txt',VALUES,'delimiter','\t','newline','pc')

    
fprintf('Ended\n')