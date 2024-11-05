% This script implements a semi-empirical method, based on the ratio transform
% algorithm first proposed by Stumpf et al. (2003), to retrieve bathymetric data
% from multispectral optical imagery.

% - Mind that all loaded files must be located in the same folder ("Current Folder" where MATLAB operates)
% - Mind not to overwrite the results by changing the file names before any new run

close all; clear; clc   % Clear workspace and command window


%% RANDOM SEED and GENERATOR

seed = 1;
rng(seed,"twister");


%% INPUT DATA AND PARAMETERS

% Any input data and parameters which have to be manually specified by the
% user are listed in this section of the script.

% Read the preprocessed multispectral raster image
[A, R] = readgeoraster('Orthophoto.tif');   % Rename with the proper file name of the multispectral orthophoto

% RGB bands to matrices correctly oriented
Red = flip(A(:,:,3));
Green = flip(A(:,:,2));
Blue = flip(A(:,:,1));

Res = R.CellExtentInWorldX;   % Spatial resolution of the raster bands (in meters)

% Read the vector file containing depth values (in meters) to use as ground
% truth points (GTPs) for the calibration and validation of the model
V = readgeotable('BathymetryPoints.shp');   % Rename with the proper file name of the ground truth bathymetric points

% Parameters for ground truth data filtering
UD = 10;   % Upper depth boundary [m]
LD = 0;    % Lower depth boundary [m]
P = 0;     % Percentage of the dataset to discard (set to 0 if you wish to work with all of the available data)

% Method to be used to assign a single depth value to each pixel of the raster image
Method = 'median';   % Options list: 'mean' / 'median'

% Ground truth points split between calibration and validation sets
P_Cal = 0.5;           % Percentage of GTPs to use for the calibration of the model
% P_Val = 1 - P_Cal;   % Percentage of GTPs to use for the validation of the model

% Reference spectral band ratio
Band_Ratio = 'Blue_Red';   % Options list: 'Blue_Red' / 'Blue_Green'


%% EPSG CODE RETRIEVAL

% Mind that the CRS of the raster image and the vector GTPs dataset must be the same

CRSstr= R.ProjectedCRS.Name;
zone = extractBetween(CRSstr, 19, 20);
emi = extractAfter(CRSstr, 20);

if emi == "S"
   rootepsg = "327"; 
else
	 rootepsg = "326";
end

epsg = strcat(rootepsg, zone);


%% FILTER GROUND TRUTH POINTS

V = V(~any(ismissing(V),2),:);   % Remove any NaN values

% Set depth range
v(:,1)  = V{1:end,2};   
vind = find(v>UD);        % Set upper depth boundary

for i = transpose(vind)
    v(i) = NaN;
end

clear('vind')
vind = find(v<LD);        % Set lower depth boundary

for i = transpose(vind)
    v(i) = NaN;
end

% Randomly discard a set percentage (P%) of the initial dataset
vsize = numel(v);
idx = transpose(randperm(vsize));

for i = 1:length(idx)
    if i <= length(idx)*P
       idx(i) = idx(i);
    elseif i > length(idx)*P
           idx(i) = NaN;
    end
end

idx(any(isnan(idx),2),:) = []; 

for i = idx
    v(i) = NaN;
end

V{1:end, 2} = v;
V = V(~any(ismissing(V),2),:);


%% RETRIEVE GROUND TRUTH POINTS COORDINATES

% Retrieve GTPs coordinates from "V"
Lat = V.Shape.Y(:,1);
Lon = V.Shape.X(:,1);
Coordinates1 = [(Lon+Res)-R.XWorldLimits(1), (Lat+Res)-R.YWorldLimits(1)];
Coordinates2 = Coordinates1./Res; 
X = floor(Coordinates2(:,1));
Y = floor(Coordinates2(:,2));

% Uncomment to display raster image with ground truth points
% pcolor(Red)
% hold on
% shading flat
% plot(Coordinates2(:,1), Coordinates2(:,2), 'or', 'MarkerSize', 4, 'markerfacecolor', 'r')
% plot(X, Y, '*g')
% axis equal
% xlim([min(X)-5 max(X)+5])
% ylim([min(Y)-5 max(Y)+5])


%% CREATE AND FILL A MATRIX OF THE PROPER SIZE WITH RGB AND DEPTH VALUES

s = length(X);
Mat_R = zeros(s,1);
Mat_G = zeros(s,1);
Mat_B = zeros(s,1);

for j = 1:s
    b =  X(j);
    a =  Y(j);
    Mat_R(j) = Red(a,b);
    Mat_G(j) = Green(a,b);
    Mat_B(j) = Blue(a,b);
end

Depth = table2array(V(:,2));

M = [X Y Mat_R Mat_G Mat_B Depth];


%% DISCARD GROUND TRUTH POINTS THAT FALL WITHIN RASTER CELLS WHERE ONE OR MORE SPECTRAL BANDS ARE NULL

M(any(M(:,3:5)==0, 2), :) = [];
Mat_R = M(:,3);
Mat_G = M(:,4);
Mat_B = M(:,5);
Depth = M(:,6);


%% ASSOCIATE A SINGLE DEPTH VALUE TO EACH RASTER CELL

[UniqueCoordinates, ind, IDs] = unique(M(:,1:2), 'rows', 'stable');
[aux, Points] = findgroups(IDs(:,1));

if strcmp(Method, 'mean') == 1 
   Depth2 = [Points splitapply(@mean, M(:,6), aux)];
elseif strcmp(Method, 'median') == 1
       Depth2 = [Points splitapply(@median, M(:,6), aux)];
end

clear('M')
M = [UniqueCoordinates Mat_R(ind) Mat_G(ind) Mat_B(ind) Depth2(:,2)];
Mat_R = M(:,3);
Mat_G = M(:,4);
Mat_B = M(:,5);
Depth = M(:,6);


%% COMPUTE SPECTRAL BANDS RATIOS (ACCORDING TO STUMPF 2003)

Ratio_BR = log(1000*Mat_B)./log(1000*Mat_R);
Ratio_BG = log(1000*Mat_B)./log(1000*Mat_G);

Map_ratio_BR = log(double(Blue)*1000)./ log(double(Red)*1000);
Map_ratio_BG = log(double(Blue)*1000)./ log(double(Green)*1000);

% Find and remove any NaN values
N = isnan(Ratio_BR);
Depth = Depth(~N);
Ratio_BR = Ratio_BR(~N);
Ratio_BG = Ratio_BG(~N);
clear('N')
N = isnan(Ratio_BG);
Depth = Depth(~N);
Ratio_BR = Ratio_BR(~N);
Ratio_BG = Ratio_BG(~N);

Sampled_Points = [Depth, Ratio_BR, Ratio_BG];


%% SPLIT GTPs BETWEEN CALIBRATION (TRAINING) AND VALIDATION (TESTING) SETS

[m, n] = size(Sampled_Points);
idx = randperm(m);
Training = Sampled_Points(idx(1:round(P_Cal*m)),:); 
Testing = Sampled_Points(idx(round(P_Cal*m)+1:end),:);

Depth_Cal = Training(:,1);
Depth_Val = Testing(:,1);

Ratio_BR_Cal = Training(:,2);
Ratio_BR_Val = Testing(:,2);

Ratio_BG_Cal = Training(:,3);
Ratio_BG_Val = Testing(:,3);


%% CALIBRATION - EXPONENTIAL REGRESSION

if strcmp(Band_Ratio, 'Blue_Red') == 1
   x = Ratio_BR_Cal;
elseif strcmp(Band_Ratio, 'Blue_Green') == 1
       x = Ratio_BG_Cal;
end

y = Depth_Cal;

% Logarithmic regression (unfiltered data)
ft1 = fittype('exp(a*x+b)+c', 'indep', 'x');   % Exponential regression: y = exp(a*x+b)+c
options1 = fitoptions('Method', 'NonlinearLeastSquares', 'Upper', [Inf Inf -1]); 
[f1, S1] = fit(x, y, ft1, options1);

% Calibration parameters of the exponential regression (unfiltered data)
a = f1.a;
b = f1.b;
c = f1.c;
yf1 = exp(a*x+b)+c;       % Fitted model
R_squared = S1.rsquare;   % R^2
N_Cal = length(x);        % Number of GTPs used for the model calibration
                          % (before the removal of the outliers)

% Outliers removal (according to the 3-sigma criterion)
I = abs(yf1 - y) > 3*std(y); 
outliers = excludedata(x,y,'indices',I);
Y = y(~outliers);
X = x(~outliers);

% Exponential regression (filtered data)
ft2 = fittype('exp(a*x+b)+c', 'indep', 'x');   % Exponential regression: Y = exp(a*X+b)+c
options2 = fitoptions('Method', 'NonlinearLeastSquares', 'Upper', [Inf Inf -1], 'StartPoint', [a b c]);
[f2, S2] = fit(X, Y, ft2, options2);

% Calibration parameters of the exponential regression (filtered data)
a_clean = f2.a;
b_clean = f2.b;
c_clean = f2.c;
R_squared_clean = S2.rsquare;   % R^2 (clean)
N_Cal_clean = length(X);        % Number of GTPs used for the model calibration
                                % (after the removal of the outliers)


%% VALIDATION - EXPONENTIAL REGRESSION

% Application of the exponential model

if strcmp(Band_Ratio, 'Blue_Red') == 1
   z = Ratio_BR_Val;
elseif strcmp(Band_Ratio, 'Blue_Green') == 1
       z = Ratio_BG_Val;
end

x = Depth_Val;
y = exp(a_clean*z+b_clean)+c_clean;     % RDB (estimated bathymetry)
RDB = y;

p = polyfit(x,y,1);   % Linear interpolation
f = polyval(p,x);

% Validation stats
N_Val = length(x);              % Number of GTPs used for the model validation
BIAS = x - y;
BIASP = BIAS > 10;
BIAS = BIAS(~BIASP);
BIASN = BIAS < -10;
BIAS = BIAS(~BIASN);
RMSE = sqrt(mean((BIAS).^2));   % Root Mean Square Error
MAE = mean(abs(BIAS));          % Mean Average Error
BIAS_AV = mean(BIAS);           % Average BIAS
BIAS_STD = std(BIAS);           % BIAS Standard Deviation
rmse = sprintf('RMSE = %0.3f\n', RMSE);
n_val = sprintf('N Validation = %0.0f\n', N_Val);
biasav = sprintf('BIAS AV = %0.3f\n', BIAS_AV);
mae = sprintf('MAE = %.3f', MAE);
biastd = sprintf('BIAS STD = %.3f', BIAS_STD );


%% RDB MAP

if strcmp(Band_Ratio, 'Blue_Red') == 1
   z = Map_ratio_BR;
elseif strcmp(Band_Ratio, 'Blue_Green') == 1
       z = Map_ratio_BG;
end

Map_RDB = flip(exp(a_clean*z+b_clean)+c_clean);

% Set depth range
idx_lower = Map_RDB < LD;
Map_RDB(idx_lower) = LD;
idx_upper = Map_RDB > UD;
Map_RDB(idx_upper) = UD;

% idx_NaN = isnan(Map_RDB);
% Map_RDB(idx_NaN) = 0;

figure(1)
imshow(Map_RDB, [LD UD])
hold on
colormap turbo % parula % jet
colorbar eastoutside

% Save RDB map in .tiff format 
Epsg = strcat('EPSG:', strcat(rootepsg, zone));
geotiffwrite('Map_SDB_exp', Map_RDB, R, CoordRefSysCode = Epsg);


%% CALIBRATION AND VALIDATION PLOTS

figure(2)

subplot(1, 5, [1 2 3])
if strcmp(Band_Ratio, 'Blue_Red') == 1
   x = Ratio_BR_Cal;
elseif strcmp(Band_Ratio, 'Blue_Green') == 1
       x = Ratio_BG_Cal;
end
y = Depth_Cal;
x_aux = linspace(0, 2*max(x), 1000);
plot(x, y, '*', 'color', [.2 .4 1])
hold on
plot(X, Y, 'w*', 'HandleVisibility','off')
scatterdensity(X, Y, 'filled', 'MarkerSize', 50)
colormap parula
if strcmp(Band_Ratio, 'Blue_Red') == 1
   plot(x_aux, f2(x_aux), 'color', [.9 0 0], 'linewidth', 1.7)
elseif strcmp(Band_Ratio, 'Blue_Green') == 1
       plot(x_aux, f2(x_aux), 'color', '#33952B', 'linewidth', 1.7)
end
plot(x_aux, f1(x_aux)+3*std(y),  '--', 'color', [.7 .7 .7], 'linewidth', 1.7, 'HandleVisibility','off')
plot(x_aux, f1(x_aux)-3*std(y),  '--', 'color', [.7 .7 .7], 'linewidth', 1.7, 'HandleVisibility','off')
grid on
xlim([0.5*floor(min(x)/0.5) 0.5*ceil(max(x)/0.5)])
ylim([0 5*ceil(max(y)/5)])
xticks(0:0.25:0.5*ceil(max(x)/0.5))
yticks(0:1:5*ceil(max(y)/5))
if strcmp(Band_Ratio, 'Blue_Red') == 1
   title('Calibration - Exponential regression')
   xlabel('Blue/Red ratio')
elseif strcmp(Band_Ratio, 'Blue_Green') == 1
       title('Calibration - Exponential regression')
       xlabel('Blue/Green ratio')
end
ylabel('Depth (m)')
legend('Data outliers', 'Data', 'Exponential regression', 'Location', 'northwest', 'Fontsize', 10, 'Linewidth', 1)
dim = [.138 .64 .1 .1];
str1 = ['Exp. regression:  y = exp(a*x+b)+c'];
str2 = ['a = ', num2str(a_clean, '%0.2f'), '     b = ', num2str(b_clean, '%0.2f'),  '     c = ', num2str(c_clean,'%0.2f')];
str3 = ['N Calibration = ', num2str(N_Cal_clean, '%0.0f')];
str4 = ['R^2 = ', num2str(R_squared_clean, '%0.3f')];
str = str1 + "\n" + str2 + "\n" + str3 + "\n" + str4;
str = compose(str);
t = annotation('textbox', dim, 'String', str, 'Fontsize', 10, 'Linewidth', 0.01);
t.BackgroundColor = [1 1 1];
t.EdgeColor = [1 1 1];
t.FaceAlpha = 0;


subplot(1,5,[4 5])
x = Depth_Val;
y = RDB;
scatterdensity(x, y, 'filled', 'MarkerSize', 50)
hold on
bisector = linspace(LD, 100, 100);
plot(bisector, bisector,'k--', 'linewidth', 1.2)
if strcmp(Band_Ratio, 'Blue_Red') == 1
   plot(x, f, '-', 'color', [.9 0 0], 'linewidth', 1.7)
   title('Validation - Exponential regression')
elseif strcmp(Band_Ratio, 'Blue_Green') == 1
       plot(x, f, '-', 'color', '#33952B', 'linewidth', 1.7)
       title('Validation - Exponential regression')
end
grid on
xlabel('Control depth (m)')
ylabel('Estimated depth (m)')
legend('Data', 'Bisector', 'Linear interpolation', 'Location', 'northwest', 'Fontsize', 10, 'Linewidth', 1)
str1 = ['N Validation = ',  num2str(N_Val, '%0.0f\n')];
str2 = ['RMSE = ', num2str(RMSE, '%0.3f\n')];
str3 = ['MAE = ', num2str(MAE, '%0.3f\n')];
str4 = ['BIAS AV = ', num2str(BIAS_AV, '%0.3f\n')];
str5 = ['BIAS STD = ', num2str(BIAS_STD, '%0.3f\n')];
str = str1 + "\n" + str2 + "\n" + str3 + "\n" + str4 + "\n" + str5;
str = compose(str);
dim = [.773 .272 .1 .1];
t = annotation('textbox', dim, 'String', str, 'Fontsize', 10, 'Linewidth', 0.01);
t.BackgroundColor = [1 1 1];
t.EdgeColor = [1 1 1];
t.FaceAlpha = 0;
axis equal
xlim([LD ceil(max(x))]) 
ylim([LD ceil(max(x))]) 
xticks(LD:1:ceil(max(x)))
yticks(LD:1:ceil(max(x)))
box on

x0 = 50;
y0 = 200;
width = 1300;
height = 400;
set(gcf, 'position', [x0, y0, width, height])

saveas(gcf, 'CAL_VAL_exp.png') 




%% Scatter density plot function

function h = scatterdensity(x, y, varargin)
c = ksdensity([x,y], [x,y]);
if nargin > 2
  % Set Marker Size
  i = find(strcmp(varargin, 'MarkerSize'),1);
  if ~isempty(i); MarkerSize = varargin{i+1}; varargin(i:i+1) = [];
  else MarkerSize = []; end
  % Plot scatter plot
  h = scatter(x, y, MarkerSize, c, varargin{:});
else
  h = scatter(x, y, [], c);
end
end
