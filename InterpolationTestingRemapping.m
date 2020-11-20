%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Script for testing of accuracy of remapped feedback signals
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Loading feedback data from *.bin-file

filename = 'C:\Users\ZeissLab\Desktop\RemappingFiles\XYpos_2000_19_11_v1.bin';

fid = fopen(filename, 'rb');
data = fread(fid, 'single');
fclose(fid);
size_data = size(data);
data = reshape(data, [(size_data(1)/2), 2]);

%% Philipps Remapping Logic
oldLinSpace = 1:size(data);
newLinSpace = linspace(1, 6150, 61500);
xInt = interp1(oldLinSpace, data(:,1), newLinSpace, 'PCHIP');
% Plotting 
figure; plot(data(:,1)); 
hold on; plot(xInt);

%% Anjas Remapping Logic
x = linspace(1, 6150, 61500);
y = linspace(1, 2, 2);
[X,Y] = meshgrid(x, y);
res = interp2(transpose(XYges), X, Y);
res = transpose(res);
figure; 
plot(res)
