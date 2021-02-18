clear
clc
%% ScanTableGenerator for Raster Capture Scan
clear all
close all
%% input parameters
outputFolder = '\\samba\p_Zeiss\4D OCT\Imaging modes\Scantables\';
addpath('D:\1 4D OCT\Matlab Functions');
nAscans = 512;
nBscans = 512;
spotsize = 23e-3; %in mm
zoom =1;
nFlybackPointsBscan = 50;
nFlybackPointsLastFirstBscan = 250;
bidirectional = 1; %1=true, 0=false
patternnumber = '20';

%% calculate A-scan positions for capture Bscan

ScanPositionsMMXApos(1:nAscans) = (0:spotsize*0.5:(nAscans-1)*spotsize*0.5);
ScanPositionsMMYApos(1:nAscans) = 0;

for i = 2:1:nBscans
    ScanPositionsMMXApos = [ScanPositionsMMXApos;ScanPositionsMMXApos(1,1:nAscans)];
    ScanPositionsMMYApos = [ScanPositionsMMYApos;ScanPositionsMMYApos(1,1:nAscans)+(i-1)*spotsize*0.5];
end

ScanPositionsMMXApos = transpose(ScanPositionsMMXApos);
ScanPositionsMMYApos = transpose(ScanPositionsMMYApos);

if bidirectional == 1
    for i = 2:2:nBscans
        ScanPositionsMMXApos(:,i) = flip(ScanPositionsMMXApos(:,i));
        ScanPositionsMMYApos(:,i) = flip(ScanPositionsMMYApos(:,i));
    end
end
%% calculate flyback between B-scans
ScanPositionsMMX = ScanPositionsMMXApos;
ScanPositionsMMY = ScanPositionsMMYApos;

ScannerDriveMMX = [];
ScannerDriveMMY = [];

for i = 1:nBscans-1
flybackX = F_CreateScannerFlyback(nFlybackPointsBscan,ScanPositionsMMX(nAscans,i), ScanPositionsMMX(nAscans-1,i), ScanPositionsMMX(1,i+1),ScanPositionsMMX(2,i+1));
flybackY = F_CreateScannerFlyback(nFlybackPointsBscan,ScanPositionsMMY(nAscans,i), ScanPositionsMMY(nAscans-1,i), ScanPositionsMMY(1,i+1),ScanPositionsMMY(2,i+1));
ScannerDriveMMX = [ScannerDriveMMX; squeeze(ScanPositionsMMX(:,i)); transpose(flybackX)];
ScannerDriveMMY = [ScannerDriveMMY; squeeze(ScanPositionsMMY(:,i)); transpose(flybackY)];
end
deltaflybackX = max(abs(diff(flybackX)));
deltaflybackY = max(abs(diff(flybackY)));

%% Calculate Flyback between last and first B-scan
lastflybackX = F_CreateScannerFlyback(nFlybackPointsLastFirstBscan,ScanPositionsMMX(nAscans,nBscans), ScanPositionsMMX(nAscans-1,nBscans), ScanPositionsMMX(1,1),ScanPositionsMMX(2,1));
lastflybackY = F_CreateScannerFlyback(nFlybackPointsLastFirstBscan,ScanPositionsMMY(nAscans,nBscans), ScanPositionsMMY(nAscans-1,nBscans), ScanPositionsMMY(1,1),ScanPositionsMMY(2,1));
deltalastflybackX = max(abs(diff (lastflybackX)));
deltalastflybackY = max(abs(diff (lastflybackY)));

ScannerDriveMMX = [ScannerDriveMMX; transpose(lastflybackX)];
ScannerDriveMMY = [ScannerDriveMMY; transpose(lastflybackY)];

figure; 
plot(ScannerDriveMMX,ScannerDriveMMY)
hold on
plot(ScanPositionsMMXApos,ScanPositionsMMYApos,'*')



%% shift middle to zero
scanPatternMM = [ScannerDriveMMX, ScannerDriveMMY];
max_position = max(scanPatternMM)*0.5;
scanPatternMM =scanPatternMM - max_position;

figure; plot(scanPatternMM(:,1),scanPatternMM(:,2))

%% merge actual scan and flyback & convert to voltages
%scan pattern in mm from center is scanPatternMM

%scan pattern in deg mirror deflection (0.42 Grad pro mm) CHECK CONVERSION FACTOR (MECHANICAL OR OPTICAL DEGREE?)!!!!
scanPatternDeg = scanPatternMM * 0.42 * zoom;

%scan pattern in V (1V per degree) CHECK CONVERSION FACTOR (MECHANICAL OR OPTICAL DEGREE?)!!!!
scanPatternV = scanPatternDeg * 1;

scanTable(:,1:2) = scanPatternV;
scanTable(:,3:4) = scanPatternMM;

%% check scan pattern
% check max voltage
if max(abs(scanPatternV(:))) > 10
    disp('Scan pattern exceeds max voltage! ScanTable was not created.')
    return
end

%% calculate field of view & size of rectangular grid
maxRadiusMM = max([max(abs(scanPatternMM(:,1))) max(abs(scanPatternMM(:,2)))]);
fovMM = maxRadiusMM * 2;

%% plots
figure;
plot(scanPatternV(:,1))
hold on
plot(scanPatternV(:,2))
hold off

%% scantable needs to be 10 columns longs
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% if saved .bin file for engine consists of less than 10 columns,
% readbin.vi needs to be adjusted
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
scanTable (:,5:10)=0;

%% save scanTable as .bin file
filename = strcat(outputFolder, num2str(patternnumber), 'BscanRaster_', num2str(nAscans), 'aScans_', num2str(nBscans), 'bScans_', num2str(nFlybackPointsBscan),'flybackBscan_', num2str(nFlybackPointsLastFirstBscan),'flybackLastFirstBscan_',num2str(fovMM),"mm",'.bin');
fid = fopen(filename,'w');
fwrite(fid, single(scanTable), 'single');
fclose(fid);