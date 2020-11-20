clear
clc
%% ScanTableGenerator for Bscan crosses
%% input parameters
outputFolder = '\\samba\p_Zeiss\4D OCT\Imaging modes\Scantables\';
nAscans =1024;
spotsize = 23e-3; %in mm
zoom =1;
nFlybackPoints = 200;
patternnumber = 101;

%% calculate A-scan positions for cross Bscan
aScanPositionsBscan1Xmm = (0:spotsize*0.5:(nAscans-1)*spotsize*0.5);  


aScanPositionsBscan2Ymm = rot90(aScanPositionsBscan1Xmm,2);
aScanPositionsBscan2Ymm = aScanPositionsBscan2Ymm - max(aScanPositionsBscan2Ymm)/2;

aScanPositionsBscan1Ymm(1:nAscans) = 0;

aScanPositionsBscan2Xmm(1:nAscans) = aScanPositionsBscan1Xmm(nAscans/2);


%% calculate flyback

dT = nFlybackPoints+1;
%% xscanner
% flyback xscanner Bscan1Bscan2

x0 = aScanPositionsBscan1Xmm(end-1);
x1 = aScanPositionsBscan1Xmm(end);
x2 = aScanPositionsBscan2Xmm(1);
x3 = aScanPositionsBscan2Xmm(2);

v1 = x1 - x0;
v2 = x3 - x2;

flybackPositionsBscan1Bscan2Xmm =zeros(nFlybackPoints,1);



for ii = 1:nFlybackPoints
    xx1 = x1;
    xx2 = (x2-x1)/dT*ii;
    xx3 = (x2-x1)/(2*pi) * sin(2*pi*ii/dT);
    xx4 = (v1+v2) / 2 * (dT/(2*pi)) * sin(2*pi*ii/dT);
    xx5 = (v1-v2) / 2 * (dT/pi) * sin(pi*ii/dT);
    flybackPositionsBscan1Bscan2Xmm(ii,1) = xx1 + xx2 - xx3 + xx4 + xx5;
end


% flyback xscanner Bscan2Bscan1


x0 = aScanPositionsBscan2Xmm(end-1);
x1 = aScanPositionsBscan2Xmm(end);
x2 = aScanPositionsBscan1Xmm(1);
x3 = aScanPositionsBscan1Xmm(2);

v1 = x1 - x0;
v2 = x3 - x2;

flybackPositionsBscan2Bscan1Xmm =zeros(nFlybackPoints,1);



for ii = 1:nFlybackPoints
    xx1 = x1;
    xx2 = (x2-x1)/dT*ii;
    xx3 = (x2-x1)/(2*pi) * sin(2*pi*ii/dT);
    xx4 = (v1+v2) / 2 * (dT/(2*pi)) * sin(2*pi*ii/dT);
    xx5 = (v1-v2) / 2 * (dT/pi) * sin(pi*ii/dT);
    flybackPositionsBscan2Bscan1Xmm(ii,1) = xx1 + xx2 - xx3 + xx4 + xx5;
end



% Merge scanpattern x scanner

aScanPositionsBscan1Xmm = transpose(aScanPositionsBscan1Xmm);
aScanPositionsBscan2Xmm = transpose(aScanPositionsBscan2Xmm);

scanPatternMM(:,1) = cat(1, aScanPositionsBscan1Xmm, flybackPositionsBscan1Bscan2Xmm, aScanPositionsBscan2Xmm, flybackPositionsBscan2Bscan1Xmm);
scanPatternMM(:,1) = scanPatternMM(:,1) - max (scanPatternMM(:,1))/2;
%% yscanner
% flyback yscanner Bscan1Bscan2

x0 = aScanPositionsBscan1Ymm(end-1);
x1 = aScanPositionsBscan1Ymm(end);
x2 = aScanPositionsBscan2Ymm(1);
x3 = aScanPositionsBscan2Ymm(2);

v1 = x1 - x0;
v2 = x3 - x2;

flybackPositionsBscan1Bscan2Ymm =zeros(nFlybackPoints,1);



for ii = 1:nFlybackPoints
    xx1 = x1;
    xx2 = (x2-x1)/dT*ii;
    xx3 = (x2-x1)/(2*pi) * sin(2*pi*ii/dT);
    xx4 = (v1+v2) / 2 * (dT/(2*pi)) * sin(2*pi*ii/dT);
    xx5 = (v1-v2) / 2 * (dT/pi) * sin(pi*ii/dT);
    flybackPositionsBscan1Bscan2Ymm(ii,1) = xx1 + xx2 - xx3 + xx4 + xx5;
end


% flyback xscanner Bscan2Bscan1


x0 = aScanPositionsBscan2Ymm(end-1);
x1 = aScanPositionsBscan2Ymm(end);
x2 = aScanPositionsBscan1Ymm(1);
x3 = aScanPositionsBscan1Ymm(2);

v1 = x1 - x0;
v2 = x3 - x2;

flybackPositionsBscan2Bscan1Ymm =zeros(nFlybackPoints,1);



for ii = 1:nFlybackPoints
    xx1 = x1;
    xx2 = (x2-x1)/dT*ii;
    xx3 = (x2-x1)/(2*pi) * sin(2*pi*ii/dT);
    xx4 = (v1+v2) / 2 * (dT/(2*pi)) * sin(2*pi*ii/dT);
    xx5 = (v1-v2) / 2 * (dT/pi) * sin(pi*ii/dT);
    flybackPositionsBscan2Bscan1Ymm(ii,1) = xx1 + xx2 - xx3 + xx4 + xx5;
end



% Merge scanpattern y scanner

aScanPositionsBscan1Ymm = transpose(aScanPositionsBscan1Ymm);
aScanPositionsBscan2Ymm = transpose(aScanPositionsBscan2Ymm);

scanPatternMM(:,2) = cat(1, aScanPositionsBscan1Ymm, flybackPositionsBscan1Bscan2Ymm, aScanPositionsBscan2Ymm, flybackPositionsBscan2Bscan1Ymm);





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
figure(1)
plot(scanPatternV(:,1))
hold on
plot(scanPatternV(:,2))
hold off

%% scantable needs to be 10 columns longs
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% if saved .bin file for engine consists of less than 10 columns,
% readbin.vi needs to be adjusted
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
scanTable (1:(nAscans+nFlybackPoints),5:10)=0;

%% save scanTable as .bin file
filename = strcat(outputFolder, num2str(patternnumber), 'BscanCross_', num2str(nAscans), 'aScans_', num2str(nFlybackPoints),'flyback_',num2str(fovMM),"mm",'.bin');
fid = fopen(filename,'w');
fwrite(fid, single(scanTable), 'single');
fclose(fid);