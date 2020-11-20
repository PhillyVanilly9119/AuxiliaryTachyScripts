% clear
% clc
%% ScanTableGenerator for spiral scans
%% input parameters
outputFolder = 'D:\1 4D OCT\Scanpattern\patternUpdate_08_26_2020\scannerfeedbackTest\1MHzcalib\';
% outputFolder = 'C:\Users\Philipp\Desktop\';

zoom = 1;
spotsize = 23E-6;  % @ 1/e^2 in m
sweepRate = 1000000; % in Hz
spectralSplittingFactor = 1;
volRate = 16; %in vol/s
nFlybackPoints = 450;
removeInnerNPoints = 1000;
patternnumber = '2000_19_11_v1';

%% pre-calculations
ascanRate = sweepRate * spectralSplittingFactor;
nAscans = floor(ascanRate/volRate); % N scans per volume
nAscans = nAscans-nFlybackPoints; % N scans per volume flyback points subtracted
dutyCycle = nAscans / (nAscans + nFlybackPoints) *100; 
deltaR = 0.5*spotsize; %distance between two spots in m

%% calculate A-scan positions
t = linspace(0,1,nAscans)';
Vcv = deltaR*nAscans; %spiral velocity 

aScanPositionsXmm = sqrt(Vcv*deltaR*t/pi).*cos(sqrt(Vcv*4*pi*t/deltaR)) * 1000; 
aScanPositionsYmm = sqrt(Vcv*deltaR*t/pi).*sin(sqrt(Vcv*4*pi*t/deltaR)) * 1000;

aScanPositionsXmm(1:removeInnerNPoints) = [];
aScanPositionsYmm(1:removeInnerNPoints) = [];
nAscans = nAscans-removeInnerNPoints;
 
% %%flip spiral (scanners start from the outside)
% aScanPositionsXmm = flip (aScanPositionsXmm);
% aScanPositionsYmm = flip (aScanPositionsYmm);

%% calculate flyback
dT = nFlybackPoints+1;

%flyback xscanner
x0 = aScanPositionsXmm(end-1);
x1 = aScanPositionsXmm(end);
x2 = aScanPositionsXmm(1);
x3 = aScanPositionsXmm(2);

v1 = x1 - x0;
v2 = x3 - x2;

flybackPositionsXmm = zeros(nFlybackPoints,1);
for ii = 1:nFlybackPoints
    xx1 = x1;
    xx2 = (x2-x1)/dT*ii;
    xx3 = (x2-x1)/(2*pi) * sin(2*pi*ii/dT);
    xx4 = (v1+v2) / 2 * (dT/(2*pi)) * sin(2*pi*ii/dT);
    xx5 = (v1-v2) / 2 * (dT/pi) * sin(pi*ii/dT);
    flybackPositionsXmm(ii,1) = xx1 + xx2 - xx3 + xx4 + xx5;
end

%flyback yscanner
y0 = aScanPositionsYmm(end-1);
y1 = aScanPositionsYmm(end);
y2 = aScanPositionsYmm(1);
y3 = aScanPositionsYmm(2);

v1 = y1 - y0;
v2 = y3 - y2;

flybackPositionsYmm = zeros(nFlybackPoints,1);
for ii = 1:nFlybackPoints
    yy1 = y1;
    yy2 = (y2-y1)/dT*ii;
    yy3 = (y2-y1)/(2*pi) * sin(2*pi*ii/dT);
    yy4 = (v1+v2) / 2 * (dT/(2*pi)) * sin(2*pi*ii/dT);
    yy5 = (v1-v2) / 2 * (dT/pi) * sin(pi*ii/dT);
    flybackPositionsYmm(ii,1) = yy1 + yy2 - yy3 + yy4 + yy5;
end


%% debug
testX = cat(1, aScanPositionsXmm, flybackPositionsXmm, aScanPositionsXmm);
diff1X = diff(testX*10,1);
diff2X = diff(testX*10,2);

testY = cat(1, aScanPositionsYmm, flybackPositionsYmm, aScanPositionsYmm);
diff1Y = diff(testY*10,1);
diff2Y = diff(testY*10,2);

figure(101)
plot(aScanPositionsXmm, aScanPositionsYmm)
hold on
plot(flybackPositionsXmm, flybackPositionsYmm)
hold off

figure(102)
plot(1:length(aScanPositionsXmm), aScanPositionsXmm, '*k')
hold on
plot(nAscans+1:nAscans+nFlybackPoints, flybackPositionsXmm, '*k')
plot(nAscans+nFlybackPoints+1:2*nAscans+nFlybackPoints, aScanPositionsXmm, '*k')
plot(2:nAscans*2+nFlybackPoints, diff1X, 'b')
plot(3:nAscans*2+nFlybackPoints, diff2X, 'r')
hold off
ylim([-0.25 0.25])
xlim([70550 70650])

figure(103)
plot(1:length(aScanPositionsYmm), aScanPositionsYmm, '*k')
hold on
plot(nAscans+1:nAscans+nFlybackPoints, flybackPositionsYmm, '*k')
plot(nAscans+nFlybackPoints+1:2*nAscans+nFlybackPoints, aScanPositionsYmm, '*k')
plot(2:nAscans*2+nFlybackPoints, diff1Y, 'b')
plot(3:nAscans*2+nFlybackPoints, diff2Y, 'r')
hold off
ylim([-0.25 0.25])
xlim([70550 70650])

%% merge actual scan and flyback & convert to voltages
%scan pattern in mm from center
scanPatternMM(:,1) = cat(1, aScanPositionsXmm, flybackPositionsXmm);
scanPatternMM(:,2) = cat(1, aScanPositionsYmm, flybackPositionsYmm);

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
maxRadiusMM = max([max(abs(aScanPositionsXmm)) max(abs(aScanPositionsYmm))]);
fovMM = maxRadiusMM * 2;

rectGridSize = ceil(fovMM / (deltaR*1000));

%% calculate original A-scan positions on rect grid
aScanPositionsXrectGridInd = (aScanPositionsXmm + maxRadiusMM) / fovMM * (rectGridSize-1) +1;
aScanPositionsYrectGridInd = (aScanPositionsYmm + maxRadiusMM) / fovMM * (rectGridSize-1) +1;

scanTable(1:nAscans,5) = aScanPositionsXrectGridInd;
scanTable(1:nAscans,6) = aScanPositionsYrectGridInd;


nAscans = 1000;
%% find destinations on rect grid for each A-scan
[squareGrid(:,1), squareGrid(:,2)] = ind2sub([rectGridSize rectGridSize], 1:rectGridSize^2);

T = delaunayn(scanTable(1:nAscans,5:6));
k = dsearchn(scanTable(1:nAscans,5:6), T, squareGrid,Inf);  


for ii = 1:nAscans
    destinations = find(k==ii);
    if ~isempty(destinations)
        for jj = 1:length(destinations)
            [scanTable(ii,6+jj*2-1), scanTable(ii,6+jj*2)] = ind2sub([rectGridSize rectGridSize], destinations(jj));
        end
    end
end

%% plots
figure(1)
plot(scanPatternMM(1:nAscans,1), scanPatternMM(1:nAscans,2), ".b")
hold on
plot(scanPatternMM(nAscans+1:end,1), scanPatternMM(nAscans+1:end,2), ".r")
hold off

figure(2)
plot(squareGrid(:,1), squareGrid(:,2), 'o')
hold on
for ii = 5%:2:size(scanTable,2)
    plot(scanTable(:,ii), scanTable(:,ii+1), '.')
end
hold off
xlim([round(rectGridSize/2)-5 round(rectGridSize/2)+5])
ylim([round(rectGridSize/2)-5 round(rectGridSize/2)+5])


%%
%check remap file 
figure;
plot(scanTable(1:nAscans,7),scanTable(1:nAscans,8),'.')
hold on
plot(scanTable(1:nAscans,9),scanTable(1:nAscans,10),'.')


%% save scanTable as .bin file

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% if saved .bin file for engine consists of less than 10 columns,
% readbin.vi needs to be adjusted
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






filename = strcat(outputFolder,num2str(patternnumber), '_Spiral_', num2str(spectralSplittingFactor), 'x_', num2str(sweepRate/1000),'kHz_', num2str(volRate),"vol_",num2str(fovMM),"mm_",num2str(nAscans),'AScans_',num2str(nFlybackPoints),'flybackpoints','.bin');
fid = fopen(filename,'w');



fwrite(fid, single(scanTable(:,1:10)), 'single');
fclose(fid);

%% save scanTable as .bin file for Philipp
% filename = strcat(outputFolder, num2str(patternnumber),'_Spiral.bin');
filename = strcat(outputFolder,'donuttest.bin');
fid = fopen(filename,'w');

fwrite(fid, scanTable(:,7:10), 'uint32');

fclose(fid);