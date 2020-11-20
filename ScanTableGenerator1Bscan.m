clear
clc

%% ScanTableGenerator for 1 bScan
%% input parameters
outputFolder = '\\samba\p_Zeiss\4D OCT\Imaging modes\Scantables\';
nAscans =900;
spotsize = 23e-3; %in mm
zoom =1;
nFlybackPoints = 200;
ascanrate = 100000;
volumerate = 16;
patternnumber = 101; 

AscansVolume = ascanrate / volumerate;

%% calculate A-scan positions
aScanPositionsXmm = (0:spotsize*0.5:(nAscans-1)*spotsize*0.5);  
aScanPositionsYmm(1:nAscans) = 0;

%% calculate flyback

dT = nFlybackPoints+1;

% flyback xscanner

x0 = aScanPositionsXmm(end-1);
x1 = aScanPositionsXmm(end);
x2 = aScanPositionsXmm(1);
x3 = aScanPositionsXmm(2);

v1 = x1 - x0;
v2 = x3 - x2;

flybackPositionsXmm =zeros(nFlybackPoints,1);



for ii = 1:nFlybackPoints
    xx1 = x1;
    xx2 = (x2-x1)/dT*ii;
    xx3 = (x2-x1)/(2*pi) * sin(2*pi*ii/dT);
    xx4 = (v1+v2) / 2 * (dT/(2*pi)) * sin(2*pi*ii/dT);
    xx5 = (v1-v2) / 2 * (dT/pi) * sin(pi*ii/dT);
    flybackPositionsXmm(ii,1) = xx1 + xx2 - xx3 + xx4 + xx5;
end
aScanPositionsXmm = transpose(aScanPositionsXmm);

%flyback yscanner

flybackPositionsYmm (1:nFlybackPoints) = 0;

flybackPositionsYmm = transpose(flybackPositionsYmm);
aScanPositionsYmm = transpose(aScanPositionsYmm);
%% merge actual scan and flyback & convert to voltages
%scan pattern in mm from center
scanPatternMM(:,1) = cat(1, aScanPositionsXmm, flybackPositionsXmm);
scanPatternMM(:,1) = scanPatternMM(:,1) - max(scanPatternMM(:,1))/2;
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
fovMM = maxRadiusMM;

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


%% add A-scans at zero position

Number0Pos = AscansVolume -nAscans -nFlybackPoints;


for    i = nAscans+nFlybackPoints+1 : AscansVolume
scanTable(i,:)= scanTable (nAscans + nFlybackPoints,:);
end


nFlybackPoints = AscansVolume - nAscans;


scanTable (1:AscansVolume,5:10)=0;



%% save scanTable as .bin file
filename = strcat(outputFolder,num2str(patternnumber), '1bScan_', num2str(nAscans), 'aScans_', num2str(nFlybackPoints),'flyback_',num2str(fovMM),"mm",'.bin');
fid = fopen(filename,'w');
fwrite(fid, single(scanTable), 'single');
fclose(fid);

