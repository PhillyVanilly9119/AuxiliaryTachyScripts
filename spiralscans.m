clear
clc
close all
%
spotsize=32E-6;  % @ 1/e^2 in m -->spot size = 32µm 
ascanrate=3600000; % in Hz
volumerate=17; %in vol/s
%%%%%%

N_AScans=ascanrate/volumerate; % N scans per volume

dis=0.5*spotsize; %distance between two spots in m

%%

% t = linspace(0,dis*(N_AScans-1),N_AScans); % linspace(x1,x2,n) n points, spacing between points: dis=(x2-x1)/(n-1)
t = linspace(0,1,N_AScans);
% xscanner = sqrt(t) .* cos(2*pi*sqrt(t));
% yscanner = sqrt(t) .* sin(2*pi*sqrt(t));

deltaR=dis; %spiral pitch
Vcv=deltaR*N_AScans; %spiral velocity 

xscanner = sqrt(Vcv*deltaR*t/pi).*cos(sqrt(Vcv*4*pi*t/deltaR)); 
yscanner = sqrt(Vcv*deltaR*t/pi).*sin(sqrt(Vcv*4*pi*t/deltaR));

%%%%
%diameter x diameter y

diameter=[(max(xscanner)-min(xscanner))*1000,(max(yscanner)-min(yscanner))*1000];
%%%%
% xscanner = xscanner / max(xscanner(:)); %x-scanner drive function
% yscanner = yscanner / max(yscanner(:)); %y-scanner drive function






%%

%flyback xscanner
x0 = xscanner(end-1);
x1 = xscanner(end);
x2 = xscanner(1);
x3 = xscanner(2);

v1 = x1 - x0; % Slote of starting points
v2 = x3 - x2; % Slope of end points

n = 99; % number of steps in the flyback pattern
% pattern = zeros(n,1); % Allocate memory for flyback pattern of number of steps = n
 
sampleI = (1:n)';
% v1 = x1 - x0; % why again?
% v2 = x3 - x2;
% A = x2 - x1; % Slope of start and end point
% B = v2 - v1; % Difference of slopes of start and end points
% A2 = A - (v1 + B/2) * (n+1); % Slope Start-End-Point - (Slope Starting-Points + (Diff. Start- & End-Slopes)/2) * Step Size
xx1 = x1 + (x2-x1)/(n+1)*sampleI; % 2nd Point + Vektor(2nd to last - last point)
xx2 = (((v1 + v2) * (n+1)) - 2*(x2-x1)) /(4*pi) .* sin(2*pi/(n+1)*sampleI); % % ... * sin(2 * pi / n=[1,...,100])
xx3 = (v1 - v2)/(2*pi)*(n+1) * sin(pi/(n+1)*sampleI); % ... * sin(pi/n=[1,...,100]) -> Half Period of sin() in xx2
xx = xx1 + xx2 + xx3;
xxdrive = [xscanner xx' ]; % xscanner and xflyback drive

   
%flyback yscanner
y0 = yscanner(end-1);
y1 = yscanner(end);
y2 = yscanner(1);
y3 = yscanner(2);

v1 = y1 - y0;
v2 = y3 - y2;

n = 99; % sampling points
% pattern = zeros(n,1);
 
sampleI = (1:n)';
% v1 = y1 - y0; % why again?
% v2 = y3 - y2;
% A = y2 - y1;
% B = v2 - v1;
% A2 = A - (v1 + B/2) * (n+1);
yy1 = y1 + (y2-y1)/(n+1)*sampleI; % = tangent of?
yy2 = (((v1 + v2) * (n+1)) - 2*(y2-y1))/(4*pi) .* sin(2*pi/(n+1)*sampleI); % constant/ amp.-magn. sine-wave w/ 100 sample-points in range of 2*pi
yy3 = (v1 - v2)/(2*pi) * (n+1) * sin(pi/(n+1)*sampleI);
yy = yy1 + yy2 + yy3; % combines entries of all 3 tangents -> Make up the y-coords of the flyback motion
yydrive = [yscanner yy' ]; % yscanner and yflyback drive




figure(1);
plot(xscanner,yscanner,'b.')
hold on
plot(xx',yy','m.')
title('in m')
% [max(xscanner)-min(xscanner),max(yscanner)-min(yscanner)]    
%m in mm

xxdrive=xxdrive*1000;
yydrive=yydrive*1000;

figure(2)
plot(xxdrive,yydrive)
title('in mm')
% [max(xxdrive)-min(xxdrive),max(yydrive)-min(yydrive)]    


%0.42 Grad pro mm, zoom factor :1.4:
xxdrive=xxdrive*0.42/1.4;
yydrive=yydrive*0.42/1.4;

figure(3)

plot(xxdrive,yydrive)
title('in degree')
% [max(xxdrive)-min(xxdrive),max(yydrive)-min(yydrive)] 
%1 Volt per 2 degree optical deflection
xxdrive=xxdrive/2;
yydrive=yydrive/2;  

figure(4)

plot(xxdrive,yydrive)    
title('in volt')
% [max(xxdrive)-min(xxdrive),max(yydrive)-min(yydrive)] 
%%

%write voltage in excel sheet
filename=strcat (num2str(ascanrate),'Hz_',num2str(volumerate),'vol_',num2str(spotsize),'m.xlsx'); 
% filename = '3200000Hz_17vol_32m.xlsx';
% header = {'xscanner'; 'yscanner'};
T = table(xxdrive',yydrive', 'VariableNames', { 'xscanner', 'yscanner'});
writetable(T,filename)




%%

% ;  %1/A-scan rate in Hz
ascantime = 1/ascanrate;
[pks,locs] = findpeaks(xscanner);
f = 1 ./ (diff(locs)*ascantime);
f = f/1000;


diffx = diff(xscanner);
diffy = diff(yscanner);

for ii = 1:length(diffx)
    stepsize(ii) = sqrt(diffx(ii)^2 + diffy(ii)^2);
end



figure(5)
line(locs(2:end)/N_AScans,f,'Color','b')
ax1 = gca; % current axes
ax1.XColor = 'b';
ax1.YColor = 'b';
ylabel('frequency [kHz]')
xlabel('fraction of volume complete [a.u.]')


ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');

ax2.XColor = 'r';
ax2.YColor = 'r';

line(locs/N_AScans,pks,'Parent',ax2,'Color','r')
ylabel('normalized amplitude [a.u.]')

amplitudefreqproduct = pks(2:end) .* f;

figure(6)
plot(locs(2:end)/N_AScans,amplitudefreqproduct)
ylabel('amplitude x frequency')
xlabel('fraction of volume complete [a.u.]')

figure(7)
plot(pks(2:end), f)
xlabel('normalized amplitude [a.u.]')
ylabel('frequency [kHz]')







