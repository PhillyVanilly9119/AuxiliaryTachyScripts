function [drivesignal] = F_CreateScannerFlyback(nFlybackPoints,scannerEnd, scannerEndm1, scannerStart,scannerStartp1 )

%calculate scanner flyback between two positions

x0 = scannerEndm1;
x1 = scannerEnd;
x2 = scannerStart;
x3 = scannerStartp1;
n = nFlybackPoints;


v1 = x1 - x0;
v2 = x3 - x2;
sampleI = (1:n)';
v1 = x1 - x0;
v2 = x3 - x2;
A = x2 - x1;
B = v2 - v1;
A2 = A - (v1 + B/2) * (n+1);
xx1 = x1 + (x2-x1)/(n+1)*sampleI;
xx2 = ((((v1 + v2) * (n+1)) - 2*(x2-x1))/(4*pi).*sin(2*pi/(n+1)*sampleI));
xx3 = ((v1 - v2)/(2*pi)*(n+1)*sin(pi/(n+1)*sampleI));
xx = xx1 + xx2 + xx3;
drivesignal = xx';
end