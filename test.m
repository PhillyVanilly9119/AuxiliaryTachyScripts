

%% feedback X data
%  
filename = 'T:\home\lvuser\ScanPattern\feedback\19_11_v1_XScannerFeedback2.bin';
fid = fopen(filename,'r');
Xdata = fread(fid, 'single');
size_Xdata = size(Xdata);
% data = reshape(Xdata,[(size_Xdata(1)/2),2]);
fclose(fid);

Xdata = Xdata / max(Xdata);

figure;
plot(Xdata)


 %% feedback Y data
 
filename = 'T:\home\lvuser\ScanPattern\feedback\19_11_v1_YScannerFeedback2.bin';
fid = fopen(filename,'r');
Ydata = fread(fid, 'single');
size_Ydata = size(Ydata);
% data = reshape(data,[(size_data(1)/2),2]);
fclose(fid);

Ydata = Ydata / max(Ydata);

hold on
plot(Ydata)
%%
%% feedback X data
%  
filename = 'D:\1 4D OCT\Scanpattern\patternUpdate_08_26_2020\scannerfeedbackTest\1006remap2_Spiral_10_27_v4.bin';

fid = fopen(filename,'r');
Xvel = fread(fid, 'single');
size_Xvel = size(Xvel);
% data = reshape(Xdata,[(size_Xdata(1)/2),2]);
fclose(fid);

% Xvel = Xvel / max(Xvel);

figure;
plot(Xvel)


 %% feedback Y data
 
filename = 'T:\home\lvuser\ScanPattern\feedback\YScannerVelocity2.bin';
fid = fopen(filename,'r');
Yvel = fread(fid, 'single');
size_Yvel = size(Yvel);
% data = reshape(data,[(size_data(1)/2),2]);
fclose(fid);

% Yvel = Yvel / max(Yvel);

hold on
plot(Yvel)

%%
%save pos vect 1 spiral
Xsave=[Xges,Yges];
filename = 'D:\1 4D OCT\Scanpattern\patternUpdate_08_26_2020\scannerfeedbackTest\gem1006\XYpos.bin'
fid = fopen (filename,'w');

fwrite(fid, Xsave, 'single');
fclose(fid);



%% original data
filename = 'D:\1 4D OCT\Scanpattern\patternUpdate_08_26_2020\scannerfeedbackTest\2000_1_Spiral_11_16_v1.bin';
fid = fopen(filename,'r');
data_expected = fread(fid, 'uint32');
size_data_expected = size (data_expected);
data_expected = reshape (data_expected,[(size_data_expected(1)/4),4]);
fclose(fid);

%normierung
data_expectedX = data_expected(:,1) / max(abs(data_expected(:,1)));
data_expectedY = data_expected(:,2) / max(abs(data_expected(:,2)));
Xdata = Xdata / max(abs(Xdata));
Ydata = Ydata / max(abs(Ydata));

figure;

plot(data_expectedX,data_expectedY,'b')
hold on
plot(Xdata(4655:10478)+0.1,Ydata(4655:10478),'r')


%%
figure; 
plot(data_expectedX)
xlim([-5000 55000])

figure;
plot(data_expectedY)
xlim([-5000 55000])


%% calibraition

% X1= Xdata(11490:17682);
% 
% X2= Xdata(23880:30072);
% 
% X3= Xdata(30070:36262);
% 
% X4= Xdata(42460:48652);
% 
% X1= Xdata(2620:5716);
% 
% X2= Xdata(5717:8813);
% 
% X3= Xdata(8814:11910);


Xges= Xdata(50285:56434);

% Xges = (X1 + X2 + X3 + X4)/4;
% Xges = (X1 + X2 + X3 )/3;

% Y1= Ydata(11490:17682);
% 
% Y2= Ydata(23880:30072);
% 
% Y3= Ydata(30070:36262);
% 
% Y4= Ydata(42460:48652);

% Y1= Ydata(2620:5716);
% 
% Y2= Ydata(5717:8813);
% 
% Y3= Ydata(8814:11910);

Yges= Ydata(50285:56434);

% Yges = (Y1 + Y2 + Y3 )/3;

XYges = [Xges,Yges];

figure;
plot(Xges,Yges)

% XYges = XYges (37:3090,:);

% interpolieren
x= linspace (1,6150,61500);
y= linspace(1,2,2);
[X,Y]= meshgrid(x,y);
res = interp2(transpose(XYges),X,Y);
res=transpose(res);
figure; 
plot(res)
% figure;
% plot(res(60:49030,1),res(60:49030,2))
% 
% res = [res(267:end,:);res(1:266,:)];
figure; 
plot(res)
hold on 
 plot((scanPatternMM(:,1)./max(scanPatternMM(:,1))))

%%

filename = 'D:\1 4D OCT\Scanpattern\patternUpdate_08_26_2020\scannerfeedbackTest\2000\XYposinterp_2000_19_11_v1.bin'
fid = fopen (filename,'w');

fwrite(fid, res, 'single');
fclose(fid);


%%
filename = 'D:\1 4D OCT\Scanpattern\patternUpdate_08_26_2020\scannerfeedbackTest\gem1006\XYposinterp.bin';
fid = fopen(filename,'r');
Ydata = fread(fid, 'single');
size_Ydata = size(Ydata);
data = reshape(Ydata,[(size_Ydata(1)/2),2]);
fclose(fid);


%% distance between  points

for i = 1: 3096
    dis(i)=sqrt(Xdata(i+2619)^2 + Ydata(i+2619)^2);
end
    
