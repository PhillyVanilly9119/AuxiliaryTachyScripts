function [] = createRasterScanTable(savePath, aScans, nBuffer, fileName, fileSuffix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Function to create a by b raster scan remapping files
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xRange = 1:aScans;
yRange = 1:nBuffer;

% File name
fullFileName = strcat(fileName, '_', num2str(aScans), 'x', num2str(nBuffer), fileSuffix, '.bin');
fullPath = fullfile(savePath, fullFileName);

% Preallocate remapping table with desired dimensions
table = zeros(aScans*nBuffer, 4); % row 3 and 4 contain just zeros

% Create table: stride through buffers and create remapping positions
for idx = 1:nBuffer
    if idx == 1
        table(1:aScans, 1) = xRange;
        table(1:aScans, 2) = idx;
    else
        table( ((idx*aScans)+1) : (idx+1)*aScans, 1) = xRange;
        table( ((idx*aScans)+1) : (idx+1)*aScans, 2) = idx;
    end
end

% Safe table as binary file
f = fopen(fullPath, "wb");
fwrite(f, table, "uint32");
fclose(f);

end