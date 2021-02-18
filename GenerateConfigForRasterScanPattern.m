% function [scanTable] = GenerateConfigForRasterScanPattern(lenA, lenFlyback, isBidirectional)

lenA = 128;
lenFlyback = 20;
isBidirectional = 1;

scanTable = zeros( (lenA + lenFlyback) * lenA, 4);

if isBidirectional == 1
    for a = 1:lenA
        if mod(a,2) == 0 % EVEN b-Scans position
            scanTable(a*(lenA + lenFlyback): a*(lenA + lenFlyback) + lenA, 1) = bScanPositions(lenA);
            scanTable(a*(lenA + lenFlyback): a*(lenA + lenFlyback) + lenA, 2) = a;
        else % UNEVEN b-Scans position
            scanTable(a*(lenA + lenFlyback): a*(lenA + lenFlyback) + lenA, 1) = inversebScanPositions(lenA);
            scanTable(a*(lenA + lenFlyback): a*(lenA + lenFlyback) + lenA, 2) = a;
        end
    end
% else
%     for a = 1:lenA
%         scanTable();
%     end
end

scanTable = uint32(scanTable);

%% Local functions for b-Scan position generation
    function [bScanLine] = bScanPositions(aLen)
        bScanLine = [1:aLen, 4];
    end

    function [inversebScanLine] = inversebScanPositions(aLen)
        inversebScanLine(aLen:1) = [aLen:1, 1];
    end

% end