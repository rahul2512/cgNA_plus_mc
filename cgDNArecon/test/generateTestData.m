% ----------------
% Copyright 2015 Jaroslaw Glowacki
% jarek (dot) glowacki (at) gmail (dot) com
%
% This file is part of cgDNArecon.
%
% cgDNArecon is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% cgDNArecon is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with cgDNArecon.  If not, see <http://www.gnu.org/licenses/>.
% ----------------
%
% The function generates a mumber of random DNA sequences, computes their
% cgDNA minimum energy configurations using a given cgDNA parameter set.
% The parameter set is used to make the generated data "DNA-like".
% The data is saved as C++ code into a requested file.
%
% Arguments:
%     nbp         The length (number of base pairs) of the sequences to
%                 compute
%     numTests    The number of random sequences to generate
%     paramset    A cgDNA parameter set; this is of little importance as it
%                 is used only to make the constructed data "DNA-like"
%     filename    Name (path) for the file to save the generated C++ code in

function generateTestData(nbp, numTests, paramset, filename)
    n = 12 * nbp - 6;
    
    bases = 'ATGC';
      
    hatShapes = cell(1, numTests);
    r1 = cell(1, numTests);
    D1 = cell(1, numTests);
    r2 = cell(1, numTests);
    D2 = cell(1, numTests);
    
    for t = 1:numTests
        % Randomize a sequence
        inds = ceil(rand(1, nbp) * 4);
        seq = bases(inds);
        
        % Reconstruct its minimum energy configuration
        hatShapes{t} = constructSeqParms(seq, paramset);
        
        confData = frames(hatShapes{t});
        r1{t} = {confData.r};
        D1{t} = {confData.D};
        r2{t} = {confData.rc};
        D2{t} = {confData.Dc};
    end
    
    outFile = fopen(filename, 'w');
    
    fprintf(outFile, ['/** This is a generated file\n' ...
        ' * any changes will be overridden */\n\n'...
        'const int nbp = %d;\n', ...
        'const int numTests = %d;\n\n' ...
        'fpType shapesCorr[][%d] = {\n'], nbp, numTests, n);   
    
    % Print the shape vector
    for t = 1:numTests
        fprintf(outFile, '{');
        
        shape = hatShapes{t};
        for k = 1:n
            fprintf(outFile, '%17.14f', shape(k));
            if k ~= n
                fprintf(outFile, ', ');
            end
        end
        
        if t == numTests
            fprintf(outFile, '}\n');
        else
            fprintf(outFile, '},\n');
        end
    end
    
    fprintf(outFile, '};\n\n// The rest is not used for timings\n#ifndef TIMING_ONLY\n\n');
    
    %Print main strand points  
    printData(outFile, 'r1Corr', r1);
    printData(outFile, 'r2Corr', r2);
    printData(outFile, 'D1Corr', D1);
    printData(outFile, 'D2Corr', D2);

    fprintf(outFile, '#endif\n');
    
    fclose(outFile);
end

function printData(outFile, name, cel)
    c = cel{1};
    [c, r] = size(c{1});
    lowestDim = c * r;

    fprintf(outFile, 'fpType %s[][%d][%d] = {\n', name, length(cel{1}), lowestDim);
    for t = 1:length(cel)
        fprintf(outFile, '{\n');
        
        c1 = cel{t};
        for k = 1:length(c1)
            c2 = c1{k};
            fprintf(outFile, '\t{');
            for l = 1:lowestDim
                fprintf(outFile, '%16.13f', c2(l));
                if l ~= lowestDim
                    fprintf(outFile, ', ');
                end
            end
            if k == length(c1)
                fprintf(outFile, '}\n');
            else
                fprintf(outFile, '},\n');
            end
        end
        
        if t == length(cel)
            fprintf(outFile, '}\n');
        else
            fprintf(outFile, '},\n');
        end
    end
    
    fprintf(outFile, '};\n\n');
end
