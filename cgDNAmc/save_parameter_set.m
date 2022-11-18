% ----------------
% Copyright 2016 Jaroslaw Glowacki
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
% The function saves a cgDNA parameter set given as argument in a text file
% formated for use with the cgDNAmc C++ Monte Carlo code.
%
% Arguments:
%     paramset    The cgDNA parameter set structure as defined by the cgDNA
%                 package
%     filename    Name (path) for the file to save the parameter set in

function save_parameter_set( paramset, filename )
    base = paramset.base;
    dimer = paramset.dimer;
    base_len = length(base);
    dimer_len = length(dimer);
    
    % Arrange everything in alphabetic order so that it is easier to compare
    % output file using simply diff
    % This is to unify the arrangement of data that is different in
    % different vesions of parameter sets
    [~, inds_d] = sort({dimer.S});
    [~, inds_b] = sort({base.S});
    
    number_format = '%16.12f ';

    file_id = fopen(filename, 'w');
    
    % Write dimer data
    for k = 1:dimer_len
        d = dimer(inds_d(k));
        if d.S ~= 'av'
            fprintf(file_id, '%s\n', d.S);
            for r = 1:18
                fprintf(file_id, number_format, d.b18(r, :));
                fprintf(file_id, '\n');
            end
            fprintf(file_id, '\n');
            fprintf(file_id, number_format, d.c18);
            fprintf(file_id, '\n\n');
        end
    end
    
    % Write base data
    for k = 1:base_len
        b = base(inds_b(k));
        if b.S ~= 'av'
            fprintf(file_id, '%s\n', b.S);
            for r = 1:6
                fprintf(file_id, number_format, b.b(r, :));
                fprintf(file_id, '\n');
            end
            fprintf(file_id, '\n');
            fprintf(file_id, number_format, b.c);
            fprintf(file_id, '\n\n');
        end
    end
    
    fclose(file_id);
end

