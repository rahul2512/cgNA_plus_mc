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
% The function saves a sequence and the provided stiffness matrix and
% shape vector in a text file formated for use with the cgDNAmc C++ Monte
% Carlo code. It allows for using any stiffness matrix and shape vector
% (not necessarily only cgDNA ones).
%
% This uses the cgDNA package function constructSeqParms.
%
% Arguments:
%     seq          The DNA sequence (a string of A, T, G or C)
%     K            The cgDNA stiffness matrix for the sequence seq
%     w_hat        The cgDNA shape vector for the sequence seq
%     do_factor    If true the provided matrix is Cholesky decomposed and
%                  the L Cholesky factor in symmetric band storage is
%                  saved; if false the matrix is saved "as is"; the main
%                  purpose of this function is to allow the use of other
%                  representations of the stiffness matrix
%                  (see the main documentation of cgDNAmc)
%     filename     Name (path) for the file to save the sequence seq and
%                  its stiffness matrix and shape vector in

function save_seq_data( seq, K, w_hat, do_factor, filename )   
    % Save first 10 meannigful digits
    number_format = '%16.9e ';

    if (do_factor)
        % Get the lower triangular matrix of the Cholesky decomposition
        L = chol(K, 'lower');
        % Put the L matrix in the band storage
        % (see LAPACK manual about band storage)
        K = fliplr(spdiags(L))';
    end
    
    [M, N] = size(K);
    
    fileId = fopen(filename, 'w');
    
    % The first line is the label of the data: the sequence
    fprintf(fileId, '%s\n', seq);
    
    % The second line gives the dimensions
    % Use column major order (as it will be used in C++)
    fprintf(fileId, '%d %d\n', N, M);

    % The matrix
    % Use column major order (as it will be used in C++)	
    for c = 1:N
        fprintf(fileId, number_format, K(:, c));
        fprintf(fileId, '\n');
    end
    fprintf(fileId, '\n');
    % The shape vector
    for c = 1:length(w_hat)
        fprintf(fileId, number_format, w_hat(c));
        fprintf(fileId, '\n');
    end    
    
    fclose(fileId);

end

