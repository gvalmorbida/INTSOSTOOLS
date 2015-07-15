function  [Lmod] = LMATN2DEG(L)
%--------------------------------------------------------------------------
%LMATN2DEG.m
%
%This function converts the matrix to 
%
%THIS IS THE CODE FOR ONE SPATIAL VARIABLE
%
%inputs - L:matrix containing the monomials with its columns being the
%           degree on each dependent variable
%
%output - Lmod:matrix containing the monomials with its colums being the
%              dependent variable to be taken the power
%
%example:
%           L = [1 1; 0 3]; u1u2, u2^3
%           Lmod = [1 2 0; 2 2 2];   
%--------------------------------------------------------------------------
%
% This file is part of INTSOSTOOLS - INTSOSTOOLS ver 1.00 a plug-in to SOSOTOOLS for integral inequalities.
%
% Copyright (C)2015  G. Valmorbida, A. Papachristodoulou

% Department of Engineering Science, University of Oxford, Oxford, U.K.
%
% Send bug reports and feedback to: giorgio.valmorbida@eng.ox.ac.uk
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------


% 06/07/15 - GV




[nm,n] = size(L);
degmax=max(sum(L,2));
Lmod = sparse(size(L,1),degmax);%allocates space for Lmod
                                %entries are zero if the degree of the monomial is not degmax
for i = 1:size(L,1)
    ind=0;
    for j = 1:n
        Lmod(i,ind+1:ind+L(i,j)) = j*ones(1,L(i,j));%fills the matrix Lmod in with the indexes of the elements on the expression
        ind = ind+L(i,j);
    end
end