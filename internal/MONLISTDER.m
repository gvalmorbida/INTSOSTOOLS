function  [monlist] = MONLISTDER(DepVar,L,Lo)
%--------------------------------------------------------------------------
%MONLISTDER.m
%
%This function takes a list of monomials on a variable and its (spatial)
%derivatives and returns a vector containing the monomials as symbolic
%variables. These are taken from the list of monomials DepVar
%
%THIS IS THE CODE FOR ONE SPATIAL VARIABLE
%
%inputs - DepVar: DepVar contains all the dependent variables and its derivatives 
%         L:      matrix containing the list of monomials that appear on the
%                 base monomials, and the indexes for the related monomials
%                 with the spatial derivatives in matrix Lo. 
%                 L in R^(nmb X (n+2))
%         Lf:     matrix containing the list of spatial derivatives on each
%                 monomial.
%                 Lo in R^(nm X dg)
%
%output - monlist:The list containing the monomials written in L,Lf and idx
%
%example:
%           DepVar = [u1 u1x u1xx; u2 u2x u2xx];
%           L = [1 1; 0 2]; u1u2, u2^2
%           Lf = [0 2;  1 1];the derivative for u1u2xx, u2x^2
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








[n,ord] = size(DepVar);
DepVar = [ones(1,ord+1); [ones(n,1) DepVar]];%introduces lines of 1's in the first row and in the first column

%CAN MAKE SPECIAL FUNCTIONS TO CONVERT FROM ONE FORMAT TO ANOTHER
%rewrites L to be in the index format

%assumes L and Lf have the same number of rows

%convert to the index format. 
Lmod = LMATN2DEG(L);

Lmod = Lmod+1;%needs to update as indexes can not be zero.
Lo = Lo+2;%needs to update as indexes can not be zero.


matInd = sub2ind(size(DepVar),Lmod,Lo(:,1:size(Lmod,2)));%gets the indexes on the matrices to perform the product. 
%it is a matrix of indexes as we are required to have a list of monomials

%the argument (:,1:size(Lmod,2)) is to account for Lf matrices that
%correspond to a reduced set of degree matrix, i.e. Lf matrices containing
%columns with only -1. 


if isa(DepVar,'polynomial')%if it is a pvar
    monlist = polynomial(sparse(size(matInd,1),1));
    for i = 1:size(matInd,1)%needs to build the vector line by line
        monlist(i) = prod(DepVar(matInd(i,:)),2);
    end
else%if it is a sym
    monlist = prod(DepVar(matInd),2);%extract the list of monomials by muliplying all columns of DepVar(matInd)
end