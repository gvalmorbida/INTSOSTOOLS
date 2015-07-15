function  [sosprog,gx] = POLYGRAM(sosprog,monlist,degG,Indvar,Depvar)
%--------------------------------------------------------------------------
%POLYGRAM.m
%
%This function builds the polynomial Gram matrix given the list of
%monomials to be in the quadratic-like form in the dependent variables. 
%
%
%
%THIS IS THE CODE FOR ONE SPATIAL VARIABLE
%
%inputs - monlist:  The list of monomials on the dependent variable and its
%                   derivatives to generate the quadratic-like expression. 
%         degG:     The degree of the polynomial gram matrix   
%                           
%         Indvar:   The set of independent variables
%
%outputs - Gx: the polynomial matrix
%
%
%example:  
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


[n,ord] = size(Depvar);
Z = [];

vecDvar = vec(Depvar)';%extracts the list of dependent variables
charvartable = converttochar([vec(Depvar)']);%obtains the string of characters to be passed to the symbolic engine
coefmon = feval(symengine,'poly2list',sum(monlist),charvartable);%calls the symbolic engine to obtain the list of monomials on the expression
coefmon = coefmon.';%this is the list of coefficients and exponents
for i = 1:size(coefmon,1)
    dummyvar = reshape(coefmon(i,:),2,1);
    Z(i,:) = double(dummyvar(2));
end

[A,ZZ] = getconstraint(Z);%ZZ contains the monomials of Z kron Z, and A is the matrix combining the elements of vec(Q) as A*vec(Q)=F and we have Z'QZ = F*ZZ
[~,dimG] = size(A);%the dimension of G is the square root of the number of elements of A
IND = find(sum(A,2)>2);%gets the indexes of rows of matrix A that have more than 2 entries. 
IND2 = find(sum(A(IND,:),1)>0);%extracts the indices of the rows having positive sum
[sosprog,g] = sospolymatrixvar(sosprog,monomials(Indvar,0:degG),[length(IND2) 1]);%creates the polynomials to be introduced in the Gram matrix
gz = sym(zeros(dimG,1));%creates a full (vectorized) matrix of zeros
gz(IND2) = g;%inserts the nonzero entries in the positions defined by IND2
sosprog = soseq(sosprog,A(IND,:)*gz);%sets the gram matrix conditions
protoG = A(IND,:)*diag(gz);%introduces the polynomials on each element of A
G = mat(sum(protoG,1),sqrt(dimG));%creates the squared matrix
monlistZ = mysympower(Depvar,Z);%it is a permutation of monlinst, required to be compatible with G
gx =  monlistZ'*(G+G.')*monlistZ;

