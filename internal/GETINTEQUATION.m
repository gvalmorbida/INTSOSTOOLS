function  [coeffs,Z,DvarDeg,DvarOrd] = GETINTEQUATION(symexpr,DepVar)
%--------------------------------------------------------------------------
%GETINTEQUATION.m
%
%THIS IS THE CODE FOR ONE SPATIAL VARIABLE
%
%inputs -           symexpr:        the kernel of the integral
%                   DepVar:         the matrix containing the dependent
%                                   variables
%                   
%
%
%outputs -          coeffs:         a vector with the coefficients of the
%                                   expression
%                   Z:              the set of monomials in the dependent
%                                   variable
%                   DvarDeg:        the matrix with the degree of the
%                                   monomials
%                   DvarOrd:        the matrix with the order in the
%                                   monomials
%
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

%extracts the monomials and the coefficients from symexpr
if isa(symexpr,'polynomial')
    cvartable =  char(DepVar.varname);%creates a string with the dependent variables

    % Collect symexpr = g(c)*h(x) where x are the vars in vartable,
    % c are the decision variables, and h(x) is a vector of monomials.
    % in the polynomial variables.  Use this to find unique
    % monomials in the polynomial variables.
    [g0,g,h] = collect(symexpr(:),setdiff(symexpr.varname,cvartable));
    g = g(:);
    if ~isequal(g0,0)
        g = [g0;g];
        h = [1; h];
    end
    [nmon,~] = size(h.degmat);
    nvar = length(cvartable);
    
    % Reorder the monomials matrix with variables in the order
    % listed in cvartable and sorted monomials
    Z = zeros(nmon,nvar);
    [~,idx]=ismember(h.varname,cvartable);%extracts the indices to replace in the column of the degrees in Z corresponding to the variable in cvartable
    Z(:,idx) = full(h.degmat);
    %Z = sortNoRepeat([],Z);
    Z = sparse(Z);
    
    [nmon,nvar] = size(Z);
    coeffs = polynomial(sparse(nmon,1));
    
    for k = 1:length(h);
        coeffs(k) = g(k);
    end
    
    
else
    Z = [];
    charvartable = converttochar([vec(DepVar)']);
    coefmon = feval(symengine,'poly2list',symexpr,charvartable);
    coefmon = coefmon.';
    coeffs = sym('coeffs',[size(coefmon,1),1]);
    for i = 1:size(coefmon,1)
        dummyvar = reshape(coefmon(i,:),2,1);
        Z(i,:) = double(dummyvar(2));
        coeffs(i,:) =  dummyvar(1) ;
    end
    Z = sparse(Z);
    
end

degmax = max(sum(Z,2));
DvarDeg = sparse(size(Z,1),n);
DvarOrd = (-1)*ones(size(Z,1),degmax);%as we put -1 in DvarOrd columns  when the degree is not equal to the maximal
for i = 1:size(Z,1)
    M = reshape(Z(i,:),n,ord);
    DvarDeg(i,:) = sum(M,2)';
    [cOrd,rOrd] = find(M');  %extracts the non-zero elements of the matrix. 
                             %needs the transpose find() orders columnwise
    vecj = [];
    for j = 1:length(cOrd)
        vecj = [vecj (cOrd(j)-1)*ones(1,M(rOrd(j),cOrd(j)))];
    end
                             
    DvarOrd(i,1:length(vecj)) = vecj;
end

%returns the ordered set of monomials
[DvarDeg,indrow] = sortrows(DvarDeg);
DvarOrd = DvarOrd(indrow,:);
totdegDvarDeg = sum(DvarDeg,2);
[~,inddeg] = sort(totdegDvarDeg);
DvarDeg = DvarDeg(inddeg,:);
DvarOrd = DvarOrd(inddeg,:);
