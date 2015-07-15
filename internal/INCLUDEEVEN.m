function  [LEVEN,LoEVEN] = INCLUDEEVEN(L,Lo)
%--------------------------------------------------------------------------
%INCLUDEEVEN.m
%
%This function includes even degree terms whenever odd degree appear in the
%expression. It increases the degree without increasing the depth. The
%resulting monomials will be considered for the integration by parts
%
%THIS IS THE CODE FOR ONE SPATIAL VARIABLE
%
%inputs - L:    matrix containing the list of monomials that appear on the
%               base monomials, and the indexes for the related monomials
%               with the spatial derivatives in matrix Lo. 
%               L in R^(nmb X (n+2))
%         Lo:   matrix containing the list of spatial derivatives on each
%               monomial.
%               Lo in R^(nm X dg)
%                           
%               nmb - number of base monomials, n - number of dependent
%               variables; nm - number of monomials; dg - degree of the
%               expression on the dependent variables. 
%
%outputs - idx: index for the list for each monomial giving the chain of
%               monomials with the spatial derivatives. This index must be
%               associated with the list of monomials in L. 
%               idx in R^(nmb X 2)
%
%          Lf:  the list of the spatial derivative to be taken with respect
%               to each position in the monomial
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


% 06/07/15 - GV




%first gets the degree
[~,n] = size(L);
degL = sum(L,2);
maxdegINIT = max(degL);%this is the maximal degree on the list
maxdeg = maxdegINIT;
oddelements = mod(degL,2);
indodd = find(oddelements);%gets the index of the odd degree terms

LEVEN = L;
LoEVEN = Lo;
for i = 1:length(indodd)
        for j=1:n
            LN = L(indodd(i),:);%extracts the indexed odd degree row of L
            LN(j)=LN(j)+1;%increase the degree of variable k
            LoN = Lo(indodd(i),:);%extracts the indexed odd degree row of Lo
            LoN = [LoN(1:sum(LN(1:j-1))) 0 LoN(sum(LN(1:j-1))+1:sum(LN(1:n))-1)];%inserts a zero corresponding to the increased degree. This does not increase the order
            if sum(LN)<=maxdegINIT
                LEVEN = [LEVEN; LN];
                LoEVEN = [LoEVEN; [LoN -1*ones(1,maxdeg-sum(LN))]];
            elseif sum(LN)>maxdegINIT
                maxdegprev = maxdeg;
                maxdeg = sum(LN);
                LEVEN = [LEVEN; LN];
                LoEVEN = [[LoEVEN  -1*ones(size(LoEVEN,1),sum(LN)-maxdegprev)]; LoN ];
            end
        end
end

[LUNIQUE] = unique([LEVEN LoEVEN] ,'rows','stable');
LEVEN = LUNIQUE(:,1:n);
LoEVEN = LUNIQUE(:,n+1:end);