function  [Lip,Loip,LipF,LoipF] = IPMONOMIALS(L,Lo)
%--------------------------------------------------------------------------
%IPMONOMIALS.m
%
%This function takes a list of monomials on a variable and its (spatial)
%derivatives and returns the list of monomials on the tree starting from
%the base monomial (without derivatives, on the same variables) which lead
%to them when differentiating to the depth of each monomial. 
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
%example:   monomials of degree 2: 
%           u1u2_xx,u1_xu2_x u2_x^2
%           L = [1 1;1 1; 0 2]; u1u2, u2^2
%           Lo = [0 2; 1 1; 1 1];the derivatives for u1u2, u1u2, u2u2
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





[nmb,n] = size(L);%nmb: number of base monomials, n: number of dependent variables 
[nm,dg] = size(Lo);%nm: number of monomials, dg: degree of the expression

Loip = [];%matrix to receive the set of monomials on the tree of monomials generating Lo
Lip = [];
for i = 1:nmb%for all monomials on the list of base monomials
    monl = [];
    dgk = sum(L(i,:));

    if sum(Lo(i,1:dgk),2)~=0%if the order of the monomial is not zero then introduce the elements on the chain of monomials for the integration by parts
        
        for j = 1:n%for each variable
            nrep = max(size(monl,1),1);
            indLo = sum(L(i,1:j));%this is the column index to be checked in Lo matrix for the maximal order on variable l
            if L(i,j)==0%if the degree of the variable is zero
                Lrep = [];
            elseif Lo(i,indLo)==0;%if the element is of order zero
                Lrep = sparse(1,L(i,j));
            else
                Lrep = sparse(mycombnk(Lo(i,indLo)-1,L(i,j)));%takes a number of elements given by the degree and obtains a combination of repeated values up to order Lo(i,indLo)-1
                %Lrep = sparse(combnk(repmat(0:Lo(i,indLo)-1,1,L(i,l)),L(i,l)));%takes a number of elements given by the degree and obtains a combination of repeated values up to order Lo(i,indLo)-1
            end
            Lodegl = sort(Lrep,2);%gives the sorted ordered matrix for variable #l, since Lo is ordered takes does not increase the depth according to what is in Lo
            %Lodegl = unique(sort(Lrep,2),'rows');%gives the sorted ordered matrix for variable #l, since Lo is ordered takes does not increase the depth according to what is in Lo
            R = repmat(Lodegl,1,nrep)';%these lines avoid the kronecker product
            Rf = reshape(R,size(Lodegl,2),nrep*size(Lodegl,1))';%repeats each row the number of times to match the exiting monl, providing all possible combinations.
            if L(i,j)~=0 %add columns to the corresponding Lo matrix only if the degree is not 0
                monl = [repmat(monl,size(Lodegl,1),1) Rf];%adds the new set of columns
            end
        end
        Loip = [Loip; [monl -1*ones(size(monl,1),(dg-dgk))]];%includes the monomials in the tree of monomials
        Lip = [Lip; repmat(L(i,:),size(monl,1),1)];%replicates the line to contain all the monomials
        
    end
end
[LUNIQUE] = unique([Lip Loip] ,'rows','stable');
Lip = LUNIQUE(:,1:n);
Loip = LUNIQUE(:,n+1:end);


nIP = size(Lip,1);

LipF = [];
LoipF = [];%matrix to receive the set of monomials on with branches from Lip
for i = 1:nIP

    Lip_i = Lip(i,:);
    Loip_i = Loip(i,:);%matrix to receive the set of monomials on with branches from Lip
    nLo = 1;
    for j = find(Lip(i,:))
        ind1 = sum(Lip(i,1:j-1))+1;
        ind2 = sum(Lip(i,1:j));
        nel = ind2-ind1+1;
        M = repmat(Loip(i,ind1:ind2),[nel,1]) +eye(nel);
        M = sort(M,2);
        MUNIQUE = unique(M,'rows','stable');
        nrows=size(MUNIQUE,1);
        Lip_i = [Lip_i; repmat(Lip(i,:),[nrows,1])];
        Loip_i = [Loip_i; repmat(Loip(i,:),[nrows,1])];
        Loip_i(nLo+1:end,ind1:ind2) = MUNIQUE;
        nLo = nLo +nrows;
        %ind = sum(Lip(i,1:j),2);
        %Loip_ij = Loip(i,:);Loip_ij(ind)=Loip_ij(ind)+1;
        %Loip_i =[Loip_i; Loip_ij];%matrix to receive the set of monomials on with branches from Lip
    end
    LipF = [LipF; Lip_i];
    LoipF = [LoipF; Loip_i];
end