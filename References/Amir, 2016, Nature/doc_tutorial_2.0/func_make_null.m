function A_null=func_make_null(A,N_null)
% THE FUNCTION ASSUMES THAT UNASSIGNED SPECIES (IN CASE OF GENUS) ARE
% ALREADY REMOVED - SO THE TOTAL ABUNDANCE OF EACH SAMPLE <=1

if nargin==1
    N_null=1;
end

[NumSpecies, NumSamples]=size(A);

NumSamples_null=N_null*NumSamples;

NotAssigned=1-sum(A); % For the case of genus level with unassigned reads

A_null=zeros(NumSpecies, NumSamples_null);
NotAssigned_null=zeros(1,NumSamples_null);
for j=1:NumSamples_null
    j_real=mod(j-1,NumSamples)+1;
    for i=1:NumSpecies
        if A(i,j_real)>0
            Inonzeros=find(A(i,:)>0);
            RAND=randi(length(Inonzeros));
            A_null(i,j)=A(i,Inonzeros(RAND));
        else
            A_null(i,j)=0;
        end
    end
    NotAssigned_null(j)=NotAssigned(randi(NumSamples));
end

% normalizes the shuffled data together with the NotAssigned
A_null=A_null./repmat(sum(A_null)+NotAssigned_null,NumSpecies,1);






