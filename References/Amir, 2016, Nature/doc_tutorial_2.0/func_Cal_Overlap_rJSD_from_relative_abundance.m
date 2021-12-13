function [Overlap,RootJSD]=func_Cal_Overlap_rJSD_from_relative_abundance(A)
% THIS FUNCTION ASSUMES A IS ALREADY NORMALIZED

[~, NumSamples]=size(A);

%% Prepare analysis
% preallocate memory 
Overlap=nan(NumSamples,NumSamples);
RootJSD=nan(NumSamples,NumSamples);

% Define Kullback-Leibler Divergence
KLD=@(x,y) sum(x.*log(x./y));

% Define Jason-Shanon Divergence
rJSD=@(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2));

%% Analize all sample pairs

% PairsCounter=0;
for i=1:NumSamples-1
    for j=i+1:NumSamples
%         PairsCounter=PairsCounter+1;
        
        % Find the shared species
        SharedSpecies=find(A(:,i)>0 & A(:,j)>0);
        
        % Calculate the overlap
        Overlap(i,j)=0.5*(sum(A(SharedSpecies,i))...
            +sum(A(SharedSpecies,j)));
        
        % Renormalize the shared species
        RenormalizedA=A(SharedSpecies,i)/sum(A(SharedSpecies,i));
        RenormalizedB=A(SharedSpecies,j)/sum(A(SharedSpecies,j));
        
        % Calculate the rJSD dissimilarity
        RootJSD(i,j)=rJSD(RenormalizedA,RenormalizedB);
        
    end % for j
end % for i



