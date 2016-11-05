function governingEqns = determineVdot()


charges = [0 0 1 1  -1 0 0 -1];
gasInput = [1 0 0 0 0 0 0 0 0];


%Stoichiometric Matrix is organized in the order of 1) Neutrals 2) Positive
%Ions, 3)Negative Ions, 4) Intermediates, 5) Electron, 6) Wall Density.
%User must input the number of neutral gases

stoichiometricMatrix = [-1	0	1	0	0	0	0	-1	0	;
-1	1	0	0	1	0	0	-1	0	;
-1	1	0	1	0	0	0	-1	0	;
-1	0	0	1	1	0	0	-1	0	;
-1	2	0	0	0	0	0	-1	0	;
-1	1	0	0	0	0	1	-1	0	;
0	-1	0	1	0	0	0	-1	0	;
0	2	-1	0	0	0	0	-1	0	;
1	-1	0	0	-1	0	0	0	0	;
-1	1	1	-1	0	0	0	0	0	;
-1	0	0	0	0	1	0	-1	0	;
1	0	0	0	0	-1	0	-1	0	;
0	0	1	0	0	-1	0	-1	0	;
0	1	0	0	1	-1	0	-1	0	;
0	2	0	0	0	-1	0	-1	0	;
1	1	0	0	-1	-1	0	0	0	;
0.5	-1	0	0	0	0	0	0	-1	;
1	0	0	0	0	-1	0	0	-1	;
1	0	-1	0	0	0	0	0	-1	;
0	1	0	-1	0	0	0	0	-1	;
1	1	-1	0	-1	0	0	0	0	;
0	3	-1	0	-1	0	0	0	0	;
0	2	0	-1	-1	0	0	0	0	];



%Specifices the number of each type of plasma species
numberOfChargedSpecies = nnz(charges)-1; %subtract one because we are not counting electrons
numberOfNeutrals = 2;
numberOfPositive = length(find(charges>0));
numberOfNegative = length(find(charges<0))-1;%not counting electrons
numberOfIntermediates = 2;
numberOfInputs = nnz(gasInput);
numberOfEqns = numberOfNeutrals + numberOfPositive+numberOfNegative + numberOfIntermediates + 1;

%Define starting indices in stoichiometric matrix for each type of species
neutralStart = 1;
neutralEnd = neutralStart + numberOfNeutrals-1;
positiveStart = neutralStart+numberOfNeutrals;
positiveEnd = positiveStart+numberOfPositive-1;
negativeStart = positiveStart+numberOfPositive;
intermediateStart = negativeStart + numberOfNegative;
numberOfGaseousSpecies = intermediateStart+numberOfIntermediates-1;

[numberOfReactions, numberOfSpecies] = size(stoichiometricMatrix);
kRates = ones(1,numberOfReactions);
speciesConcentrations = [2 1 2 1 2];
reducedMatrix = zeros(numberOfReactions,numberOfSpecies,numberOfGaseousSpecies);

%Flag recombination rxns
recReactionIndices = [];
wallReactionIndices = [];
for k = 1:numberOfReactions
    nonZeroElements = find(stoichiometricMatrix(k,1:numberOfGaseousSpecies)<0); %not including ne or wall
    f = 0;
    flagged = 0;
    while flagged == 0 && f < length(nonZeroElements)
        f = f + 1;
        if charges(nonZeroElements(f))>0
            for g = 1:length(nonZeroElements)
                if charges(nonZeroElements(g))<0; 
                    recReactionIndices = [recReactionIndices k];
                    flagged = 1;
                end
            end
        end
    end
end

for k = 1:numberOfReactions
        if stoichiometricMatrix(k,end)==-1
            wallReactionIndices = [wallReactionIndices k];
        end
end

numberOfWallReactions = length(wallReactionIndices);
numberOfRecombinationReactions = length(recReactionIndices);
for i = 1:length(recReactionIndices)
    recMatrix(i,:) = stoichiometricMatrix(recReactionIndices(i),:);
end
for i = 1:length(wallReactionIndices)
    wallMatrix(i,:) = stoichiometricMatrix(wallReactionIndices(i),:);
end

removalIndices = [recReactionIndices'; wallReactionIndices'];
stoichiometricMatrix(removalIndices,:) = [];

%add matrix back in
stoichiometricMatrix = [stoichiometricMatrix;recMatrix];
stoichiometricMatrix = [stoichiometricMatrix;wallMatrix];
startingRecIndex = numberOfReactions - numberOfWallReactions - numberOfRecombinationReactions+1;
endingRecIndex = numberOfReactions - numberOfWallReactions;
startingWallIndex = numberOfReactions - numberOfWallReactions +1;

%{'a2*a5*k2 - 2*a1^2*k1'}
for i = 1:numberOfReactions 
    for j = 1:numberOfGaseousSpecies 
            if stoichiometricMatrix(i,j)<0
                 for k = 1:numberOfSpecies 
                    if stoichiometricMatrix(i,k)<0
                         reducedMatrix(i,k,j) = abs(stoichiometricMatrix(i,k));
                    else
                         reducedMatrix(i,k,j) = 0;
                    end
                 end
             end
            if stoichiometricMatrix(i,j)>0
                for q = 1:numberOfSpecies 
                    if stoichiometricMatrix(i,q)<0
                        reducedMatrix(i,q,j) = abs(stoichiometricMatrix(i,q));
                    else
                        reducedMatrix(i,q,j) = 0;
                    end
                end
            end
            
        
    end
end



%kinetic balance equations
volumeSpeciesConcentrations = sym('n', [1 numberOfSpecies]);
peakSpeciesConcentrations = sym('np', [1 numberOfSpecies]);
gasInputs = sym('Q', [1 numberOfInputs]);
kRates = sym('k', [1 numberOfReactions]);
governingEqns = cell(numberOfEqns,1); %plus one for power eqn


%Write particle balance eqns

%All reactions except for recombination reactions (uses volume avgs for
%species concentrations)
for r = 1:numberOfGaseousSpecies %not doing reactions for walls or ne
    for rxn = 1:numberOfReactions-(numberOfRecombinationReactions + numberOfWallReactions)
        term = char(stoichiometricMatrix(rxn,r)*kRates(rxn));
        for species = 1:numberOfSpecies
            m1 = char(volumeSpeciesConcentrations(species)^reducedMatrix(rxn,species,r));
            if ~strcmp(m1,'1')
             term = strcat(term,[' * ' m1]);
            end
        end
            if ~isempty(governingEqns{r})&&~strcmp(term,'0')
                governingEqns{r} = strcat(governingEqns{r},[' + ' term]);
            elseif isempty(governingEqns{r})&&~strcmp(term,'0')
                governingEqns{r} = term;
            end
    end
    if charges(r) == 0
        governingEqns{r} = strcat(governingEqns{r},[' - Kpump*' char(volumeSpeciesConcentrations(r))]);
    end
end

%Add recombination reaction terms (uses peak species concentrations)
for r = 1:numberOfGaseousSpecies
    for rxn = startingRecIndex:endingRecIndex
        term = char(stoichiometricMatrix(rxn,r)*kRates(rxn));
        for species = 1:numberOfSpecies
            m1 = char(peakSpeciesConcentrations(species)^reducedMatrix(rxn,species,r));
            if ~strcmp(m1,'1')
             term = strcat(term,[' * ' m1]);
            end
        end
            if ~isempty(governingEqns{r})&&~strcmp(term,'0')
                governingEqns{r} = strcat(governingEqns{r},[' + ' term ' * Vrec/V']);
            elseif isempty(governingEqns{r})&&~strcmp(term,'0')
                governingEqns{r} = strcat(term,'*Vrec/V');
            end
    end

end

%Adds wall reaction terms (also uses peak species concentrations)
for r = 1:numberOfGaseousSpecies
    for rxn = startingWallIndex:numberOfReactions
        term = char(stoichiometricMatrix(rxn,r)*kRates(rxn));
        for species = 1:numberOfSpecies
            m1 = char(peakSpeciesConcentrations(species)^reducedMatrix(rxn,species,r));
            if ~strcmp(m1,'1')
             term = strcat(term,[' * ' m1]);
            end
        end
            if ~isempty(governingEqns{r})&&~strcmp(term,'0')
                governingEqns{r} = strcat(governingEqns{r},[' + ' term]);
            elseif isempty(governingEqns{r})&&~strcmp(term,'0')
                governingEqns{r} = term;
            end
    end

end

% Adds molecular flow rates (uses input gases)
for r = 1:numberOfGaseousSpecies
    if gasInput(r) ==1
      term = char(gasInputs(r));
            if ~isempty(governingEqns{r})&&~strcmp(term,'0')
                governingEqns{r} = strcat(governingEqns{r},[' + ' term '/V']);
            elseif isempty(governingEqns{r})&&~strcmp(term,'0')
                governingEqns{r} = strcat(term, '/V');
            end
    end

end


%Write power balance eqn

%flags ionization reactions
collReactionIndices = [];
for k = 1:numberOfReactions
    for i = neutralStart:neutralEnd
            if stoichiometricMatrix(k,i)<0
                if stoichiometricMatrix(k,i+numberOfNeutrals)>0 && nnz(stoichiometricMatrix(k,i+1:numberOfGaseousSpecies))==1 %looks at corresponding positive ion to neutral
                    collReactionIndices = [collReactionIndices k];
                end
            end
    end 
end

%flag wall reactions with positive ions
positiveWallReactionIndices = [];
for k = numberOfReactions-numberOfWallReactions+1:numberOfReactions
    for j = positiveStart:positiveEnd
        if stoichiometricMatrix(k,j) < 0
            positiveWallReactionIndices = [positiveWallReactionIndices k];
        end
    end
end

%Write collisional energy loss terms
numberOfPositiveWallReactions = length(positiveWallReactionIndices);
numberOfCollReactions = length(collReactionIndices);
collEnergyLoss = sym('Ec_', [1 numberOfCollReactions]);
wallEnergyLoss = sym('Ew_', [ 1 numberOfPositiveWallReactions]);
powerEqn = {'pabs'};
 for i = 1:numberOfCollReactions %these need to correspond to neutrals (neutral wall reactions)
            coll_term = char(collEnergyLoss(i)*kRates(collReactionIndices(i)));
            for species = 1:numberOfSpecies
                m1 = char(volumeSpeciesConcentrations(species)^reducedMatrix(collReactionIndices(i),species,i));
                if ~strcmp(m1,'1')
                    coll_term = strcat(coll_term,[' * ' m1])
                end
            end
 powerEqn = strcat(powerEqn,[' - ' ,coll_term])
 end
 
 %Write wall terms
 wallIndex = 0;
  for i = positiveStart:positiveEnd
            wallIndex = wallIndex+1;
            wall_term = char(wallEnergyLoss(wallIndex)*kRates(positiveWallReactionIndices(wallIndex)));
            for species = 1:numberOfSpecies
                m1 = char(peakSpeciesConcentrations(species)^reducedMatrix(positiveWallReactionIndices(wallIndex),species,i));
                if ~strcmp(m1,'1')
                    wall_term = strcat(wall_term,[' * ' m1])
                end
            end
 powerEqn = strcat(powerEqn,[' - ' ,wall_term])
 end

%Add power eqn to group of eqns
governingEqns{end} = powerEqn;
























% 
% myValues = ['k2 '; 'k3 '; 'k4 '; 'k5 '; 'k6 '; 'k7 '; 'k8 '; 'k9 '; 'k10'; 'k11'; 'k12'; 'k13'; 'k14'; 'k15'; 'k16'; 'k17'; 'k18'; 'k19'; 'k20'; 'k21'; 'k22'; 'k23';'k1 ']
% cellMyValues = cellstr(myValues);
% paperValues = ['Katt '; 'Kiz4 '; 'Kiz3 '; 'Kdiss'; 'Kdiss'; 'Kiz2 '; 'Kei  '; 'Kdet '; 'Kch  '; 'Kex  '; 'Kdeex'; 'Kizm '; 'Kattm'; 'Kdism'; 'Krec4'; 'vO   '; 'vO2* '; 'vO2+ '; 'vO+  '; 'Krec '; 'Krec2'; 'Krec3';'Kiz1 '];
% cellPaperValues = cellstr(paperValues);
% mySpecies = ['n1 '; 'n2 '; 'n3 ';'n4 '; 'n5 '; 'n6 '; 'n7 '; 'n8 ';'np3';'np4';'np5'];
% cellMySpecies = cellstr(mySpecies);
% paperSpecies = ['nO2  '; 'nO   '; 'nO2+ ';'nO+  '; 'nO-  '; 'nO2* '; 'nO*  '; 'ne   ';'nO2+p';'nO+p ';'nO-p '];
% cellPaperSpecies = cellstr(paperSpecies);
% for j = 1:numberOfSpecies-1
%     for i=1:numberOfReactions
%         eqntest{j} =strrep(eqntest{j},cellMyValues(i),cellPaperValues(i));
%     end
% end

% for j = 1:numberOfSpecies-1
%     for i=1:length(mySpecies)
%         eqntest{j} =strrep(eqntest{j},cellMySpecies(i),cellPaperSpecies(i));
%     end
% end


%mass balance eqn 

% 
% %quasineutrality condition 
% %[Cl2 Cl Cl2+ Cl+ Cl- ne]
% 
% quasineutralityCondition = [];
% for i=1:numberOfSpecies
%     q = char(speciesConcentrations(i)*charges(i))
%     if ~isempty(quasineutralityCondition)&&~strcmp(q,'0')
%         quasineutralityCondition = strcat(quasineutralityCondition,[' + ' q]);
%     elseif isempty(quasineutralityCondition)&&~strcmp(q,'0')
%         quasineutralityCondition = q;
%     end
% end

%bohm velocity eqn




%power balance equations
