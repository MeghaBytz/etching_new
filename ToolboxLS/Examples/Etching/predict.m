clear all
close all

%this method uses Bayesian inference to infer ion sputtering yield,
%reaction probabilty of neutrals and ion stimulated desorption (all
%hard to measure constants in the surface kinetics model for etching). 

%declare global variables
global proposalLB
global proposalUB
global expData
global noUnknowns
global priorLB
global priorUB
global expRows
global expColumns
global noExperiments
global proposalMean
global proposalSD
global numberOfIons
global numberOfRadicals
global expConditions
global y
global noTrainingData
global noTestData
global noExperimentalMetrics
%global algernonEtchRatesWithNoise
global levelSetUnknowns
    global numberOfMaterials 

noExperimentalMetrics = 2;
numberOfIons = 2;
numberOfRadicals = 2;
levelSetUnknowns = 2;
noUnknowns = numberOfIons*2+numberOfRadicals*2+levelSetUnknowns+3;
proposalLB = [0 0 1 1 0 0 0 0 .8 0 1 1 .001];
proposalUB = [2 2 10 10 1 1 1 1 5 1 100 50 .1];
proposalMean = [1 1 5 5 .5 .5 .5 .5 3 .5 60 7 .05];
proposalSD =  [.2 .2 1 1 .1 .1 .1 .1 1 .2 5 3 .01];
priorLB = [0 0 1 1 0 0 0 0 .8 0 40 1 .001];
priorUB = [2 2 8 8 1 1 1 1 5 1 80 10 .1];
numberOfMaterials = 1;


%Bayesian model method
load allSynExpConditions
load algernonTrenchDim
%Make set of training/test data %no added noise
noTrainingData = 5;
noTestData = 10;
k = noTrainingData + noTestData;
y = randsample(length(algernonTrenchDim),k);
%data = algernonEtchRatesWithNoise(y(1:noTrainingData));
expConditions = allSynExpConditions(y(1:noTrainingData),:);
expData = algernonTrenchDim(y(1:noTrainingData),:);
%realParameters = [.05 .05 .05 .05 0.25]
%Make initial guess for unknown parameters
current = ones(noUnknowns,1);

%Peform MH iterations

N = 50;
theta = zeros(N*noUnknowns,noUnknowns);
acc = zeros(1,noUnknowns);
PosteriorCurrent = Posterior(current,1);
% for i = 1:burnin    % First make the burn-in stage
%     for j=1:noUnknowns
%       [alpha,t, a,prob, PosteriorCatch] = MetropolisHastings(theta,current,PosteriorCurrent,j);
%     end
% end
tic
for cycle = 1:N  % Cycle to the number of samples
    ind = 1;
     for m=1:noUnknowns %assumingknown experimental noise % Cycle to make the thinning
            [alpha,t, a,prob, PosteriorCatch] = MetropolisHastings(theta,current,PosteriorCurrent,m);
            theta((cycle-1)*noUnknowns+m,:) = t;        % Samples accepted
            AlphaSet(cycle,m) = alpha;
            current = t;
            PosteriorCurrent = PosteriorCatch;
            acc(m) = acc(m) + a;  % Accepted ?
     end
     cycle
end
 
elaspedTime = toc;
accrate = acc/N;     % Acceptance rate,. 

%calculate final etch rate from Bayes using mean values
for i =1:noUnknowns
    inferred(i) = mean(theta(:,i));
end
for k=1:length(expData)
        [data, g, data0, inferredTrenchDim(k,1), inferredTrenchDim(k,2)] = etchingWithMask_new('high', 'contour', 'varyAlpha',inferred,k);
end

%Compare real versus prediction
exp = linspace(1,length(expData),length(expData));
figure (1)
scatter(exp,inferredTrenchDim(1:5,1)');
hold on
scatter(exp,expData(:,1),'r');
title('Real versus predicted Width (based on mean parameter values');
xlabel('Exp No');
ylabel('Width(um)');

%Compare real versus prediction
exp = linspace(1,length(expData),length(expData));
figure (2)
scatter(exp,inferredTrenchDim(1:5,2));
hold on
scatter(exp,expData(:,2),'r');
title('Real versus predicted Height (based on mean parameter values');
xlabel('Exp No');
ylabel('Height (um)')
figure;
    for i =1:noUnknowns
        subplot(4,4,i);
        outputTitle = sprintf('Unknown %d',i);
        hist(theta(:,i)); 
 %       f.Normalization = 'probability';
        %counts = f.Values;
%         hold on
%         line([real(i) real(i)],[0 max(counts)], 'Color', 'r')
%         hold on
%         xmin = min(theta(:,i));
%         xmax = max(theta(:,i));
%         x = xmin:1:xmax;
%         hold on
%         line([priorLB(i) priorLB(i)],[0 max(counts)], 'Color', 'g')
%         hold on
%         line([priorUB(i) priorUB(i)],[0 max(counts)], 'Color', 'g')
        title(outputTitle);
        xlabel('Value');
        ylabel('Frequency');
    end



figure; 
for k =1:noUnknowns
    outputTitle = sprintf('Unknown %d',k);
    subplot(4,4,k);
    plot(theta(:,k));
    hold on
    line([0 N],[real(k) real(k)], 'Color', 'r')
    hold on
    line([0 N],[proposalLB(i) proposalLB(i)], 'Color', 'g')
    line([0 N],[proposalUB(i) proposalUB(i)], 'Color', 'g')
    title(outputTitle);
    xlabel('Cycle #');
    ylabel('Value');
end

%genPlots(theta)
