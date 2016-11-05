function F = Prior(current,index)
global priorLB
global priorUB
global priorRecord

priorTime = tic;
prior = 0; %for when you are not doing cwise MH

 prior = log(unifpdf(current(index),priorLB(index),priorUB(index))) + prior;

F = prior;
priorRecord = [priorRecord prior];
priorTimeElapsed = toc(priorTime);
assignin('base', 'priorTimeElapsed', priorTimeElapsed);
end