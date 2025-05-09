%%
clc
clear all
close all

X = -10:0.01:10;
Y = -15*(sin(2*X).^2)-((X-2).^2) + 160;

[a, b] = size(X);
r4 = randperm(b);
XX = X(:, r4);
x = XX(1:10);

n = 1000;
selectionType = 2;  % 1 for Roulette Wheel, 2 for Tournament Selection
crossoverRate = 0.8;
mutationRate = 0.01;

FF = zeros(1, n);
FF_hat = zeros(1, n);

for i = 1:n
    fitness = -15*(sin(2*x).^2)-((x-2).^2) + 160;
    
    best = find(fitness == max(fitness));
    best1 = x(best);
    best2 = best1(1);

    FF(i) = max(fitness);

    % Selection mechanism
    if selectionType == 1  % Roulette Wheel Selection
        totalFitness = sum(fitness);
        prob = fitness / totalFitness;

        cumProb = cumsum(prob);
        newPop = zeros(1, 10);
        
        for j = 1:10
            randVal = rand;
            newPop(j) = x(find(cumProb >= randVal, 1));
        end
        x = newPop;
        
    elseif selectionType == 2  % Tournament Selection
        for k = 1:10
            r4 = randperm(10);
            candidates = x(r4(1:2));
            fitnessCandidates = fitness(r4(1:2));
            
            [~, winnerIdx] = max(fitnessCandidates);
            x(k) = candidates(winnerIdx);
        end
    end

    % Crossover (with crossover rate)
    for k = 1:2:7
        if rand < crossoverRate
            r5 = rand;
            r = -1:0.1:1;
            rr = randperm(21);
            sigma = r(:, rr);

            x(k) = x(k) + sigma(1) * (x(k) - x(k+1));
            x(k+1) = x(k+1) - sigma(1) * (x(k) - x(k+1));
        end
    end
    
    % Mutation (with mutation rate)
    if rand < mutationRate
        s = 0.01;
        rrr = randperm(10);
        x(1) = x(1) + s * rrr(1);
    end
    
    fitness_hat = -15*(sin(2*x).^2)-((x-2).^2) + 160;
    FF_hat(i) = max(fitness_hat);
    
    bad = find(fitness_hat == min(fitness_hat));
    x(bad) = best2;
end

% Output 
figure(1)
subplot(211), plot(X, Y)
xlabel('X')
ylabel('Fitness Value')
title('Fitness Function')

subplot(212), plot(FF_hat)
xlabel('Generation')
ylabel('Fitness Value')
title('Fitness Progression over Generations')

% Display the best solution found
bestSolution = x(find(fitness_hat == max(fitness_hat), 1));
disp(['最佳解 x = ', num2str(bestSolution)])
disp(['對應適應度 = ', num2str(max(FF_hat))])
