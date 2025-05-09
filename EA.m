%%
clc; clear; close all;

X = -10:0.02:10;
Y = -15 * (sin(2 * X).^2) - (X - 2).^2 + 160;
figure(1); subplot(2,1,1); plot(X, Y);
xlabel('X'); ylabel('Fitness'); title('適應度函數');

% === 參數設定 ===
pop_size = 10;
n_gen = 1000;
crossover_rate = 0.8;
mutation_rate = 0.01;
selection_type = 2;  % 1: Tournament, 2: Roulette Wheel

% === 初始族群 ===
x = X(randperm(length(X), pop_size));
Fitness_f = zeros(1, n_gen);

for gen = 1:n_gen
    % === 計算適應度 ===
    fitness = -15 * (sin(2 * x).^2) - (x - 2).^2 + 160;
    [max_fit, idx_best] = max(fitness);
    best_individual = x(idx_best);
    Fitness_f(gen) = max_fit;

    % === 選擇 ===
    new_x = zeros(1, pop_size);
    switch selection_type
        case 1  % Tournament Selection
            for i = 1:pop_size
                c = randperm(pop_size, 2);
                if fitness(c(1)) > fitness(c(2))
                    new_x(i) = x(c(1));
                else
                    new_x(i) = x(c(2));
                end
            end
        case 2  % Roulette Wheel Selection
            prob = fitness / sum(fitness);
            cum_prob = cumsum(prob);
            for i = 1:pop_size
                r = rand;
                idx = find(cum_prob >= r, 1);
                new_x(i) = x(idx);
            end
    end
    x = new_x;

    % === Crossover（實數交配，交配率控制）===
    for k = 1:2:pop_size-1
        if rand < crossover_rate
            alpha = rand;  % 線性交配比例
            temp1 = x(k);
            temp2 = x(k+1);
            x(k) = alpha * temp1 + (1 - alpha) * temp2;
            x(k+1) = alpha * temp2 + (1 - alpha) * temp1;
        end
    end

    % === Mutation（突變率控制）===
    for i = 1:pop_size
        if rand < mutation_rate
            delta = ((rand * 2) - 1) * 0.1;  % [-0.1, 0.1]
            x(i) = x(i) + delta;
        end
    end

    % === 最差替換 ===
    fitness = -15 * (sin(2 * x).^2) - (x - 2).^2 + 160;
    [~, idx_worst] = min(fitness);
    x(idx_worst) = best_individual;
end

% === 繪圖 ===
subplot(2,1,2); plot(1:n_gen, Fitness_f);
xlabel('Generation'); ylabel('Best Fitness');
title('每代最佳適應度');

% === 輸出最佳結果 ===
final_fitness = -15 * (sin(2 * x).^2) - (x - 2).^2 + 160;
[best_fit_value, idx] = max(final_fitness);
disp(['最佳解 x = ', num2str(x(idx))]);
disp(['對應適應度 = ', num2str(best_fit_value)]);
