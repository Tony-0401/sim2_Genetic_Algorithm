%%
clc; clear; close all;

X = -10:0.01:10;
Y = -15*(sin(2*X).^2)-((X-2).^2) + 160;

% === GA 參數 ===
pop_size = 10;
gene_len = 10;
n_gen = 1000;

crossover_rate = 0.8;   % 交配率
mutation_rate = 0.01;   % 突變率
selection_type = 2;     % 1 = Tournament, 2 = Roulette Wheel

% 初始化族群
bit = randi([0 1], pop_size, gene_len);
decode = @(b) sum(b .* 2.^(gene_len-1:-1:0), 2) * 20 / (2^gene_len - 1) - 10;
x = decode(bit);

FF = zeros(1, n_gen);
xx = zeros(1, n_gen);

for i = 1:n_gen
    fitness = -15*(sin(2*x).^2) - (x-2).^2 + 160;

    % 精英保存
    [maxFit, idx_best] = max(fitness);
    best_x = x(idx_best);
    best_gene = bit(idx_best, :);
    FF(i) = maxFit;
    xx(i) = best_x;

    % === 選擇機制 ===
    new_bit = zeros(size(bit));
    switch selection_type
        case 1  % Tournament
            for k = 1:pop_size
                c = randperm(pop_size, 2);
                [~, winner] = max(fitness(c));
                new_bit(k,:) = bit(c(winner), :);
            end
        case 2  % Roulette Wheel
            fit_sum = sum(fitness);
            prob = cumsum(fitness / fit_sum);
            for k = 1:pop_size
                r = rand;
                selected = find(prob >= r, 1, 'first');
                new_bit(k,:) = bit(selected, :);
            end
    end
    bit = new_bit;

    % === 單點交配（交配率控制） ===
    for k = 1:2:pop_size-1
        if rand < crossover_rate
            point = randi([2, gene_len-1]);
            bit([k k+1], point:end) = bit([k+1 k], point:end);
        end
    end

    % === 突變（每個基因位元以突變率變異） ===
    for k = 1:pop_size
        for j = 1:gene_len
            if rand < mutation_rate
                bit(k,j) = ~bit(k,j);
            end
        end
    end

    x = decode(bit);
    fitness = -15*(sin(2*x).^2) - (x-2).^2 + 160;
    [~, idx_worst] = min(fitness);
    x(idx_worst) = best_x;
    bit(idx_worst,:) = best_gene;
end

% === 繪圖與結果 ===
figure;
subplot(2,1,1);
plot(X, Y);
xlabel('X'); ylabel('Fitness value'); title('Fitness Function');

subplot(2,1,2);
plot(FF);
xlabel('Generation'); ylabel('Best Fitness'); title('Best Fitness per Generation');

[~, best_gen] = max(FF);
disp(['最佳解 x = ', num2str(xx(best_gen))]);
disp(['對應適應度 = ', num2str(FF(best_gen))]);
