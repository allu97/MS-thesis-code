% Predictor Importance score
% [idx, scores] = fscmrmr(X_ree(1:41,:), subclass_s_ree(1:41));
% [idx, scores] = fscmrmr(X_ree(41:end,:), subclass_s_ree(41:end,:));
% [idx, scores] = fscmrmr(X_ree, subclass_s_ree);
[idx, scores] = fscmrmr(X_ree, class_s_ree);

fig = figure;
bar(scores(idx))
xlabel('Predictor Rank')
ylabel('Predictor Importance Score')
title("Predictor Importace Ranking")
xticks(1:length(idx))
xticklabels(X_ree_var(idx))
% Idea: give weights to the main components of rapakivi granite
%% remove these from fingerprint step 2.
idx3 = [1 7 15 16 18];
X_ree(:,idx3) = '';
X_ree_var(idx3) = '';
%%
[coeff,score,latent,tsquared,explained] = pca(X_ree(42:end,:));

cumsum(explained)