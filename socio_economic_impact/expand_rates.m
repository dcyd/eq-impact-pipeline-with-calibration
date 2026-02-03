function rates = expand_rates(hazus_idx, R)
% 把每栋楼的 Hazus 索引（uint8/16; 0=未映射）展开为 20 列伤亡率
% R: 36×20 伤亡率矩阵（与 hazusOrder 对齐）
    n = numel(hazus_idx);
    rates = nan(n, size(R,2));
    v = hazus_idx > 0;
    idx = double(hazus_idx(v));      % 1..36
    rates(v,:) = R(idx, :);
end