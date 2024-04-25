function [edges0, N0] = EdgesN(P0)
% this function generates histogram versions of the Monte Carlo p-curves
edges = linspace(0,1, 401);
[N0,edges0] = histcounts(P0, edges, 'Normalization','pdf');
edges0 = edges0(2:end) - (edges0(2)-edges0(1))/2;
end