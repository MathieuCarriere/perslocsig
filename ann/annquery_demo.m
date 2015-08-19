function annquery_demo()
%ANNQUERY_DEMO A demo to show how to search and plot with ANN Wrapper
%
% [ History ]
%   - Created by Dahua Lin, on Aug 10, 2007
%

%%  Prepare Data Points
%ref_pts = rand(2, 300);
ref_pts = [1 2 3 4 5 6; 5 1 3 6 2 4];
%query_pts = rand(2, 100);
query_pts = [1 4 2; 7 5 9];

%% Do ANN query
k = 1;
[nnidx,dists] = annquery(ref_pts, query_pts, k);

nnidx
dists.^2

%% Plot the results
anngplot(ref_pts, query_pts, nnidx);
axis([0 10 0 10]);
