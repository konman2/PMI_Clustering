% Compute stability-based communities for a weighted graph.
% Input: 
%   filepath: the folder where edgelist.csv resides. 
%   time: the Markov time range 
% Output: save the communities found at all markov time. 


function run_weighted_stability(filepath, time)
% Time: 1:100 (100 pts) for linspace. 
% Corresponding numpy: range(1, 101)

% Default parameters.
if ~exist('filepath', 'var')
    filepath = 'E:\Synchronized_Files\Codes\PyCharmProjects\MarkovStability\graphs\airport_w\';
end
if ~exist('time', 'var')
    time = 1:100;
end

input_file = [filepath,'edgelist.csv'];
edgelist = load(input_file);

% Add symmetric edges in the graph to make it undirected. 
% Mark all weights to be 1. 
% Relable the node ids from 0.
m = length(edgelist);
sym_edgelist = zeros(2*m,3);
for i = 1:m
    sym_edgelist(i,:) = [edgelist(i,1:2)-1, edgelist(i,3)];
    sym_edgelist(i+m,:) = [fliplr(edgelist(i,1:2))-1,edgelist(i,3)];
end
% Run stability. 
[S, N, VI, C] = stability_discrete(sym_edgelist, time, 'noVI', 'plot','full');

% Save the communities. 
save([filepath, 'Stability_clustering_tau_1_100_prediction.mat'], 'C')

