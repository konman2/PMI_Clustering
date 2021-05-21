% Compute stability-based communities for both weighted and unweighted
% graphs. 
% Input: 
%   filepath: the folder where edgelist.csv resides. 
%   time: the Markov time range 
% Output: save the communities found at all markov time. 


function run_discrete_stability(filepath, time)
% Time: 1:100 (100 pts) for linspace. 
% Corresponding numpy: range(1, 101)
tic
% Default parameters.
if ~exist('filepath', 'var')
    % filepath = '/Users/zexi/Codes/PyCharmProjects/MarkovStability/graphs/airport_w/';
    filepath = 'E:\Synchronized_Files\Codes\PyCharmProjects\MarkovStability\graphs\blog\';
    % filepath = '/Users/zexi/Codes/PyCharmProjects/MarkovStability/graphs/snp500ll/preprocess/alpha=0.55/';
end
if ~exist('time', 'var')
    time = 1:100;
end

input_file = [filepath,'edgelist.csv'];
edgelist = load(input_file);

% Add symmetric edges in the graph to make it undirected. 
% Mark all weights to be 1 for unweighted graph. 
% Relable the node ids from 0.
[m,n] = size(edgelist);
sym_edgelist = zeros(2*m,3);
if n == 2 % unweighted:    
    for i = 1:m
        sym_edgelist(i,:) = [edgelist(i,:)-1,1];
        sym_edgelist(i+m,:) = [fliplr(edgelist(i,:))-1,1];
    end
else
    for i = 1:m
        sym_edgelist(i,:) = [edgelist(i,1:2)-1, edgelist(i,3)];
        sym_edgelist(i+m,:) = [fliplr(edgelist(i,1:2))-1,edgelist(i,3)];
    end
end

    
% Run stability. 
[S, N, VI, C] = stability_discrete(sym_edgelist, time, 'noVI', 'plot','full');

% Save the results.
% result = [time;S;N];
% table_full = array2table(result');
% writetable(table_full,[filepath,'stability_1e-5_1e5.csv'],'Delimiter',' ','WriteVariableNames',false);

% Save the communities. 
save([filepath, 'stability_clustering_tau_1_100_prediction.mat'], 'C')

toc

