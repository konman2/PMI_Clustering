% Generate stability for both unweighted and weighted graph.
% Input: 
%   filepath: the folder where edgelist.csv resides. 
%   time: the Markov time range 
% Output: save a stability.csv in the filepath folder with space-delimited 
% time, stability, # of communities. 

function run_stability(filepath, time)
% Time: 10.^[-3:0.1:3] (61 pts) for logspace and 1:100 (100 pts) for linspace. 
% Corresponding numpy: np.logspace(-3,3,61) and np.arange(1,101).

% Default parameters.
if ~exist('filepath', 'var')
    % filepath = '/Users/zexi/Codes/PyCharmProjects/MarkovStability/graphs/airport_w/';
    filepath = '../Graphs/airport_ww/';
    % filepath = '/Users/zexi/Codes/PyCharmProjects/MarkovStability/graphs/snp500ll/preprocess/alpha=0.55/';
end
if ~exist('time', 'var')
    time = 10.^[-5:0.1:5];
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
[S, N, VI, C] = stability(sym_edgelist, time, 'noVI', 'plot','full');

% Save the results.
result = [time;S;N];
table_full = array2table(result');
writetable(table_full,[filepath,'stability_1e-5_1e5.csv'],'Delimiter',' ','WriteVariableNames',false);

% Also save the results for between 1e-5 to 1e1.5
% table_650 = array2table(result(:,1:651)');
% writetable(table_650,[filepath,'stability_650.csv'],'Delimiter',' ','WriteVariableNames',false);

