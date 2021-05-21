% Some used codes for processing.
% To be formalized later.

% The result is suitable for gephi drawing.
table = [Time;C];
ids = [0, 1:108]';
table = [ids table];
T = array2table(table);
writetable(T,'hier18_coms.csv');

% This generates the markov_time.csv and partitions.csv
table = array2table(C');
writetable(table,'partitions.csv');
table = array2table(Time);
writetable(table,'markov_time.csv');





