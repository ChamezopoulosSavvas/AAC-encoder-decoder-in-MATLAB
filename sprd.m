%% step 0
vars = {'B219a','B219b'};
load('Tableb219.mat', vars{:}); 

%incr first 3 columns by 1 so they match matlab array indexes
B219a(:,1:3) = B219a(:,1:3) + 1;
B219b(:,1:3) = B219b(:,1:3) + 1;

sprda = zeros(length(B219a));
sprdb = zeros(length(B219b));

for i=1:length(B219a)
    for j=1:length(B219a)
        sprda(i,j) = spreadingFunction(i,j, B219a(i,5), B219a(j,5));
    end
end

for i=1:length(B219b)
    for j=1:length(B219b)
        sprdb(i,j) = spreadingFunction(i,j, B219b(i,5), B219b(j,5));
    end
end

vars = {'sprda','sprdb'};
save('spreads.mat', vars{:});