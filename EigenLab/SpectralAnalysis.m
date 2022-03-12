%Read data
filename = 'DS19hH2_dk0_FTIR_Spectra_instant_coffee.csv';
data = csvread(filename, 1.0);
fn = 'EigenVectors.csv';
eigenvectors = csvread(fn);

projections = zeros(56, 8);

%Project spectra onto 8 Eigenvectors
for i=1:8
    for j=1:56
        proj = dot(eigenvectors(:, i), data(j, :));
        projections(j, i) = proj;
    end
end

%I forgot you could use a for loop so yepp. No time to make it better now
%but it works so it's fine
figure();
tiledlayout(7, 4)
nexttile
hold on;
scatter(projections(1:29, 1), projections(1:29, 2), 'or');
scatter(projections(30:56, 1), projections(30:56, 2), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 1), projections(1:29, 3), 'or');
scatter(projections(30:56, 1), projections(30:56, 3), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 1), projections(1:29, 4), 'or');
scatter(projections(30:56, 1), projections(30:56, 4), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 1), projections(1:29, 5), 'or');
scatter(projections(30:56, 1), projections(30:56, 6), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 1), projections(1:29, 6), 'or');
scatter(projections(30:56, 1), projections(30:56, 6), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 1), projections(1:29, 7), 'or');
scatter(projections(30:56, 1), projections(30:56, 7), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 1), projections(1:29, 8), 'or');
scatter(projections(30:56, 1), projections(30:56, 8), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 2), projections(1:29, 2), 'or');
scatter(projections(30:56, 2), projections(30:56, 2), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 2), projections(1:29, 3), 'or');
scatter(projections(30:56, 2), projections(30:56, 3), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 2), projections(1:29, 4), 'or');
scatter(projections(30:56, 2), projections(30:56, 4), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 2), projections(1:29, 5), 'or');
scatter(projections(30:56, 2), projections(30:56, 5), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 2), projections(1:29, 6), 'or');
scatter(projections(30:56, 2), projections(30:56, 6), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 2), projections(1:29, 7), 'or');
scatter(projections(30:56, 2), projections(30:56, 7), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 2), projections(1:29, 8), 'or');
scatter(projections(30:56, 2), projections(30:56, 8), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 3), projections(1:29, 3), 'or');
scatter(projections(30:56, 3), projections(30:56, 3), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 3), projections(1:29, 4), 'or');
scatter(projections(30:56, 3), projections(30:56, 4), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 3), projections(1:29, 5), 'or');
scatter(projections(30:56, 3), projections(30:56, 5), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 3), projections(1:29, 6), 'or');
scatter(projections(30:56, 3), projections(30:56, 6), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 3), projections(1:29, 7), 'or');
scatter(projections(30:56, 3), projections(30:56, 7), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 3), projections(1:29, 8), 'or');
scatter(projections(30:56, 3), projections(30:56, 8), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 4), projections(1:29, 4), 'or');
scatter(projections(30:56, 4), projections(30:56, 4), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 4), projections(1:29, 5), 'or');
scatter(projections(30:56, 4), projections(30:56, 5), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 4), projections(1:29, 6), 'or');
scatter(projections(30:56, 4), projections(30:56, 6), 'ob');
nexttile;
hold on;
scatter(projections(1:29, 4), projections(1:29, 7), 'or');
scatter(projections(30:56, 4), projections(30:56, 7), 'ob');