function plotWobblingEnergies(EnergyData, MAX_ANGLE, MIN_LENGTH, MAX_LENGTH, ITERATIONS)
%EnergyData should be a matrix where each row is of
% the form "middleLength, angle, eigenvalue, hf_energy"
%And the procession of the data goes:
%   for each length:
%       for each angle:
%           middleLength angle eigenvalue hf_energy

angles = zeros(MAX_ANGLE,1);
lengths = zeros(ITERATIONS,1);

datapoint = 1;
EigenValues = zeros(MAX_ANGLE, ITERATIONS);
HartreeFockEnergies = zeros(MAX_ANGLE, ITERATIONS);
for i = 1:MAX_ANGLE
    angles(i) = (i-1)*10;
    for j = 1:ITERATIONS
        lengths(j) = MIN_LENGTH+(j-1)*(MAX_LENGTH-MIN_LENGTH)/ITERATIONS;
        EigenValues(i, j) = EnergyData(datapoint, 3);
        HartreeFockEnergies(i, j) = EnergyData(datapoint, 4); 
        datapoint = datapoint + 1;
    end
end

% Surface plot
figure;
surf(lengths, angles, EigenValues)
hold on;
surf(lengths, angles, HartreeFockEnergies)
xlabel('length'); ylabel('angle');

% Contour plot
%figure;
% Plot J_vals as 15 contours spaced logarithmically between 0.01 and 100
%contour(lengths, angles, EigenValues, logspace(-2, 3, 20))
%xlabel('length'); ylabel('angle');
%hold on;
%plot(theta(1), theta(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end