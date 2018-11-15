function plotData(EnergyData)
%EnergyData should be a one col, vector

MAX_ANGLE = 3;
MIN_LENGTH = 0.2;
MAX_LENGTH = 4;
ITERATIONS = 25;

angles = zeros(MAX_ANGLE,1);
lengths = zeros(ITERATIONS,1);

datapoint = 1;
EnergyMatrix = zeros(MAX_ANGLE, ITERATIONS);
for i = 1:MAX_ANGLE
    angles(i) = (i-1)*10*3.14159/180;
    for j = 1:ITERATIONS
        lengths(j) = MIN_LENGTH+(j-1)*(MAX_LENGTH-MIN_LENGTH)/ITERATIONS;
        EnergyMatrix(i, j) = EnergyData(datapoint);
        %EnergyMatrix(i, j) = datapoint;
        datapoint = datapoint + 1;
    end
end

% Surface plot
figure;
surf(lengths, angles, EnergyMatrix)
xlabel('length'); ylabel('angle');

% Contour plot
figure;
% Plot J_vals as 15 contours spaced logarithmically between 0.01 and 100
contour(lengths, angles, EnergyMatrix, logspace(-2, 3, 20))
xlabel('length'); ylabel('angle');
%hold on;
%plot(theta(1), theta(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end