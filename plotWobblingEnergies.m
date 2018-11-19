function plotWobblingEnergies(EnergyData, MIN_ANGLE, MAX_ANGLE, MIN_LENGTH, MAX_LENGTH, ITERATIONS)
%EnergyData should be a matrix where each row is of
% the form "middleLength, angle, eigenvalue, hf_energy"
%And the procession of the data goes:
%   for each length:
%       for each angle:
%           middleLength angle eigenvalue hf_energy
[max_iter, datacols] = size(EnergyData);
max_iter = max_iter/(MAX_ANGLE-MIN_ANGLE);
angles = zeros(MAX_ANGLE-MIN_ANGLE,1);
lengths = zeros(max_iter,1);

datapoint = 1;
EigenValues = zeros(max_iter, MAX_ANGLE-MIN_ANGLE);
HartreeFockEnergies = zeros(max_iter, MAX_ANGLE-MIN_ANGLE);

for i=1:max_iter
    lengths(i) = MIN_LENGTH+(i-1)*(MAX_LENGTH-MIN_LENGTH)/ITERATIONS;
    for j=1:MAX_ANGLE-MIN_ANGLE
        angles(j) = (j-1+MIN_ANGLE)*10;
        EigenValues(i, j) = EnergyData(datapoint, 3);
        HartreeFockEnergies(i, j) = EnergyData(datapoint, 4); 
        datapoint = datapoint + 1;
    end
end
        


%for i = 1:MAX_ANGLE
 %   angles(i) = (i-1)*10;
  %  for j = 1:ITERATIONS
   %     lengths(j) = MIN_LENGTH+(j-1)*(MAX_LENGTH-MIN_LENGTH)/ITERATIONS;
    %    EigenValues(i, j) = EnergyData(datapoint, 3);
     %   HartreeFockEnergies(i, j) = EnergyData(datapoint, 4); 
      %  datapoint = datapoint + 1;
    %end
%end

EigenValues = EigenValues';
HartreeFockEnergies = HartreeFockEnergies';
angles = angles';
lengths = lengths';

% Surface plot
figure;
surf(lengths, angles, EigenValues)
hold on;
surf(lengths, angles, HartreeFockEnergies)
xlabel('length'); ylabel('angle');zlabel('energy');

% Contour plot
%figure;
% Plot J_vals as 15 contours spaced logarithmically between 0.01 and 100
%contour(lengths, angles, EigenValues, logspace(-2, 3, 20))
%xlabel('length'); ylabel('angle');
%hold on;
%plot(theta(1), theta(2), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
end