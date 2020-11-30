clear variables;
a = 3.441;  % a-axies of the lattice
b = a;      % b-axies of the lattice
gamma = 120;% angle of <a, b>
sizeLattice = 51;  % Scale of the system (should be odd)

%% Initialize of the variables
centerOrder = (sizeLattice + 1) ./ 2;
lattice.x = zeros(sizeLattice);
lattice.y = zeros(sizeLattice);
lattice.distance = zeros(sizeLattice);
% Difference in x and y in the primitive cell
deltaA = [a, 0];
deltaB = [b*cosd(gamma), b*sind(gamma)];
% Set the zero of the axies
lattice.x(sizeLattice, 1) = 0;
lattice.y(sizeLattice, 1) = 0;
% Initialize the position of the center
lattice.x(centerOrder, centerOrder) = ((centerOrder - 1)*deltaA(1) + (sizeLattice - centerOrder)*deltaB(1));
lattice.y(centerOrder, centerOrder) = ((centerOrder - 1)*deltaA(2) + (sizeLattice - centerOrder)*deltaB(2));
%% Calculate the distance from the center
for row = sizeLattice: -1: 1
    for column = 1: sizeLattice
        lattice.x(row, column) = ((column - 1)*deltaA(1) + (sizeLattice - row)*deltaB(1));
        lattice.y(row, column) = ((column - 1)*deltaA(2) + (sizeLattice - row)*deltaB(2));
        lattice.hoppingA(row, column) = column - centerOrder;
        lattice.hoppingB(row, column) = centerOrder - row;
        lattice.distance(row, column) = sqrt((lattice.x(row, column) - lattice.x(centerOrder, centerOrder)).^2 + (lattice.y(row, column) - lattice.y(centerOrder, centerOrder)).^2);
    end
end
%% Get the N^(th) neighbour
output.distanceTemp = lattice.distance;
i = 1;
while min(min(output.distanceTemp)) ~= Inf
    minValue = min(min(output.distanceTemp));
    output.neighbour{i} = (output.distanceTemp - minValue) < 1e-6;
    output.neighbourDistance{1, i} = output.distanceTemp(output.neighbour{i});
    output.neighbourDistance{2, i} = length(output.neighbourDistance{1, i});
    output.hopping{i}(:, 1) = lattice.hoppingA(output.neighbour{i});
    output.hopping{i}(:, 2) = lattice.hoppingB(output.neighbour{i});
    output.distanceTemp(output.neighbour{i}) = Inf;
    i = i+1;
end

%% Plot the neighour
color = lines(2^10);
fig = figure;
limit = [lattice.x(centerOrder, centerOrder) - 50, lattice.x(centerOrder, centerOrder) + 50, lattice.y(centerOrder, centerOrder) - 50, lattice.y(centerOrder, centerOrder) + 50];
scat = scatter(lattice.x(:), lattice.y(:), [], ones(length(lattice.x).^2, 3), 'filled');
axis(limit);
fig.PaperPositionMode = 'auto';
for i = 1: length(output.neighbour)
    scat.CData(output.neighbour{i}, :) = repmat(color(i, :), length(scat.CData(output.neighbour{i}, 1)), 1);
    %saveas(fig, ['C:\Users\SCES\Desktop\temp\img\', num2str(i),'.png']);
end