function  [mask,g] = maskTetraPore1(g)

% function  [mask,g] = maskTetraPore1(g)

h1 = 0.25;

%throat radii
r1 = 0.14;
r2 = 0.14;
r3 = 0.14;
r4 = 0.14;

p(:,1) = [g.min(1);r1];
p(:,2) = [-h1;h1];
p(:,3) = [-r2; g.max(2)];

[mask,g] = maskCorner(g,p);

p(:,1) = [r2; g.max(2)];
p(:,2) = [h1;h1];
p(:,3) = [g.max(1);r3];

[mask1,g] = maskCorner(g,p);
mask = max(mask,mask1);

p(:,1) = [g.max(1);-r3];
p(:,2) = [h1;-h1];
p(:,3) = [r4;g.min(2)]

[mask1,g] = maskCorner(g,p);
mask = max(mask,mask1);

p(:,1) = [-r4;g.min(2)];
p(:,2) = [-h1;-h1];
p(:,3) = [g.min(1); -r1];

[mask1,g] = maskCorner(g,p);
mask = max(mask,mask1);