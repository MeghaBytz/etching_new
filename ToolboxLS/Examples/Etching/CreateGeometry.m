
R1 = [3,4,-1,1,1,-1,0.5,0.5,-0.75,-0.75]';
R2 = [3,4,-0.2,0.2,.2,-.2,0.5, 0.5, 0,0]';
gm = [R1,R2];
sf = 'R1-R2';
ns = char('R1','R2');
ns = ns';
g = decsg(gm,sf,ns);
figure;
pdegplot(g);
model = createpde;
geometryFromEdges(model,g);
mesh = generateMesh(model);
pdeplot(model);

%compute viewfactors of all mesh triangles from coordinates above
vfactorinfo = zeros(length(mesh.Elements),4);
%vfactor info = [vfactor12,vfactor21,area1,area2]
for i=1:length(mesh.Elements)%get coordinates of nodes attached to each element
    elem = mesh.Elements(:,i);
    nodeCoord1 = [mesh.Nodes(:,elem(1))' 0];
    nodeCoord2 = [mesh.Nodes(:,elem(2))' 0];
    nodeCoord3 = [mesh.Nodes(:,elem(3))' 0];
   [vfactorinfo(i,1),vfactorinfo(i,2),vfactorinfo(i,3),vfactorinfo(i,4)]=viewfactor([nodeCoord1;nodeCoord2;nodeCoord3],[[-1,0.6,0];[0,0.8,0];[1,0.6,0]],6)
end
    

%[vfactor12,vfactor21,area1,area2]=viewfactor([[0,0,0];[1,0,0];[0,1,0]],[[0,0,1];[2,0,1];[0,2,1]],6)