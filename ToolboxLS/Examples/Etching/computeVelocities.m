function [velXGrid, velYGrid] = computeVelocities(g,init_level,current,i)
    global numberOfMaterials 
    [startMatInd endMatInd centMatInd map] = getMaterialMap(g);
    vel{1,1} = zeros(size(g.xs{1,1}));
    vel{2,1} = zeros(size(g.xs{1,1}));
    vy = [5 1]
    vx = [.3 .3]
    for m = 1:numberOfMaterials
        k = find(map==m);
        [velocityX, velocityY] = plasma(current,i)
        velocityY = vy(m);
        velocityX = vx(m);
        for q = 1:length(k)
            dy = abs(g.xs{2,1} - init_level);
            dx = abs(g.xs{1,1}-centMatInd(k(q)));
            lambda1 = current(9); %Make into current parameter to be optimized
            lambda2 = current(10);
              %make matrix with rows that decay exponentially symmetrically across 0
            assignin('base', 'lambdaInFunction', lambda2)
            v_sy = exp(-dy/lambda2)*-1*velocityY; %multiply by negative one to put in right direction
            v_sx = exp(-dx/lambda1)*velocityX;
            %make matrix with rows that decay exponentially symmetrically across 0
            assignin('base', 'dy', dy)
            assignin('base', 'dx', dx)
            vel{1,1}(g.xs{1,1}>=startMatInd(1,k(q)) & g.xs{1,1} <= centMatInd(k(q))& g.xs{2,1} >= startMatInd(2,k(q)) & g.xs{2,1} <= endMatInd(2,k(q))) = v_sx(g.xs{1,1}>=startMatInd(1,k(q)) & g.xs{1,1} <= centMatInd(k(q))&g.xs{2,1} >= startMatInd(2,k(q)) & g.xs{2,1} <= endMatInd(2,k(q)));
            vel{1,1}(g.xs{1,1}>=centMatInd(k(q)) & g.xs{1,1} <= endMatInd(1,k(q))& g.xs{2,1} >= startMatInd(2,k(q)) & g.xs{2,1} <= endMatInd(2,k(q))) = -v_sx(g.xs{1,1}>=centMatInd(k(q)) & g.xs{1,1} <= endMatInd(1,k(q))& g.xs{2,1} >= startMatInd(2,k(q)) & g.xs{2,1} <= endMatInd(2,k(q)));
            vel{2,1}(g.xs{2,1} >= startMatInd(2,k(q)) & g.xs{2,1} <= endMatInd(2,k(q))& g.xs{1,1}>=startMatInd(1,k(q)) & g.xs{1,1} <= endMatInd(1,k(q))) = v_sy(g.xs{2,1} >= startMatInd(2,k(q)) & g.xs{2,1} <= endMatInd(2,k(q))& g.xs{1,1}>=startMatInd(1,k(q)) & g.xs{1,1} <= endMatInd(1,k(q)));
        end
            assignin('base', 'v_sx', v_sx)
            assignin('base', 'v_sy', v_sy)
%     vel{1,1}(g.xs{1,1}>-init_level & g.xs{1,1} < 0) = -v_sx(g.xs{1,1}>-init_level & g.xs{1,1} < 0);
%     vel{1,1}(g.xs{1,1}>= 0 & g.xs{1,1} < init_level) = v_sx(g.xs{1,1}>= 0 & g.xs{1,1} < init_level);
%     vel{2,1}(g.xs{2,1} > -init_level) = v_sy(g.xs{2,1} > -init_level);
        assignin('base', 'v_sx_assigned', vel{1,1})
        assignin('base', 'v_sy_assigned', vel{2,1})
        assignin('base', 'dy', dy)
        assignin('base', 'init_level', init_level)
    gridSize = (g.xs{1,1});
    end
    velXGrid = vel{1,1};
    velYGrid = vel{2,1};
end
