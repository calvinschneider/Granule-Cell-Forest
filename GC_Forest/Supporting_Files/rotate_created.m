% Unrotate Trees then Rotate Based on Mean Transverse and Longitudinal

function tree3 = rotate_created(input,R1s,R2s)

tree1 = input;
tree2 = cell(1,length(tree1));
tree3 = cell(1,length(tree1));

for i = 1:length(tree1)
    % Determine offset from origin
    offsetx = tree1{i}.X(1,1);
    offsety = tree1{i}.Y(1,1); 
    offsetz = tree1{i}.Z(1,1);

    % Move trees stem to origin
    tree1{i}.X = tree1{i}.X - offsetx;
    tree1{i}.Y = tree1{i}.Y - offsety;
    tree1{i}.Z = tree1{i}.Z - offsetz;

    unrotated_pts    = [tree1{i}.X tree1{i}.Y tree1{i}.Z] * transpose(R2s{i});
    unrotated_pts    = unrotated_pts * transpose(R1s{i});
    
    tree1{i}.X = unrotated_pts(:,1);
    tree1{i}.Y = unrotated_pts(:,2);
    tree1{i}.Z = unrotated_pts(:,3);    
    tree2{i} = rot_tree(tree1{i}, [90 0 0]);
    
    %{
    % Determine mean in transverse plane and rotate to y-axis
    mXY     = [mean(tree2{i}.X) mean(tree2{i}.Y)];
    mXY     = mXY / sqrt (sum (mXY.^2));
    vector  = [0 1];
    angle   = acosd(dot(mXY,vector)/(norm(mXY)*norm(vector)));
    if mean(tree2{i}.X)>0
        angle = -angle;
    end
    tree3{i} = rot_tree(tree2{i}, [0 0 angle]);
    %}
    
    % Rotate Based on Max X-Distance
    T = T_tree(tree2{i});
    Tx = tree2{i}.X(T);
    Ty = tree2{i}.Y(T);
    xy_slope  = (Ty(Tx == max(Tx))-Ty(Tx == min(Tx)))/(Tx(Tx == max(Tx)) - Tx(Tx == min(Tx)));
    xy_slope2 = -inv(xy_slope);
    XY      = [1 xy_slope2];
    XY      = XY / sqrt (sum (XY.^2));
    if xy_slope2 < 0 
        vector  = [0 -1];
    else
        vector  = [0 1];
    end
    xy_angle  = acosd(dot(XY,vector)/(norm(XY)*norm(vector)));
    if xy_slope2 < 0
        xy_angle = abs(xy_angle);
    else
        xy_angle = -abs(xy_angle);
    end
    tree3{i} = rot_tree(tree2{i}, [0 0 xy_angle]);
    
end