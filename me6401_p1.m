home_theta_list=[1.1;0.8;0.5;0.4;0.2;0.5;0];

all_screw_axis = transpose([0 0 1 0 0 0;
                  0 1 0 -0.333 0 0;
                  0 0 1 0 0 0;
                  0 -1 0 0.649 0 -0.088;
                  0 0 1 0 0 0;
                  0 -1 0 1.033 0 0;
                  0 0 -1 0 0.088 0]);

zero_T = [[1 0 0; 0 -1 0; 0 0 -1] [0.088;0;0.926]; 0 0 0 1];

Theta_0 = home_theta_list;
space_T0 = space_forward_kinematics(all_screw_axis,Theta_0,zero_T);
space_J0 = space_jacobian(all_screw_axis,Theta_0);
Z0=space_T0(1:3,3);
Qz = acos(Z0(3))/pi*180;
if (Z0(2) >=0)
    Qx = (acos(Z0(1)))/pi*180;
else 
    Qx = 2*pi - acos(Z0(1))/pi*180;
end

% 1.2m max radius, 1x1x1cm boxes, 45 degree per quadrant
raw_map(1:240,1:240,1:240,1:32) = 0;

count = 0;
visited_count = 0;
% joint limits
% (-166,166),(-101,101),(-166,166),(-176,-4),(-166,166),(-1,215),(-166,166)
while(count < 1000000000 && visited_count < 100000)
    k = rand(7,1);
    theta_cur = k .* [332; 202; 332; 172; 332; 216; 332] + [-166; -101; -166; -176; -166; -1; -166];
    theta_cur = theta_cur / 180 * pi;
    space_Tcur = space_forward_kinematics(all_screw_axis,theta_cur,zero_T);
    %space_Jcur = space_jacobian(all_screw_axis,Theta_cur);
    
    % Find approach direction
    Z0=space_Tcur(1:3,3);
    Qz = acos(Z0(3))/pi*180;
    if (Z0(2) >=0)
        Qx = (acos(Z0(1)))/pi*180;
    else 
        Qx = 360 - acos(Z0(1))/pi*180;
    end

    position = space_Tcur(1:3,4) * 100;
    position = fix(position/1)+120; % Change to centimeter and integer only, origin at 120,120,120
    polar_quadrant = fix(Qz/45);
    azimuthal_quadrant = fix(Qx/45) + 1;
    if polar_quadrant == 4
        final_quadrant = polar_quadrant * 8;
    else
        final_quadrant = polar_quadrant * 8 + azimuthal_quadrant;
    end

    if raw_map(position(1),position(2),position(3),final_quadrant) == 0
        raw_map(position(1),position(2),position(3),final_quadrant) = 1;
        visited_count = 0;
    else
        visited_count = visited_count + 1;
    end

    count = count + 1;
    if mod(count,10000)==9999
        disp(count);
    end
end

final_map = sum(raw_map,4)/32; % Reachibility score at each position
save('final_map.mat','final_map');

% volumeViewer

function R3_matrix = R3_to_matrix(W)
    R3_matrix = [0 -W(3) W(2); W(3) 0 -W(1); -W(2) W(1) 0];
end

% adjoint_T is the 6x6 adjoint matrix of 4x4 homogeneous transformation
%   matrix T
function adjoint_T = adjoint_T(T)
    R = T(1:3,1:3);
    P = T(1:3,4);
    adjoint_T = [R eye(3)*0;
                 R3_to_matrix(P)*R R];
end


function exp_of_screw = exp_screw(screw_axis,theta)
    omega = screw_axis(1:3);
    v = screw_axis(4:6);
    if sqrt(sum(omega.*omega)) == 1
        omega_mat = R3_to_matrix(screw_axis(1:3));
        G_theta = eye(3)*theta + (1-cos(theta))*omega_mat + (theta - sin(theta))*omega_mat*omega_mat;
        exp_of_screw = [expm(omega_mat*theta) G_theta*v; 0 0 0 1];
    end

    if sqrt(sum(omega.*omega)) == 0 && sqrt(sum(v.*v)) == 1
        exp_of_screw = [eye(3) theta*v; 0 0 0 1];
    end
end

% S_list is 6xn, where each column is the screw axis for each link
% theta_list = nx1, each row is the joint displacement from 'zero' position
% zero_T is the EE transformation when at 'zero' position
function space_forward_T = space_forward_kinematics(S_list,theta_list,zero_T)
    n = length(theta_list);
    space_forward_T = eye(4);
    for k = 1:n
        space_forward_T = space_forward_T * exp_screw(S_list(:,k),theta_list(k));
    end
    space_forward_T = space_forward_T * zero_T;
end

% space_J is 6xn matrix, multiplying with nx1 theta_dot_vector to give
%  space 6x1 twist of EE.
function space_J = space_jacobian(S_list,theta_list)
    n = length(theta_list);
    J1 = S_list(:,1);
    space_J = J1;
    for k = 2:n
        T = eye(4);
        for i = 1:(k-1)
            T = T*exp_screw(S_list(:,i),theta_list(i));
        end
        space_J = [space_J adjoint_T(T)*S_list(:,k)];
    end
end

