%project1
%author: Michael Luo

clear %clear the memory
n = 2; %number of bodies
g = 0.2; %gravitational constant
m = [1000, 1000]; %masses
x = [-5, 5]; %initial positions
y = [0, 0];
z = [0, 0];
u = [0, 0];
v = [2, -2];
%calculate and subtract the velocity of the center of mass as a constant
Vcm = ((m(1)*v(1)+m(2)*v(2))/(m(1)+m(2)));
v = [2 - Vcm, -2 - Vcm];
w = [0, 0];
tmax = 100; %length of simulation
dt = 0.005; %time step
clockmax = ceil(tmax/dt);

%v = [sqrt(g * m(1)/(x(2)-(x(1)))), -sqrt(g * m(2)/(x(2)-(x(1))))]; %initial velocity

distanceX = x(2) - x(1);
distanceY = y(1) - y(2);
distance = sqrt(distanceX^2 + distanceY^2);


%initialization for plotting
%set (gcf, 'double plotting', 'or')
subplot(3,2,1), hxy = plot(x,y,'bo');
    axis ([-5 5 -5 5])
    %axis ([-10 10 -10 10 -10 10])
    axis equal
subplot(3,2,2), hzy = plot(z,y,'bo');
    axis ([-5 5 -5 5])
    %axis ([-10 10 -10 10 -10 10])
    axis equal
subplot(3,2,3), hxz = plot(x,z,'bo');
    axis ([-5 5 -5 5])
    %axis ([-10 10 -10 10 -10 10])
    axis equal
for clock = 1: clockmax
    set (hxy,'xdata',x,'ydata',y)
    set (hzy,'xdata',z,'ydata',y)
    set (hxz,'xdata',x,'ydata',z)
    t = clock*dt;
    %summation for the mathematical methods is implicit in the for loops
    for i = 1:n
        for j = 1:n
            if (j ~= i)
                %change in x y and z
                dx = x(j)-x(i);
                dy = y(j) - y(i);
                dz = z(j) - z(i);
                %the radius or r_ij
                r = sqrt(dx^2 + dy^2 + dz^2);
                %calculating respective velocities
                u(i) = u(i) + dt*g*m(j)*dx/r^3;
                v(i) = v(i) + dt*g*m(j)*dy/r^3;
                w(i) = w(i) + dt*g*m(j)*dz/r^3;
            end
        end
    end
    for i = 1:n
        %new positions
        x(i) = x(i) + dt*u(i);
        y(i) = y(i) + dt*v(i);
        z(i) = z(i) + dt*w(i);
        %distance values to be used for graphing and determining
        %the maximum radius vectors
        distanceX = x(2) - x(1);
        distanceY = y(1) - y(2);
        distance(i) = sqrt(distanceX^2 + distanceY^2);
    end
    
    subplot(3,2,4), hxyz = plot3(x,y,z,'bo');
    axis ([-5 5 -5 5 -5 5])
    %axis ([-15 15 -15 15 -10 15])
    set (hxyz,'xdata',x,'ydata',y, 'zdata', z)
    
    subplot(3,2,5), hxyz = plot(x,distance,'b.');
    axis ([-5 5 -5 5 -5 5])
    %axis ([-15 15 -15 15 -10 15])
    set (hxyz,'xdata',x,'ydata',y, 'zdata', z)
    hold on
    drawnow
end