%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Control of Redundant manipulator (planar 2 dof) 
%%%%% in the presence of moving obstacles, using a method 
%%%%% that involves the maximization of a criterion (distance)
close all;
clear all;
%clear q*; clear x*; clear y*; clear H*;

%% Parameters initialisation
kc = 15;  %%gain for control of redundant degrees of freedom (joints' null space)
%kc= 0;   %%no consideration of redundacies for collision avoidance
%kc= 20;
v_obst = 1.6; %1.5;

% Lengths of robot links
l(1) = 1.0;  %% in dm
l(2) = 1.0;


%% *** sampling period ***
%% *** for the robot motion, kinematic simulation:
dt = 0.001; %dt = 0.001; i.e. 1 msec) 

%% *** Create (or load from file) reference signals ***
%% *** DESIRED MOTION PROFILE - TASK SPACE ***
Tf=1; 	%  duration of motion (in secs)
t=0:dt:Tf;
kmax=Tf/dt + 1;

%%Initial configuration
q0(1) = 60*pi/180; q0(2) = -50*pi/180;
pobst_y(1) = 2;
pobst_x(1) = 0.75;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the 1st priority subtask: 
% xd = 1.5 = constant
xd = 1.5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2nd subtask: will be defined according to a criterion V(q)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ****** KINEMATIC SIMULATION - Main loop ******
disp('Kinematic Simulation ...'); %%
disp(' '); %%


%% ***** INVESRE KINEMATICS  -->  DESIRED MOTION - JOINT SPACE *****
%% compute the reference joint-motion vectors:
%% {qd(k,i), i=1,...,n (num of degrees of freedom), with k=1,..., kmax,}
%% and reference joint (angular) velocities {qd_1(k,i)}


%%% Inverse geometric model %%%%
%rd2 = xd(:).^2 + yd(:).^2;
%%qd(:,2) = acos( (rd2(:)-l(1)^2-l(2)^2)./(2*l(1)*l(2)) ); %% 1st solution: elbow down
%qd(:,2) = -acos( (rd2(:)-l(1)^2-l(2)^2)./(2*l(1)*l(2)) ); %% 2nd solution: elbow up
%ss2 = sin(qd(:,2));
%qd(:,1) = atan2(yd(:),xd(:)) - asin(l(2)*ss2(:)./sqrt(rd2(:)));


tt=0; tk=1;
qd(tk,1)=q0(1); qd(tk,2)=q0(2);
while (tt<=Tf)

   %% Forward geometric model
   s1 = sin(qd(tk,1)); s12 = sin(qd(tk,1)+qd(tk,2));
   c1 = cos(qd(tk,1)); c12 = cos(qd(tk,1)+qd(tk,2));
   p1x(tk) = l(1)*c1; 
   p1y(tk) = l(1)*s1;
   p2x(tk) = p1x(tk) + l(2)*c12; 
   p2y(tk) = p1y(tk) + l(2)*s12;
 
   %% Compute Jacobian
   Jac1(1,2) = - l(2)*s12;
   Jac1(1,1) = Jac1(1,2) - l(1)*s1;
   
   %%Pseudo-inverse computation of Jacobian matrix Jac1
   Jac1_psinv = Jac1'*inv(Jac1*Jac1');
   
   %% ---- Subtask 1 ----
   kp1 = 10;
   xd_(tk) = kp1*(xd-p2x(tk));
   task1=Jac1_psinv*[xd_(tk)];
   
   %% ---- Subtask 2 ----
   %% Motion of the obstacle
   pobst_y(tk+1) = pobst_y(tk) - v_obst*dt;
   pobst_x(tk+1) = pobst_x(tk);
   %% Distance criterion --> gradient vector (dqc) computation
   d1 = pobst_y(tk) - p1y(tk);
   dqc(1) = -d1*l(1)*c1;
   dqc(2) = 0;
   H2=eye(2)*kc;
   task2 = (eye(2)-Jac1_psinv*Jac1)*H2*dqc';
   nn=norm(task2);
   if (nn>15), 
       task2n(tk,:) = task2/nn; 
       task2 = task2n(tk,:)';
   end
   
   % debug:
   crit(tk) = dqc(1);
   ctask(tk,:) = task2;
   
   %angular velocity
   qd_(tk,:) = task1'+task2';
   
   %% numerical integration --> kinematic simulation
   qd(tk+1,1) = qd(tk,1) + dt*qd_(tk,1);
   qd(tk+1,2) = qd(tk,2) + dt*qd_(tk,2);
   
   
   tk=tk+1;
   tt=tt+dt;
end



%% ***** FORWARD KINEMATICS  JOINT MOTION -->  CARTESIAN POSITIONS *****

%%(xd1, yd1) : cartesian position of the 1st link's local reference frame
xd1 = l(1)*cos(qd(:,1));  
yd1 = l(1)*sin(qd(:,1));

%%(xd2, yd2) : cartesian position of the 2nd link's local reference frame
xd2 = xd1 + l(2)*cos( qd(:,1)+qd(:,2) );  
yd2 = yd1 + l(2)*sin( qd(:,1)+qd(:,2) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% *** SAVE and PLOT output data ***
%%** use functions plot(...)
save;  %% --> save data to 'matlab.mat' file

fig1 = figure;

subplot(2,1,1);
plot(t,qd(:,1));
ylabel('qd1 (rad)');
xlabel('time t (sec)');


subplot(2,1,2);
plot(t,qd(:,2));
ylabel('qd2 (rad)');
xlabel('time t (sec)');



%%*** stick diagram --> animate robot motion ... (**optional**)
%% within a for (or while) loop, use periodic plot(...) functions to draw the geometry (current pos) 
%% of the robot, and thus animate its motion ...

fig2 = figure;
axis([-1 3 -2 2]) %%set xy plot axes (caution: square axes, i.e. dx=dy)
axis on
hold on

xlabel('x (cm)');
ylabel('y (cm)');
%plot(xd,yd,'rs');
plot([xd,xd], [-1.5,1.5],'r:');
dtk=50; %% plot robot position every dtk samples, to animate its motion


plot([0],[0],'o');

%% OBSTACLE
%fill([0.4,1.1,1.1,0.4],[0.0,0.0,-0.4,-0.4],[0.8,0.3,0.3]);

for tk=1:dtk:kmax,

   %%% --- redraw each frame --- %%%%
   clf;
   axis([-1 3 -2 2]) %%set xy plot axes (caution: square axes, i.e. dx=dy)
    axis on
    hold on

    xlabel('x (cm)');
    ylabel('y (cm)');
    %plot(xd,yd,'rs');
    plot([xd,xd], [-1.5,1.5],'r:');
    dtk=50; %% plot robot position every dtk samples, to animate its motion
    plot([0],[0],'o');
    %%% ---------------------

   plot([0,xd1(tk)],[0,yd1(tk)]);					
   plot([xd1(tk)],[yd1(tk)],'o');
   plot([xd1(tk),xd2(tk)],[yd1(tk),yd2(tk)]);
   plot([xd2(tk)],[yd2(tk)],'o');
   
   plot([pobst_x(tk)],[pobst_y(tk)],'*');
   %plot([xd(tk)],[yd(tk)],'g+');
    
	pause(0.5);	%% pause motion to view successive robot configurations
    
end





















