%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% KINEMATIC CONTROL of REDUNDANT ROBOT MANIPULATOR %%%%%
%%%%%  (PLANAR 3 DOF MANIPULATOR WITH 1 REDUNDANT DOF) %%%%%

close all;
clear all;

%%Gain for control of redundant degrees of freedom 
kc = 20;  
%kc= 0;   
%kc= 20;

%% Length of the links
l(1) = 1.0;  %% in dm
l(2) = 1.0;
l(3) = 0.3;


%% *** Sampling period ***
%% *** for the robot motion, kinematic simulation:
dt = 0.001; %dt = 0.001; i.e. 1 msec) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% *** Create (or load from file) reference signals ***
%% *** DESIRED MOTION PROFILE - TASK SPACE ***
Tf=1.0; 	%  duration of motion (in secs)
t=0:dt:Tf;


% Example of desired trajectory : linear segment (x0,y0)-->(xf,yf); Time duration: Tf;
disp('Initialising Desired Task-Space Trajectory (Motion Profile) ...'); %%
disp(' '); 

%% *** Initial configuration *** 
qd1_0 = 20*pi/180; qd2_0 = 30*pi/180; qd3_0 = 20*pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *** 1st subtask: desired trajectory *** 
% *** (xd2_0,yd2_0),(xd3_0,yd3_0),(xde_0,yde_0): initial positions 
%     for positions of joints 2, 3, and end-effector, respecitvely
% *** (xde_f,yde_f): final end-effector position
s1 = sin(qd1_0); 
s12 = sin(qd1_0+qd2_0); 
s123 = sin(qd1_0+qd2_0+qd3_0);
c1 = cos(qd1_0); 
c12 = cos(qd1_0+qd2_0); 
c123 = cos(qd1_0+qd2_0+qd3_0);

xd2_0=l(1)*c1;	%%% initial positions
yd2_0=l(1)*s1;
xd3_0=l(1)*c1+l(2)*c12;
yd3_0=l(1)*s1+l(2)*s12;
xde_0=l(1)*c1+l(2)*c12+l(3)*c123;
yde_0=l(1)*s1+l(2)*s12+l(3)*s123;

xde_f=xde_0;	%%% final position
yde_f=0;

%% desired trajectory (1st subtask): position (xd,yd); velocity (xd_,yd_) 
lambda_x = (xde_f-xde_0)/Tf;
lambda_y = (yde_f-yde_0)/Tf;
xd(1) = xde_0;
yd(1) = yde_0;
xd_(1)=lambda_x;
yd_(1)=lambda_y;

Nmax=Tf/dt + 1;
for k=2:Nmax;
   xd(k) = xd(k-1) + lambda_x*dt;
   yd(k) = yd(k-1) + lambda_y*dt;
   
   xd_(k)= lambda_x;
   yd_(k)= lambda_y;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2nd subtask: reference trajectory
qr1 = 45*pi/180; qr2 = -70*pi/180; qr3 = 0;
p2d = [qr1;qr2;qr3];
safe_pos=[qr1; qr2; qr3]; %%%***delete
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%% FIGURE FOR ANIMATION %%%%%%%%%%%%
%%%%%%% (moving robot stick diagram) %%%%%%%
fig0 = figure;

%%% GUInterface: control buttons to move obstacle
d_safe = 20;
h1 = uicontrol('Style', 'pushbutton', 'String', 'Up',...
   'Position', [10 150 45 40],...
   'Callback', 'box_y=box_y+0.1;'); %%deletee
h2 = uicontrol('Style', 'pushbutton', 'String', 'Down',...
   'Position', [10 110 45 40],...
   'Callback', 'box_y=box_y-0.1;'); %%delete

%%%%%%%%% START PLOTING FIGURE %%%%%%%%%%%%%%%
figure(fig0);
axis([-0.5, 2.5 -1.5 1.5])
axis on
hold on

xlabel('x (m)');
ylabel('y (m)');
plot(xd, yd, 'm:');
dtk=50; %% interval between consecutive plots <---***

plot([0], [0], 'o');

%% OBSTACLE
fill([0.5, 0.8, 0.8, 0.5], [-0.2, -0.2, -0.5, -0.5], [127/255 1 212/255]);

plot([0, xd2_0], [0, yd2_0],'r');
plot([xd2_0], [yd2_0], '*');
plot([xd2_0, xd3_0], [yd2_0, yd3_0],'m');
plot([xd3_0], [yd3_0], '*');
plot([xd3_0, xde_0], [yd3_0, yde_0],'r');
plot([xde_0], [yde_0], 'y*');
plot([xde_0], [yde_0], 'g+');

box_x=0.8;  %%current box position 
box_y=-0.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ****** KINEMATIC SIMULATION - Main loop ******
disp('Kinematic Simulation ...'); %%
disp(' '); %%


%% ***** INVESRE KINEMATICS  -->  DESIRED MOTION - JOINT SPACE *****
%% compute the reference joint-motion vectors:
%% {qd(k,i), i=1,...,n (num of degrees of freedom), with k=1,..., Nmax,}
%% and reference joint (angular) velocities {qd_1(k,i)}

%%% Inverse geometric model %%%%
%rd2 = xd(:).^2 + yd(:).^2;
%%qd(:,2) = acos( (rd2(:)-l(1)^2-l(2)^2)./(2*l(1)*l(2)) ); %% 1st solution: elbow down
%qd(:,2) = -acos( (rd2(:)-l(1)^2-l(2)^2)./(2*l(1)*l(2)) ); %% 2nd solution: elbow up
%ss2 = sin(qd(:,2));
%qd(:,1) = atan2(yd(:),xd(:)) - asin(l(2)*ss2(:)./sqrt(rd2(:)));

tt=0; tk=1;
qd(tk,1)=qd1_0; qd(tk,2)=qd2_0; qd(tk,3)=qd3_0;
while (tt<=Tf)
   
   %% Compute Jacobian
   s1 = sin(qd(tk,1)); 
   s12 = sin(qd(tk,1)+qd(tk,2)); 
   s123 = sin(qd(tk,1)+qd(tk,2)+qd(tk,3));
   c1 = cos(qd(tk,1)); 
   c12 = cos(qd(tk,1)+qd(tk,2)); 
   c123 = cos(qd(tk,1)+qd(tk,2)+qd(tk,3));
	Jac1(1,3) = -l(3)*s123;
   Jac1(1,2) = Jac1(1,3) - l(2)*s12;
   Jac1(1,1) = Jac1(1,2) - l(1)*s1;
	Jac1(2,3) = l(3)*c123;
   Jac1(2,2) = Jac1(2,3) + l(2)*c12;
   Jac1(2,1) = Jac1(2,2) + l(1)*c1;
   
   %%Pseudo-inverse computation of Jacobian matrix Jac1
   Jac1_psinv = Jac1'*inv(Jac1*Jac1');
   
   %Subtask 1
   task1=Jac1_psinv*[xd_(tk);yd_(tk)];
   
   %Subtask 1
   %%p2d=...
   p2=qd(tk,:)';
   H2=eye(3)*kc;
   task2 = (eye(3)-Jac1_psinv*Jac1)*H2*(p2d-p2);
   
   %angular velocity
   qd_(tk,:) = task1'+task2';
   

   %% ***** FORWARD KINEMATICS  JOINT MOTION -->  CARTESIAN POSITIONS *****
   %% ***** store successive positions for links 1, 2, and 3 (end-effector)
		%%(xd1, yd1) : cartesian position of the 1st link's local reference frame
		xd1(tk) = l(1)*c1;   
		yd1(tk) = l(1)*s1; 
		%%(xd2, yd2) : cartesian position of the 2nd link's local reference frame
		xd2(tk) = xd1(tk) + l(2)*c12; 
		yd2(tk) = yd1(tk) + l(2)*s12; 
		%%(xd3, yd3) : cartesian position of the 3rd link's local reference frame
		xd3(tk) = xd2(tk) + l(3)*c123; 
		yd3(tk) = yd2(tk) + l(3)*s123; 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %% numerical integration --> kinematic simulation
   qd(tk+1,1) = qd(tk,1) + dt*qd_(tk,1);
   qd(tk+1,2) = qd(tk,2) + dt*qd_(tk,2);
   qd(tk+1,3) = qd(tk,3) + dt*qd_(tk,3);
   
   
  	%% ***** PLOT MOTION - ANIMATE THE ROBOT SIMULATION *****
   %% (also, display moving box)
   if (mod(tk,dtk)==0),
      cla						%% clear plot at each time step (also possible to consider double-buffering?...)
      plot(xd, yd, 'm:');  %% replot desired trajectory (if given)
      
      fill([0.5, 0.8, 0.8, 0.5], [box_y, box_y, box_y+0.3, box_y+0.3], [0.2, 0.2, 0.7]);
      fill([0.5, 0.8, 0.8, 0.5], [box_y+1.3, box_y+1.3, box_y+1.3+0.3, box_y+1.3+0.3], [0.2, 0.2, 0.7]);

		plot([0], [0], 'o');
  		plot([0, xd1(tk)], [0, yd1(tk)],'r');
	  	plot([xd1(tk)], [yd1(tk)], '*');
  		plot([xd1(tk), xd2(tk)], [yd1(tk), yd2(tk)],'m');
	  	plot([xd2(tk)], [yd2(tk)], '*');
  		plot([xd2(tk), xd3(tk)], [yd2(tk), yd3(tk)],'b');
	  	plot([xd3(tk)], [yd3(tk)], 'y*');
     	plot([xd3(tk)], [yd3(tk)], 'g+');
       
      pause(0.4); %% pause motion to view successive robot configurations
   end
    


   
   tk=tk+1;	 %step increment, and
   tt=tt+dt; %time increment (for the simulation)
end




%% *** SAVE and PLOT output data ***
%%** use functions plot(...)
save;  %% --> save data to 'matlab.mat' file

fig1 = figure;

subplot(2,2,1);
plot(t,xd);
ylabel('xd (cm)');
xlabel('time t (sec)');


subplot(2,2,2);
plot(t,yd);
ylabel('yd (cm)');
xlabel('time t (sec)');


subplot(2,2,3);
plot(t,qd(:,1));
ylabel('qd1 (rad)');
xlabel('time t (sec)');


subplot(2,2,4);
plot(t,qd(:,2));
ylabel('qd2 (rad)');
xlabel('time t (sec)');






















