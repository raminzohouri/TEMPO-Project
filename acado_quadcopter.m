clc;
clear all;
close all;

Ts = 0.01;
EXPORT = 1;

DifferentialState  phi theta psi p q rr u v w;
Control u1 u2 u3 u4;

%Physical parameters
m = 2.1;
Ixx = 0.011;
Iyy = 0.012;
Izz = 0.013;
Jr = 0.085;
l = 0.28;
g = 9.81;
d = 2e-7;
%% Differential Equation
f_expl = acado.DifferentialEquation();
f_expl.add(dot(phi) == p);
f_expl.add(dot(theta) == q);
f_expl.add(dot(psi) == rr);
f_expl.add(dot(p) == q*rr*((Iyy-Izz)/Ixx) + q*(Jr/Ixx)*u1 + (l/Ixx)*u2);
f_expl.add(dot(q) == p*rr*((Izz-Ixx)/Iyy) - p*(Jr/Iyy)*u1 + (l/Iyy)*u3);
f_expl.add(dot(rr) == q*p*((Ixx-Iyy)/Izz) + (d/Izz)*u4);
f_expl.add(dot(u) == ( (cos(phi)*sin(theta)*cos(psi))+(sin(phi)*sin(psi)) )*(u1/m));
f_expl.add(dot(v) == ( (cos(phi)*sin(theta)*cos(psi))-(sin(phi)*cos(psi)) )*(u1/m));
f_expl.add(dot(w) == g - (cos(phi)*cos(theta))*(u1/m));
%% SIMexport

acadoSet('problemname', 'sim');
%controllaw = acado.RealTimeAlgorithm;
%controller = acado.Controller( controllaw [, reference] );
numSteps = 6;
sim = acado.SIMexport( Ts );
sim.setModel(f_expl);
sim.set( 'INTEGRATOR_TYPE',             'INT_RK4' );
sim.set( 'NUM_INTEGRATOR_STEPS',        numSteps        );

if EXPORT
    sim.exportCode( 'export_SIM' );
    
    cd export_SIM
    make_acado_integrator('../integrate_pendulum')
    cd ..
end

%% MPCexport
acadoSet('problemname', 'mpc');

N = 40;

      
r = [ 0; 0; 0; 0; 0; 0; 0; 0; 0];
tStart = 0.0;
tEnd   = N*Ts;

ocp = acado.OCP(tStart, tEnd, N );

h =[ phi; theta; psi; p; q; rr; u; v; w; u1; u2; u3; u4]; 
hN = [ phi; theta; psi; p; q; rr; u; v; w];

W = (eye(length(h)))*0.01;
WN = (eye(length(hN)))*0.0001;

ocp.minimizeLSQ(W,h);

ocp.minimizeLSQEndTerm( WN, hN);
ocp.subjectTo( -90.0 <= u1 <= 90.0 );
ocp.subjectTo( -90.0 <= u2 <= 90.0 );
ocp.subjectTo( -90.0 <= u3 <= 90.0 );
ocp.subjectTo( -90.0 <= u4 <= 90.0 );
 

xmin = -2; xmax = 2;
%Fmin = -20; Fmax = 20;

%ocp.subjectTo( xmin <= x1 <= xmax );
%ocp.subjectTo( Fmin <= F <= Fmax );
% ocp.subjectTo( 'AT_END', [x1 theta v1 dtheta] == 0 );

ocp.setModel(f_expl);

mpc = acado.OCPexport( ocp );
mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'      );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
mpc.set( 'SPARSE_QP_SOLUTION',          'FULL_CONDENSING_N2');
mpc.set( 'INTEGRATOR_TYPE',             'INT_RK4'           );  % RK4 method
mpc.set( 'NUM_INTEGRATOR_STEPS',        20                 );
mpc.set( 'QP_SOLVER',                   'QP_QPOASES'    	);
if EXPORT
    mpc.exportCode( 'export_MPC' );
    
   global ACADO_;
    copyfile([ACADO_.pwd '/../../external_packages/qpoases'], 'export_MPC/qpoases')
    
    cd export_MPC
    make_acado_solver('../acado_MPCstep')
    cd ..
end

%% PARAMETERS SIMULATION
X0 = zeros(1,9)+1;
Xref = zeros(1,9);
simulationStart =  0.0;
simulationEnd   =  16.0;
input.x = repmat(Xref,(N+1),1);
%input.od = [];

Uref = [1,0,0,0];
input.u = zeros(N,4);

input.y = [repmat(Xref,N,1) repmat(Uref,N,1)];
input.yN = Xref.';

input.W = (eye(length(h)));
input.WN =(eye(length(hN)));

input.shifting.strategy =1;

%% SIMULATION LOOP
display('------------------------------------------------------------------')
display('               Simulation Loop'                                    )
display('------------------------------------------------------------------')

iter = 0; time = 0;
Tf = 10;
KKT_MPC = []; INFO_MPC = [];
controls_MPC = [];
state_sim = X0;


while time(end) < Tf
    tic
    % Solve NMPC OCP
    input.x0 = state_sim(end,:).';
    output = acado_MPCstep(input);
    
    % Save the MPC step
    INFO_MPC = [INFO_MPC; output.info];
    KKT_MPC = [KKT_MPC; output.info.kktValue];
    controls_MPC = [controls_MPC; output.u(1,:)];

    input.x = output.x;
    input.u = output.u;
    
    % Simulate system
    sim_input.x = state_sim(end,:).';
    sim_input.u = 0*output.u(1,:).';
    %sim_input.od = 0.2;
    states = integrate_pendulum(sim_input);
    state_sim = [state_sim; states.value'];
    
    iter = iter+1;
    nextTime = iter*Ts; 
    disp(['current time: ' num2str(nextTime) '   ' char(9) ' (RTI step -- QP error: ' num2str(output.info.status) ',' ' ' char(2) ' KKT val: ' num2str(output.info.kktValue,'%1.2e') ',' ' ' char(2) ' CPU time: ' num2str(round(output.info.cpuTime*1e6)) ' µs)'])
    time = [time nextTime];
    
    %visualize(time, state_sim, Xref, xmin, xmax, l); 
    pause(0.75*abs(Ts-toc));
end

figure;
subplot(2,3,1);
title('Controls')
plot(time, state_sim(:,1)); hold on;
plot(time, state_sim(:,2)); hold on;
plot(time, state_sim(:,3)); hold on;

plot([0 time(end)], [0 0], 'r:');
plot([0 time(end)], [xmin xmin], 'g--');
plot([0 time(end)], [xmax xmax], 'g--');
xlabel('time(s)');
ylabel('Euler Angles');

subplot(2,3,2);
plot(time, state_sim(:,4)); hold on;
plot(time, state_sim(:,5)); hold on;
plot(time, state_sim(:,6)); hold on;
plot([0 time(end)], [0 0], 'r:');
xlabel('time(s)');
ylabel('Angular Velocity');

subplot(2,3,3);
plot(time, state_sim(:,7)); hold on;
plot(time, state_sim(:,8)); hold on;
plot(time, state_sim(:,9)); hold on;
plot([0 time(end)], [0 0], 'r:');
xlabel('time(s)');
ylabel('Translational Velocity');


subplot(2,3,[4 6]);
stairs(time(1:end-1), controls_MPC(:,1),'r'); hold on;
stairs(time(1:end-1), controls_MPC(:,2),'b'); hold on;
stairs(time(1:end-1), controls_MPC(:,3),'g'); hold on;
stairs(time(1:end-1), controls_MPC(:,4)); hold on;
plot([0 time(end)], [0 0], 'r:');
% plot([0 time(end)], [Fmin Fmin], 'g--');
% plot([0 time(end)], [Fmax Fmax], 'g--');
xlabel('time(s)');
ylabel('Controls');
