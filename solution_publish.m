%% Solution for Exercise 10: ACADO code generation for Nonlinear MPC
%    TEMPO Summer School on Numerical Optimal Control and Embedded Optimization 
%    University of Freiburg, July 27 - August 7, 2015 
%    Rien Quirynen, Dimitris Kouzoupis and Moritz Diehl 

%% Swing-up of an inverted pendulum
clc;
clear all;
close all;

Ts = 0.05;
EXPORT = 1;

DifferentialState x1 theta v1 dtheta;
Control F;

M = 1;
m = 0.1;
g = 9.81;
l = 0.8;

% Differential Equation
    
f_expl = [   dot(x1) == v1; ...
             dot(theta) == dtheta; ...
             dot(v1) == (- l*m*sin(theta)*dtheta^2 + F + g*m*cos(theta)*sin(theta))/(M + m - m*cos(theta)^2); ...
             dot(dtheta) == (- l*m*cos(theta)*sin(theta)*dtheta^2 + F*cos(theta) + g*m*sin(theta) + M*g*sin(theta))/(l*(M + m - m*cos(theta)^2)) ];

% SIMexport
acadoSet('problemname', 'sim');

numSteps = 2;
sim = acado.SIMexport( Ts );
sim.setModel(f_expl);
sim.set( 'INTEGRATOR_TYPE',             'INT_IRK_GL2' );
sim.set( 'NUM_INTEGRATOR_STEPS',        numSteps        );

if EXPORT
    sim.exportCode( 'export_SIM' );
    
    cd export_SIM
    make_acado_integrator('../integrate_pendulum')
    cd ..
end

% MPCexport
acadoSet('problemname', 'mpc');

N = 40;
ocp = acado.OCP( 0.0, N*Ts, N );

h = [x1 theta v1 dtheta F];
hN = [x1 theta v1 dtheta];

W = acado.BMatrix(eye(length(h)));
WN = acado.BMatrix(eye(length(hN)));

ocp.minimizeLSQ( W, h );
ocp.minimizeLSQEndTerm( WN, hN );

xmin = -2; xmax = 2;
Fmin = -20; Fmax = 20;

ocp.subjectTo( xmin <= x1 <= xmax );
ocp.subjectTo( Fmin <= F <= Fmax );
% ocp.subjectTo( 'AT_END', [x1 theta v1 dtheta] == 0 );

ocp.setModel(f_expl);

mpc = acado.OCPexport( ocp );
mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'      );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
mpc.set( 'SPARSE_QP_SOLUTION',          'FULL_CONDENSING_N2');
mpc.set( 'LEVENBERG_MARQUARDT',          1e-5               );
mpc.set( 'INTEGRATOR_TYPE',             'INT_RK4'           );
mpc.set( 'NUM_INTEGRATOR_STEPS',        2*N                 );
mpc.set( 'QP_SOLVER',                   'QP_QPOASES'    	);

if EXPORT
    mpc.exportCode( 'export_MPC' );
    
    global ACADO_;
    copyfile([ACADO_.pwd '/../../external_packages/qpoases'], 'export_MPC/qpoases')
    
    cd export_MPC
    make_acado_solver('../acado_MPCstep')
    cd ..
end

% PARAMETERS SIMULATION
X0 = [0 pi 0 0];
Xref = [0 0 0 0];
input.x = [repmat(X0,N/2,1); repmat(Xref,N/2+1,1)];
input.od = [];

Uref = zeros(N,1);
input.u = Uref;

input.y = [repmat(Xref,N,1) Uref];
input.yN = Xref.';

input.W = diag([5e-1 1 2e-3 2e-3 1e-4]);
input.WN = diag([5e-1 1 2e-3 2e-3]);

input.shifting.strategy = 1;

% SIMULATION LOOP
display('------------------------------------------------------------------')
display('               Simulation Loop'                                    )
display('------------------------------------------------------------------')

iter = 0; time = 0;
Tf = 4;
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
    sim_input.u = output.u(1,:).';
    sim_input.od = 0.2;
    states = integrate_pendulum(sim_input);
    state_sim = [state_sim; states.value'];
    
    iter = iter+1;
    nextTime = iter*Ts; 
    disp(['current time: ' num2str(nextTime) '   ' char(9) ' (RTI step -- QP error: ' num2str(output.info.status) ',' ' ' char(2) ' KKT val: ' num2str(output.info.kktValue,'%1.2e') ',' ' ' char(2) ' CPU time: ' num2str(round(output.info.cpuTime*1e6)) ' µs)'])
    time = [time nextTime];
    
    visualize(time, state_sim, Xref, xmin, xmax, l);
end

figure;
subplot(2,2,1);
plot(time, state_sim(:,1)); hold on;
plot([0 time(end)], [0 0], 'r:');
plot([0 time(end)], [xmin xmin], 'g--');
plot([0 time(end)], [xmax xmax], 'g--');
xlabel('time(s)');
ylabel('x');
ylim([1.5*xmin 1.5*xmax])

subplot(2,2,2);
plot(time, state_sim(:,2)); hold on;
plot([0 time(end)], [0 0], 'r:');
xlabel('time(s)');
ylabel('theta');

subplot(2,2,[3 4]);
stairs(time(1:end-1), controls_MPC); hold on;
plot([0 time(end)], [0 0], 'r:');
plot([0 time(end)], [Fmin Fmin], 'g--');
plot([0 time(end)], [Fmax Fmax], 'g--');
xlabel('time(s)');
ylabel('F');
ylim([1.5*Fmin 1.5*Fmax])


figure;
semilogy(time(1:end-1), KKT_MPC, ':bx');
xlabel('time(s)')
ylabel('KKT')


%% Generalized tangential predictor of RTI for Nonlinear MPC
clc;
clear all;
close all;

% Define parameters:
N = 40;
Xref = [0 0 0 0];
X0 = [1.5 -0.2 0 0].';
input.x = repmat(X0.',N+1,1);
input.x0 = X0;

Uref = zeros(N,1);
input.u = Uref;

input.y = [repmat(Xref,N,1) Uref];
input.yN = Xref.';

input.W = diag([5e-1 1 5e-2 5e-2 2e-2]);
input.WN = diag([5e-1 1 5e-2 5e-2]);

% Converge the SQP method as the linearization point:
for k = 1:30
    output = acado_MPCstep(input);
    input.x = output.x;
    input.u = output.u;
end
X_RTI = output.x; U_RTI = output.u;

cost_RTI = []; u_RTI = [];
cost_SQP = []; u_SQP = [];
p_values = 1.5:0.001:1.8;
for p = p_values
    X0 = [p -0.2 0 0].';
    input.x0 = X0;
    
    % Perform one RTI step:
    input.x = X_RTI; input.u = U_RTI; % NOTE: we use the same, original linearization point
    output = acado_MPCstep(input);
    cost_RTI = [cost_RTI; output.info.objValue];
    u_RTI = [u_RTI; output.u(1)];
    if output.info.status ~= 0
       error(['output status = ' num2str(output.info.status) ' for p = ' num2str(p)]) 
    end
    
    % Perform a fixed amount of 30 SQP steps:
    for k = 1:30
        output = acado_MPCstep(input);
        input.x = output.x;
        input.u = output.u;
    end
    cost_SQP = [cost_SQP; output.info.objValue];
    u_SQP = [u_SQP; output.u(1)];
    if output.info.status ~= 0
       error(['output status = ' num2str(output.info.status) ' for p = ' num2str(p)]) 
    end
end

figure;
plot(p_values, u_RTI, ':rx'); hold on;
plot(p_values, u_SQP, '-b');
plot(p_values(1), u_SQP(1), 'ko', 'MarkerSize', 10)
xlabel('position');
ylabel('feedback control')
legend('RTI step', 'exact', 'linearization point');
xlim([min(p_values)-0.01 max(p_values)+0.01])
ylim([min(u_SQP)-1 max(u_SQP)+1]);
title('Illustration of the generalized tangential predictor')


