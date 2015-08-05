
clc;
clear all;
close all;

EXPORT = 1;
Ts = 0.01;      % sampling time

%%Physical parameters
m = 0.1;
Ixx = 0.0011;
Iyy = 0.0012;
Izz = 0.0013;
Jr = 0.0085;
l = 0.05;
g = 9.81;
%d = 7e-7;

DifferentialState phi theta psi p q rr x y z u v w;      % provide proper names for all differential states and controls etc.
Control u1 u2 u3 u4;

%% Differential Equations

ode = [ dot(phi) == p ,...
        dot(theta) == q , ...
        dot(psi) == rr, ...
        dot(p) == q*rr*((Iyy-Izz)/Ixx) + q*(Jr/Ixx)*u4 + (l/Ixx)*u2,...
        dot(q) == p*rr*((Izz-Ixx)/Iyy) - p*(Jr/Iyy)*u4 + (l/Iyy)*u3, ...
        dot(rr) == q*p*((Ixx-Iyy)/Izz) + (u4/Izz), ...
        dot(x) == u, ...
        dot(y) == v, ...
        dot(z) == w, ...
        dot(u) == ( (cos(phi)*sin(theta)*cos(psi))+(sin(phi)*sin(psi)) )*(u1/m), ...
        dot(v) == ( (cos(phi)*sin(theta)*sin(psi))-(sin(phi)*cos(psi)) )*(u1/m), ...
        dot(w) == -g + (cos(theta)*cos(phi))*(u1/m)];          % <-- TODO: define the set of differential equations
    

%% Export of a simulation routine:
acadoSet('problemname', 'sim');

sim = acado.SIMexport( Ts );

sim.setModel(ode);      % pass the ODE model

sim.set( 'INTEGRATOR_TYPE',             'INT_BDF'   );  % RK4 method
sim.set( 'NUM_INTEGRATOR_STEPS',        4           );

if EXPORT
    sim.exportCode( 'export_SIM' );
    
    cd export_SIM
    make_acado_integrator('../simulate_system')
    cd ..
end


%% Export of an optimization routine:
acadoSet('problemname', 'mpc');

N = 40;     % number of shooting intervals
ocp = acado.OCP( 0.0, N*Ts, N );

h = [phi; theta; psi; p; q; rr; x; y; z; u; v; w; u1; u2; u3; u4] ;   % <-- TODO: residual function for the stage cost
hN = [phi; theta; psi; p; q; rr; x; y; z; u ; v; w] ;   % <-- TODO: residual function for the stage cost
W = acado.BMatrix(eye(16));  % <-- TODO: weighting matrix for the stage cost
WN = acado.BMatrix(eye(12));

ocp.minimizeLSQ( W, h );            % stage cost
ocp.minimizeLSQEndTerm( WN, hN );   % terminal cost

vtransmin=-10;
vtransmax=10;

%constraints
%ocp.subjectTo( ode );

ocp.subjectTo( vtransmin <= v <= vtransmax );
ocp.subjectTo( -90.0 <= u1 <= 90.0 );
ocp.subjectTo( -90.0 <= u2 <= 90.0 );
ocp.subjectTo( -90.0 <= u3 <= 90.0 );
ocp.subjectTo( -90.0 <= u4 <= 90.0 );

%Acado settings
ocp.setModel(ode);      % pass the ODE model

mpc = acado.OCPexport( ocp );
mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'      );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
mpc.set( 'SPARSE_QP_SOLUTION',          'FULL_CONDENSING_N2');
mpc.set( 'INTEGRATOR_TYPE',             'INT_BDF'           );  % RK4 method
mpc.set( 'NUM_INTEGRATOR_STEPS',         20                 );
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
X0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; % <-- TODO: initial state (downward position)
Xref = [0, 0, 0.05, 0, 0, 0, 0, 0, 5, 0, 0, 0];     % <-- TODO: reference point (upward position)



Uref = [m*g, 0, 0, 0];
input.u = zeros(N,4);     % <-- TODO: initialization of the control trajectory
input.x = repmat(X0,N+1,1);      % <-- TODO: initialization of the state trajectory

input.y = [repmat(Xref,N,1) repmat(Uref,N,1)];   %  <-- TODO: reference trajectory for the stage cost
input.yN = Xref.';  % <-- TODO: reference trajectory for the terminal cost  

input.shifting.strategy = 0;    % shifting is optional but highly recommended with RTI!
                                %      1: use xEnd, 2: integrate

% SIMULATION LOOP
display('------------------------------------------------------------------')
display('               Simulation Loop'                                    )
display('------------------------------------------------------------------')

iter = 0; time = 0;
Tf = 10;
INFO_MPC = [];
controls_MPC = [];
state_sim = X0;
input.W = diag([10 10 1000 10 10 100 10 10 10 1 1 1 1 1 1 1]);
input.WN = eye(12);
input.x0 = X0.';
output = acado_MPCstep(input);

%MPC iteration
while time(end) < Tf
    tic
    % Solve NMPC OCP
    input.x0 = state_sim(end,:).';
    output = acado_MPCstep(input);
    
    % Save the MPC step
    INFO_MPC = [INFO_MPC; output.info];
    controls_MPC = [controls_MPC; output.u(1,:)];
    input.x = output.x;
    input.u = output.u;
    
    % Simulate system
    sim_input.x = state_sim(end,:).';
    sim_input.u = output.u(1,:).';
    states = simulate_system(sim_input);
    state_sim = [state_sim; states.value'];
    
    iter = iter+1;
    nextTime = iter*Ts; 
    disp(['current time: ' num2str(nextTime) '   ' char(9) ' (RTI step -- QP error: ' num2str(output.info.status) ',' ' ' char(2) ' KKT val: ' num2str(output.info.kktValue,'%1.2e') ',' ' ' char(2) ' CPU time: ' num2str(round(output.info.cpuTime*1e6)) ' µs)'])
    time = [time nextTime];
    
    %isualize(time, state_sim, Xref, xmin, xmax); 
end


%% 
figure;
subplot(2,4,1);
plot(time, state_sim(:,1)); hold on;
plot(time, state_sim(:,2)); hold on;
plot(time, state_sim(:,3)); hold on;

plot([0 time(end)], [0 0], 'r:');
xlabel('time(s)');
ylabel('Euler angles');

subplot(2,4,2);
plot(time, state_sim(:,4)); hold on;
plot(time, state_sim(:,5)); hold on;
plot(time, state_sim(:,6)); hold on;
plot([0 time(end)], [0 0], 'r:');
xlabel('time(s)');
ylabel('Euler rates');

subplot(2,4,3);
plot(time, state_sim(:,7)); hold on;
plot(time, state_sim(:,8)); hold on;
plot(time, state_sim(:,9)); hold on;
plot([0 time(end)], [0 0], 'r:');
xlabel('time(s)');
ylabel('Translational states');

subplot(2,4,4);
plot(time, state_sim(:,10)); hold on;
plot(time, state_sim(:,11)); hold on;
plot(time, state_sim(:,12)); hold on;
plot([0 time(end)], [0 0], 'r:');
xlabel('time(s)');
ylabel('Translational velocities');


subplot(2,4,[5 7]);
stairs(time(1:end-1), controls_MPC(:,1),'r'); hold on;
stairs(time(1:end-1), controls_MPC(:,2),'b'); hold on;
stairs(time(1:end-1), controls_MPC(:,3),'g'); hold on;
stairs(time(1:end-1), controls_MPC(:,4)); hold on;
plot([0 time(end)], [0 0], 'r:');
% plot([0 time(end)], [Fmin Fmin], 'g--');
% plot([0 time(end)], [Fmax Fmax], 'g--');
xlabel('time(s)');
ylabel('Inputs');




