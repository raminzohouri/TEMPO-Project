
clc;
clear all;
close all;

EXPORT = 0;
Ts = 0.01;      % sampling time
Nx = 12;
Nu = 4;

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

sim.set( 'INTEGRATOR_TYPE',             'INT_IRK_GL4'   );  % RK4 method
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
W = acado.BMatrix(eye(Nx+Nu));  % <-- TODO: weighting matrix for the stage cost
WN = acado.BMatrix(eye(Nx));

ocp.minimizeLSQ( W, h );            % stage cost
ocp.minimizeLSQEndTerm( WN, hN );   % terminal cost



% %constraints to tranlational velocity
% vtransmin=-7.5; vtransmax=7.5;
% ocp.subjectTo( vtransmin <= v <= vtransmax );
% ocp.subjectTo( vtransmin <= u <= vtransmax );
% ocp.subjectTo( vtransmin <= w <= vtransmax );
% 
% euler_rates_min=-2*pi; euler_rates_max=2*pi;
% ocp.subjectTo( euler_rates_min <= p <= euler_rates_max );
% ocp.subjectTo( euler_rates_min <= q <= euler_rates_max );
% ocp.subjectTo( euler_rates_min <= rr <= euler_rates_max );
% 
% u1_min = -10; u1_max = 10; 
% u2_min = -10; u2_max = 10; 
% u3_min = -10; u3_max = 10; 
% u4_min = -10; u4_max = 10; 
% 
% ocp.subjectTo( u1_min <= u1 <= u1_max );
% ocp.subjectTo( u2_min <= u2 <= u2_max );
% ocp.subjectTo( u3_min <= u3 <= u3_max );
% ocp.subjectTo( u4_min <= u4 <= u4_max );

%Acado settings
ocp.setModel(ode);      % pass the ODE model

mpc = acado.OCPexport( ocp );
mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'      );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
mpc.set( 'SPARSE_QP_SOLUTION',          'FULL_CONDENSING_N2');
mpc.set( 'INTEGRATOR_TYPE',             'INT_IRK_GL4'           );  % RK4 method
mpc.set( 'NUM_INTEGRATOR_STEPS',         2*N                 );
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
Xref = [0, 0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0];     % <-- TODO: reference point (upward position)



Uref = [m*g, 0, 0, 0];
input.u = zeros(N,4);     % <-- TODO: initialization of the control trajectory
input.x = repmat(X0,N+1,1);      % <-- TODO: initialization of the state trajectory

input.y = [repmat(Xref,N,1) repmat(Uref,N,1)];   %  <-- TODO: reference trajectory for the stage cost
input.yN = Xref.';  % <-- TODO: reference trajectory for the terminal cost  

input.shifting.strategy = 1;    % shifting is optional but highly recommended with RTI!
                                %      1: use xEnd, 2: integrate

% SIMULATION LOOP
display('------------------------------------------------------------------')
display('               Simulation Loop'                                    )
display('------------------------------------------------------------------')

iter = 0; time = 0;
Tf = 15;
INFO_MPC = [];
controls_MPC = [];
state_sim = X0;
input.W = diag([1e-2 1e-2 1e-2 1e-3 1e-3  1e-3 1000 1000 1000 1e-2 1e-2 1e-2 1e-5 1e-5 1e-5 1]);
input.WN = diag([1 1 1 1e-2 1e-2  1e-2 1000 1000 1000 1e-2 1e-2 1e-2]);
input.x0 = X0.';
output = acado_MPCstep(input);
ref_traj =input.y;

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

    % shift reference:
    ref_traj = [ref_traj; [input.yN.' Uref]];
    input.y = [input.y(2:end,:); [input.yN.' Uref]];
%     input.yN(end-3) = sin(time(end));
%     input.yN(end-4) = sin(time(end));
    if (time(end) > 1)
    input.yN(end-5) = 1;
    end
    if (time(end) > 5)
    input.yN(end-4) = 1;
    end
    if (time(end) > 10)
    input.yN(end-3) = 1;
    end
    
    % Simulate system
    sim_input.x = state_sim(end,:).';
    sim_input.u = output.u(1,:).';
    states = simulate_system(sim_input);
    state_sim = [state_sim; states.value'];
    
    iter = iter+1;
    nextTime = iter*Ts; 
    disp(['current time: ' num2str(nextTime) '   ' char(9) ' (RTI step -- QP error: ' num2str(output.info.status) ',' ' ' char(2) ' KKT val: ' num2str(output.info.kktValue,'%1.2e') ',' ' ' char(2) ' CPU time: ' num2str(round(output.info.cpuTime*1e6)) ' µs)'])
    time = [time nextTime];
    
end


%% 
figure;
subplot(2,4,1);
plot(time, state_sim(:,1),'r'); hold on;
plot(time, state_sim(:,2),'g'); hold on;
plot(time, state_sim(:,3),'b'); hold on;

plot([0 time(end)], [0 0], 'r:');
xlabel('time(s)');
ylabel('Euler angles');

subplot(2,4,2);
plot(time, state_sim(:,4),'r'); hold on;
plot(time, state_sim(:,5),'g'); hold on;
plot(time, state_sim(:,6),'b'); hold on;
plot([0 time(end)], [0 0], 'r:');
xlabel('time(s)');
ylabel('Euler rates');

subplot(2,4,3);
plot(time, state_sim(:,7),'r'); hold on;
plot(time, state_sim(:,8),'g'); hold on;
plot(time, state_sim(:,9),'b'); hold on;
plot([0 time(end)], [0 0], 'r:');
xlabel('time(s)');
ylabel('Translational states');

subplot(2,4,4);
plot(time, state_sim(:,10),'r'); hold on;
plot(time, state_sim(:,11),'g'); hold on;
plot(time, state_sim(:,12),'b'); hold on;
plot([0 time(end)], [0 0], 'r:');
xlabel('time(s)');
ylabel('Translational velocities');


subplot(2,4,5);
stairs(time(1:end-1), controls_MPC(:,1),'r'); hold on;
subplot(2,4,6);
stairs(time(1:end-1), controls_MPC(:,2),'b'); hold on;
subplot(2,4,7);
stairs(time(1:end-1), controls_MPC(:,3),'g'); hold on;
subplot(2,4,8);
stairs(time(1:end-1), controls_MPC(:,4)); hold on;
plot([0 time(end)], [0 0], 'r:');
xlabel('time(s)');
ylabel('Inputs');


figure(2)
subplot(3,1,1);
plot(time, state_sim(:,7),'b'); hold on;
plot(time, ref_traj(N:end,7),'r--');

subplot(3,1,2);
plot(time, state_sim(:,8),'b'); hold on;
plot(time, ref_traj(N:end,8),'r--');

subplot(3,1,3);
plot(time, state_sim(:,9),'b'); hold on;
plot(time, ref_traj(N:end,9),'r--');





