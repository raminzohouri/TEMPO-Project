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

