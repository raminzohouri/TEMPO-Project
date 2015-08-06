A = eye(9);
% Generate the geometry used to draw the quadcopter
r = .5; d = 1.25; h = .25; %inches: rotor dia., quad motor distance from 
% center of mass, and rotor height above arms (entirely cosmetic)
a = 1; b = 1; c = 0.2; % Top , bottom, and sides lengths of hub (?)

% Construct rotor representations
N = [d  0 h].';% m1 rotor center [X Y Z]
E = [0 -d h].';% m4 rotor center
W = [0  d h].';% m2 rotor center
S = [-d 0 h].';% m3 rotor center
Nr = circlePoints(N, r, 10); Nr = [Nr Nr(:,1)]; % Rotor blade circles
Er = circlePoints(E, r, 10); Er = [Er Er(:,1)];
Wr = circlePoints(W, r, 10); Wr = [Wr Wr(:,1)];
Sr = circlePoints(S, r, 10); Sr = [Sr Sr(:,1)];
% Motors connecting to center of blade circles
mN = [d,d;
      0,0;
      h,0];
mE = [0,0;
     -d,-d;
      h,0];
mW = [0,0;
      d,d;
      h,0];
mS = [-d,-d;
       0,0;
       h,0];
% Construct body plot points
bNS = [ d, -d;
        0,  0;
        0,  0]; %For drawing the body "X" shape
bEW = [ 0,  0;
        d, -d;
        0,  0];
% Body (HUB) Squares
Top = [ a/2,   0,-a/2,   0;
          0, b/2,   0,-b/2;
        c/2, c/2, c/2, c/2];
Bot = vertcat(Top(1:2,:),-Top(3,:)); % Bot is same as top just below the center of mass
NEB = [ a/2, a/2,   0,   0;
          0,   0, b/2, b/2;
        c/2,-c/2,-c/2, c/2]; % By the way... north east is actually north west from above since x is north and y is east :P
NWB = [ a/2, a/2,   0,   0;
          0,   0,-b/2,-b/2;
        c/2,-c/2,-c/2, c/2];
SEB = -NWB;
SWB = -NEB;

phi = A(1,4);
the = A(1,5);
psi = A(1,6);

% ROTATION MATRIX --- ZYX ROTATION (R = Rib)
R = [cos(psi)*cos(the) cos(psi)*sin(the)*sin(phi)-sin(psi)*cos(phi) cos(psi)*sin(the)*cos(phi)+sin(psi)*sin(phi);
       sin(psi)*cos(the) sin(psi)*sin(the)*sin(phi)+cos(psi)*cos(phi) sin(psi)*sin(the)*cos(phi)-cos(psi)*sin(phi);
       -sin(the)         cos(the)*sin(phi)                            cos(the)*cos(phi)];

% Rotate body frame velocity vector
U = A(:,7);
V = A(:,8);
W = A(:,9);
Vi = zeros(length(A),3);
MvMax= max(sqrt(U.^2+V.^2+W.^2)); 
Vb = 3/MvMax*[U(1), V(1), W(1)]'; % Scale velocity
Vi(1,:) = R*Vb; % Rotate velocity vector to inertial frame for plotting

% Support for X-configuration nifty trick (Requires that quadModel
% structure is in base workspace)
% if (~evalin('base','quadModel.plusConfig'))
%     Rz = [ sqrt(2)/2, sqrt(2)/2, 0;
%           -sqrt(2)/2,sqrt(2)/2, 0;
%                    0,          0, 1];
%     Nr = Rz*Nr;
%     Er = Rz*Er;
%     Wr = Rz*Wr;
%     Sr = Rz*Sr;
%     mN = Rz*mN;
%     mE = Rz*mE;
%     mW = Rz*mW;
%     mS = Rz*mS;
%     bNS = Rz*bNS;
%     bEW = Rz*bEW;
%     Top = Rz*Top;
%     Bot = Rz*Bot;
%     NEB = Rz*NEB;
%     NWB = Rz*NWB;
%     SWB = Rz*SWB;
%     SEB = Rz*SEB;
% end

% Rotate body parts Via Initialized R
NrR = R*Nr;
ErR = R*Er;
WrR = R*Wr;
SrR = R*Sr;
mNr = R*mN;
mEr = R*mE;
mWr = R*mW;
mSr = R*mS;
bNSR = R*bNS;
bEWR = R*bEW;
TopR = R*Top;
BotR = R*Bot;
NEBR = R*NEB;
NWBR = R*NWB;
SWBR = R*SWB;
SEBR = R*SEB;

% Plot the box rotation and ang. velocity and inertial frame velocity
% vector
figure();
plot3(bNSR(1,:),bNSR(2,:),bNSR(3,:),'b','LineWidth',3) % Body Arm
hold on
plot3(bEWR(1,:),bEWR(2,:),bEWR(3,:),'b','LineWidth',3) % Body Arm
plot3(mNr(1,:),mNr(2,:),mNr(3,:),'k','LineWidth',4) % Motor
plot3(mEr(1,:),mEr(2,:),mEr(3,:),'k','LineWidth',4) % Motor
plot3(mWr(1,:),mWr(2,:),mWr(3,:),'k','LineWidth',4) % Motor
plot3(mSr(1,:),mSr(2,:),mSr(3,:),'k','LineWidth',4) % Motor
plot3(NrR(1,:),NrR(2,:),NrR(3,:),'g') % Motor 1 blades
plot3(ErR(1,:),ErR(2,:),ErR(3,:),'k') % Motor 4 blades
plot3(WrR(1,:),WrR(2,:),WrR(3,:),'k') % Motor 2 blades
plot3(SrR(1,:),SrR(2,:),SrR(3,:),'g') % Motor 3 blades
grey = [0.5 0.5 0.5];
top = fill3(TopR(1,:),TopR(2,:),TopR(3,:),'r'); alpha(top,0.8); % Top Surface
bot = fill3(BotR(1,:),BotR(2,:),BotR(3,:),'g'); alpha(bot,0.8); % Bottom Surface
ne  = fill3(NEBR(1,:),NEBR(2,:),NEBR(3,:),'c'); alpha(ne,0.8); % North East surface
nw  = fill3(NWBR(1,:),NWBR(2,:),NWBR(3,:),grey); alpha(nw,0.8); % North West surface
sw  = fill3(SWBR(1,:),SWBR(2,:),SWBR(3,:),grey); alpha(sw,0.8); % South West surface
se  = fill3(SEBR(1,:),SEBR(2,:),SEBR(3,:),grey); alpha(se,0.8); % South East surface
axis square
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([-3 3])
ylim([-3 3])
zlim([-3 3])
% view(handles.AZval,handles.ELval)
grid on

omi = zeros(length(A),3); % Initialize omega inertiaF points array (L(A)x3)
P = A(:,1); Q = A(:,2); Rw = A(:,3);
MombMax = max(sqrt(P.^2+Q.^2+Rw.^2)); % Calculate max magnitude of omb
omb = 3/MombMax*[P,Q,Rw].'; % Store and scale current omega bodyF
omi(1,:) = R*omb(:,1); % Rotate omegab to inertiaF store in omegai array
qp1 = quiver3(0,0,0,omi(1,1),omi(1,2),omi(1,3),'ro');
qp2 = quiver3(0,0,0,Vi(1,1),Vi(1,2),Vi(1,3),'k');
hold off