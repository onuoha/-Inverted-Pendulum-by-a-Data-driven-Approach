%% Closed-loop Subspace Algorithmus 1/2
% Dieses Skript erstellt:
% - das Modell des inversen Pendels
% - den stabilisierenden Regler
% - Definitionen fuer den closed-loop Subspace Algorithmus
close all; clear
%% Parameter des Pendels
m_p = 0.329; % Masse des Pendelstabes                       [kg]
m_s = 3.2;   % Masse des Schlitten                          [kg]
l_p = 0.44;  % Abstand Pendelschwerpunkt zur Aufhaengung    [m]
R_s = 6.2;   % Reibungskonstante fuer den Schlitten         [kg/s]
R_p = 0.009; % Reibungskonstante fuer den Pendelstab        [kgm^2/s]
g   = 9.81;  % Erdbeschleunigung                            [m/s^2]
J_p = 0.072; % Auf die Drehachse bezogenes Traegheitsmoment [kgm^2]

%% Zustandsraummodell Parameter
% x1 := Positon des Wagens (Position of the cart)
% x2 := Geschwindigkeit des Wagens (velocity of the cart)
% x3 := Winkel des Pendels (Angle of the pendulum)
% x4 := Winkelgeschwindigkeit des Pendels (Angular velocity of the pendulum)
% u  := Am Pendel angreifende Kraft (Force acting on the pendulum)

% Systemmatrizen (kontinuierlich)
a_22 = -(J_p*R_s)/(J_p*(m_s+m_p)-m_p^2*l_p^2);
a_23 = -(m_p^2*l_p^2*g)/(J_p*(m_s+m_p)-m_p^2*l_p^2);
a_24 = (m_p*l_p*R_p)/(J_p*(m_s+m_p)-m_p^2*l_p^2);
a_42 = (m_p*l_p*R_s)/(J_p*(m_s+m_p)-m_p^2*l_p^2);
a_43 = ((m_p+m_s)*m_p*l_p*g)/(J_p*(m_s+m_p)-m_p^2*l_p^2);
a_44 = -((m_p+m_s)*R_p)/(J_p*(m_s+m_p)-m_p^2*l_p^2);
A = [0    1       0       0;
     0    a_22    a_23    a_24;
     0    0       0       1
     0    a_42    a_43    a_44];
B = [0;
    J_p/(J_p*(m_s+m_p)-m_p^2*l_p^2);
    0;
    -m_p*l_p/(J_p*(m_s+m_p)-m_p^2*l_p^2)];
C = [1   0   0   0;
     0   1   0   0;
     0   0   1   0];

D = [0; 0; 0];
K = 0.02*(2*rand(4,4)-1);                       
h = 0.02;                               % Abtastrate
sys_ss = ss(A,B,C,D);          % Kontinuierliches Zustandsraummodell
sys_ss_d = c2d(sys_ss, h);            % Diskretes Zustandsraummodell
[A_d,B_d,C_d,D_d] = ssdata(sys_ss_d);   % Diskrete Systemmatrizen

%% Inverses Pendel Definitionen
m = 1; % Number of Inputs
l = 4; % Number of Outputs
n = 4; % Number of States (nicht notwendig)

%% Closed-loop Subspace Algorithmus Definitionen
s = 1000;       % Anzahl Daten
M = 70;         % backward Horizon, Anzahl Reihen (past)
N = 20;         % forward Horizon, Anzahl Reihen (future)
j = s-M-N+1;    % Anzahl Spalten der Subspace Matrizen
r_f = 0*ones((l+m-1)*N,1); r_k = 0; % Referenzsignal

%% Erstelle stabilisierenden Regler
Q = eye(n);
Q(1,1) = 200;       % Wagenposition
Q(3,3) = 2000;      % Pendelwinkel
R = 0.2*eye(m);     % Eingangsgewichtung
K_d = dlqr(A_d,B_d,Q,R);
H_dd = [0 0 0 0]';
D_dd = [0 0 0 0]';
syse = ss(A_d,[B_d randn(n,1)],eye(4),[D_dd H_dd],h);
[kest,L,P] = kalman(syse,1,1e-5*eye(l));
rlqg = lqgreg(kest,K_d);
[A_c,B_c,C_c,D_c] = ssdata(rlqg);
C_c = -C_c; D_c = -D_c;
