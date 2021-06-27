function out=PID_SSE(pid, it, maxiter)

% system model parameters - CSTR: continously stirred tank reactor - water level control
Da1=3; Da2=0.5; Da3=1;  
Ts=0.1; % sampling time
t=120;  % duration simulation
N=t/Ts; 
d=3;    % system state vector dimension
x0=[0.1 0.1 0.1]'; % initial states

Kp=pid(1); Ki=pid(2); Kd=pid(3);
        
dot_yref_yaw=zeros(1, N);   % derivative of the output reference 
y_yaw=zeros(1, N);          % measured yaw output - measured by IMU
dot_y_yaw=zeros(1, N);      % measured deriveative of yaw output which is yaw velocity
e_yaw=zeros(1, N);          % tracking error
dot_e_yaw=zeros(1, N);      % derivative of tracking error
intgrl_e_yaw=zeros(1, N);   % integral/summation of tracking error
u_yaw=zeros(1, N);              % control input signal to the plant

SNR=40; % signal to noise ratio - dB
d2_fixed=1; noisy_case=0; sinusiodal_case=0;
yref_yaw=zeros(1, N);   % reference output - which is desired to be tracked
n=1:N;
if sinusiodal_case==1
    T=25; % 25
    yref_yaw=0.5+0.1*sin(2*pi*(1/T)*Ts.*n);
else
    block=t/8; % 50 s if t=300 s   -linspace(-0.15, -0.1, rft);
    rft=2/Ts; % rise-fall time - 1 s
    yref_yaw(block*0+1:block*1/Ts-rft)=0.5*ones(1, block/Ts-rft);
    yref_yaw(block*1/Ts-rft+1:block*1/Ts)=linspace(0.5, 0.55, rft);
    yref_yaw(block*1/Ts+1:block*2/Ts-rft)=0.55*ones(1, block/Ts-rft);
    yref_yaw(block*2/Ts-rft+1:block*2/Ts)=-linspace(-0.55, -0.45, rft);
    yref_yaw(block*2/Ts+1:block*3/Ts-rft)=0.45*ones(1, block/Ts-rft);
    yref_yaw(block*3/Ts-rft+1:block*3/Ts)=-linspace(-0.45, -0.40, rft);
    yref_yaw(block*3/Ts+1:block*4/Ts-rft)=0.40*ones(1, block/Ts-rft);
    yref_yaw(block*4/Ts-rft+1:block*4/Ts)=linspace(0.40, 0.50, rft);
    yref_yaw(block*4/Ts+1:block*5/Ts-rft)=0.50*ones(1, block/Ts-rft);
    yref_yaw(block*5/Ts-rft+1:block*5/Ts)=linspace(0.50, 0.55, rft);
    yref_yaw(block*5/Ts+1:block*6/Ts-rft)=0.55*ones(1, block/Ts-rft);
    
    yref_yaw(block*6/Ts-rft+1:block*6/Ts)=-linspace(-0.55, -0.45, rft);
    yref_yaw(block*6/Ts+1:block*7/Ts-rft)=0.45*ones(1, block/Ts-rft);
    yref_yaw(block*7/Ts-rft+1:block*7/Ts)=-linspace(-0.45, -0.40, rft);
    yref_yaw(block*7/Ts+1:block*8/Ts)=0.40*ones(1, block/Ts);
    
end

if d2_fixed==1
    d2=ones(1, N);
else
    d2=1+0.1*sin(Ts.*n);
end

stdy=std(yref_yaw);
if noisy_case==1
    pwr=SNR/10;
    rstdytstdv=10^(pwr/2);
    wgnoise=(1/rstdytstdv)*stdy*randn(1, N);
end

x=x0;  % states of the system
xlast=x;
e_yaw_last=0; intgrl_e_yaw_last=0; dot_e_yaw_last=0;
trcke=zeros(N, 1); % tracking errors

% control begins..
for n=1:N
    % PID paramertreleri ile rota takibi - kanat açýsý gibi de olabilir
    % kendi kodun

    % produce the control input using PID
    u_yaw(n) = Kp * e_yaw_last + Ki * intgrl_e_yaw_last + Kd * dot_e_yaw_last;
    
    % apply the control signal to the plant
    % system model
    % forward Euler discretization to obtain the vealues of the states
    x(1) = xlast(1) + Ts * (1-xlast(1)-Da1*xlast(1)+Da2*xlast(2)^2);
    x(2) = xlast(2) + Ts * (-xlast(2)+Da1*xlast(1)-Da2*xlast(2)^2-Da3*d2(n)*xlast(2)^2 + u_yaw(n));
    dot_y_yaw(n) = -xlast(3)+Da3*d2(n)*xlast(2)^2;
    x(3) = xlast(3) + Ts * dot_y_yaw(n);
    xlast=x;
    y_yaw(n)=[0 0 1]*x;

    e_yaw(n) = yref_yaw(n) - y_yaw(n);              % tracking error - y(n) - IMU ile ölçülür
    e_yaw_last=e_yaw(n);
    dot_e_yaw(n) = dot_yref_yaw(n) - dot_y_yaw(n);  % derivative of tracking error - dot_y(n) - IMU ile ölçülür - açýsal hýz
    dot_e_yaw_last=dot_e_yaw(n);
    intgrl_e_yaw(n) = sum(e_yaw(1:n));              % integral of tracking error
    intgrl_e_yaw_last=intgrl_e_yaw(n);

    trcke(n)=yref_yaw(n)-y_yaw(n);
    str=sprintf('iter: %d/%d, Tracking error: %g', it, maxiter, trcke(n));
    disp(str)

end

out=trcke'*trcke;