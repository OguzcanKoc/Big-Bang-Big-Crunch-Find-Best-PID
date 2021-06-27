clear all
close all
warning off
clc

PID_min=zeros(3, 1);
PID_max=ones(3, 1) * 5;
objfunc='PID_SSE';
popnum=10;              % population number
ds=10;                  % decrease scale of the spread range
sc=10;                  % stop ciriterion

% bestPID: the P, I, D values which give the best tracking performance
% besttrckngprfmnc: best tracking performance with minimum sum of squared
% errors - SSE
[bestPID, besttrckngprfmnc]=BBBC(PID_min, PID_max, objfunc, popnum, ds, sc);
Kp=bestPID(1); Ki=bestPID(2); Kd=bestPID(3);
Kp
Ki
Kd
save 'bestPID.mat' bestPID