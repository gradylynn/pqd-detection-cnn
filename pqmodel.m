%This file includes a Matlab/Octave function (pqmodel) that implements the power quality model found in 
%R.Igual, C.Medrano, F.J.Arcega, G.Mantescu, "Integral mathematical model
%of power quality disturbances"
%    Copyright (C) 2017  R.Igual,C.Medrano,F.J.Arcega,G.Mantescu

%    This function is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%    Please, CITE THE PAPER referenced above if you use or modify this
%    program

function [out]=pqmodel(ns, fs, f, n, A)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This function generates distorted signals following several power quality
    %models. This function is the implementation of the integral model found in:
    %R.Igual, C.Medrano, F.J.Arcega, G.Mantescu, "Integral mathematical model of power quality disturbances"
    %
    %In this version 29 different classes of distortions are implemented:
    % Class 1 - Pure sinusoidal
    % Class 2 - Sag
    % Class 3 - Swell
    % Class 4 - Interruption
    % Class 5 - Transient/Impulse/Spike
    % Class 6 - Oscillatory transient
    % Class 7 - Harmonics
    % Class 8 - Harmonics with Sag
    % Class 9 - Harmonics with Swell
    % Class 10 - Flicker
    % Class 11 - Flicker with Sag
    % Class 12 - Flicker with Swell
    % Class 13 - Sag with Oscillatory transient
    % Class 14 - Swell with Oscillatory transient
    % Class 15 - Sag with Harmonics
    % Class 16 - Swell with Harmonics
    % Class 17 - Notch
    % Class 18 - Harmonics with Sag with Flicker
    % Class 19 - Harmonics with Swell with Flicker
    % Class 20 - Sag with Harmonics with Flicker
    % Class 21 - Swell with Harmonics with Flicker
    % Class 22 - Sag with Harmonics with Oscillatory transient
    % Class 23 - Swell with Harmonics with Oscillatory transient
    % Class 24 - Harmonics with Sag with Oscillatory transient
    % Class 25 - Harmonics with Swell with Oscillatory transient
    % Class 26 - Harmonics with Sag with Flicker with Oscillatory transient
    % Class 27 - Harmonics with Swell with Flicker with Oscillatory transient
    % Class 28 - Sag with Harmonics with Flicker with Oscillatory transient
    % Class 29 - Swell with Harmonics with Flicker with Oscillatory transient
    %
    %Input Parameters
    %   ns: Number of total signals/samples to generate of each class (default
    %   value: 10). Range of possible values: 1-1000000
    %   fs: Sampling frequency of the signals(default value: 16kHz). Range
    %   of possible values:200-30000Hz
    %   f: Fundamental frequency of the electrical signal (default value:
    %   50Hz). Range of possible values:40-100Hz
    %   n: Number of total cycles (periods) of the fundamental frequency
    %   contained in each sample (default value: 10). Range of possible
    %   values:3-100
    %   A: Amplitude (default value: per unit). Range of possible
    %   values:0.1-400000
    %
    %Output
    %   out: Matrix of dimensions(Number of samples per class, Number of points per signal, Number of classes)
    %   
    %Use example
    %   This function can be called from the Command Window of
    %   Matlab/Octave. 
    %   e.g.: out=pqmodel(10,16000,50,10,1);
    %   A representation of a particular data in the out matrix can be obtained then by typing:
    %   e.g.: plot(out(2,:,18)), all points of the 2nd sample of the 18th
    %   class (Harmonics with Sag with Flicker) are represented
    %%%%%%%%%%%%%%%%%%%%------------------------------------%%%%%%%%%%%%%%%%%%%%%
    
%USER CONFIGURABLE PARAMETERS. If not defined, the default values are provided
    if nargin < 1 || isempty(ns)
      ns = 10;%Default number of samples per class
    end
    if nargin < 2 || isempty(fs)
      fs = 16000;%Default sampling frequency
    end
    if nargin < 3 || isempty(f)
      f = 50;%Default fundamental frequency
    end
    if nargin < 4 || isempty(n)
      n = 10; %Default number of cycles per sample
    end
    if nargin < 5 || isempty(A)
      A = 1; %Default amplitude value, per unit
    end
    
    %Check that configurable parameters are within permitted limits
    if ns<1 || ns>1000000
        disp('Error in function parameter: The number of samples per class must be a value between 1 and 1.000.000')
        out=0;
        return
    end
    if fs<200 || fs>30000
        disp('Error in function parameter: The sampling frequency must be a value between 200Hz and 30kHz')
        out=0;
        return
    end
    if f<40 || f>100
        disp('Error in function parameter: The fundamental frequency must be a value between 40Hz and 100Hz (recommended values are 50Hz or 60Hz)')
        out=0;
        return
    end
    if n<3 || n>100
        disp('Error in function parameter: The number of cycles per sample must be a value between 3 and 100')
        out=0;
        return
    end
    if A<0.1 || A>400000
        disp('Error in function parameter: The normal amplitude must be a value between 0.1V/A and 400kV/kA')
        out=0;
        return
    end
   
%RANDOM PARAMETERS OF THE MODEL
    %General
    phi_min=-pi;
    phi_max=pi;
    %Oscillatory transient and harmonics related
    theta_min=-pi;
    theta_max=pi;    
    %Sag, swell, interruption related
    periodMax=n-1;
    periodMin=1;
    alpha_min=0.1;
    alpha_max=0.9;
    beta_min=0.1;
    beta_max=0.8;
    rho_min=0.9;
    rho_max=1;
    %Transient related
    taPeriod_min=1;
    taPeriod_max=n-1;
    psi_min=0.222;
    psi_max=1.11;
    Onems=round(0.001/(1/fs));%Number of points equivalent to 1ms
    %Flicker related
    ff_min=8;
    ff_max=25;
    lambda_min=0.05;
    lambda_max=0.1;
    %Oscillatory transient related
    tau_min=0.008;
    tau_max=0.04;
    fn_min=300;
    fn_max=900;
    periodMaxOT=n/3.33;
    periodMinOT=0.5;
    pointsfithpartperiod=round(fs/(5*f));
    %Harmonics related
    alpha1=1;
    alpha3_min=0.05;
    alpha3_max=0.15;
    alpha5_min=0.05;
    alpha5_max=0.15;
    alpha7_min=0.05;
    alpha7_max=0.15;
    %Notch related
    k_min=0.1;
    k_max=0.4;
    c=[1,2,4,6];%number of notches per cycle (possible values: 1,2,4 and 6)
    td_min=0;
    tc_min=0;
    tdminustc_min=0.01*(1/f);
    tdminustc_max=0.05*(1/f);
    
%General variables
    N_classes=29; %Total number of different disturbances to simulate
    PointsPerSignal=(fs/f)*n;%number of points contained in a signal to be generated
    t=0:(1/fs):((1/f)*n-(1/fs)); %time vector for all signals
    AllSignals=zeros(ns, PointsPerSignal, N_classes); %output matrix containg the distortions to be generated

%GENERATION OF THE DISTORTIONS
    %% %%%%%%%%%%%Class 1 - Pure sinusoidal%%%%%%%%%%%%%%%%%
    ClassNumber=1; %ID of the class
    for i=1:1:ns %In each iteration a sample is generated
      phi= phi_min+(phi_max-phi_min)*rand; %Select randomly a number between phi_min and phi_max. In general, you can generate N random numbers in the interval (a,b) with the formula r = a + (b-a).*rand(N,1). 
      AllSignals(i,:,ClassNumber)= A*sin(2*pi*f*t-phi);%Generate the signal and store it in the output matrix
    end
    
    %% %%%%%%%%%%%Class 2 - Sag%%%%%%%%%%%%%%%%%
    ClassNumber=2;
    for i=1:1:ns
      [alpha,phi]=selectAmpPhase(alpha_min,alpha_max,phi_min,phi_max);%Select the sag parameters randomly
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);%Select the interval (of points) in which the sag will be applied
      AFinal=A*(1-alpha*u); %The points of the sag will be multiplied by the factor alpha (their amplitude will decrease). The rest, no. 
      AllSignals(i,:,ClassNumber)= AFinal.*sin(2*pi*f*t-phi); %Final distorted sample
    end

    %% %%%%%%%%%%%Class 3 - Swell%%%%%%%%%%%%%%%%%
    ClassNumber=3;
    for i=1:1:ns
      [beta,phi]=selectAmpPhase(beta_min,beta_max,phi_min,phi_max);
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);
      AFinal=A*(1+beta*u); %The points of the sag will be multiplied by the factor beta (their amplitude will increase). The rest, no. 
      AllSignals(i,:,ClassNumber)= AFinal.*sin(2*pi*f*t-phi); %Final distorted sample
    end


    %% %%%%%%%%%%%Class 4 - Interruption%%%%%%%%%%%%%%%%%
    ClassNumber=4;
    for i=1:1:ns
      [rho,phi]=selectAmpPhase(rho_min,rho_max,phi_min,phi_max);
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);
      AFinal=A*(1-rho*u); %The points of the sag will be multiplied by the factor rho (their amplitude will increase much). The rest, no. 
      AllSignals(i,:,ClassNumber)= AFinal.*sin(2*pi*f*t-phi); %Final distorted sample
    end

    %% %%%%%%%%%%%Class 5 - Transient/Impulse/Spike%%%%%%%%%%%%%%%%%
    ClassNumber=5; 
    for i=1:1:ns
      phi=phi_min+(phi_max-phi_min)*rand;% Random number in the range phi_min-phi_max 
      psi=psi_min+(psi_max-psi_min)*rand;
      %select transient interval (ocurrence point: tb-ta)
      taPeriod=taPeriod_min+(taPeriod_max-taPeriod_min)*rand; 
      ta=taPeriod*(1/f);
      pointts=round(taPeriod*(fs/f));
      u1=[zeros(1,pointts) ones(1,PointsPerSignal-pointts)]; %Step funcion, translated pointts to the right
      u2=[zeros(1,pointts+Onems) ones(1,PointsPerSignal-(pointts+Onems))]; %Step funcion, translated (pointts+Onems) points to the right
      u=u1-u2;%Unitary function with all 0 except the interval in which the transient occurs (all 0 except tb-ta)
      %Final distorted sample
      exp_sub=(exp(-750*(t-ta)))-(exp(-344*(t-ta)));
      AFinal=A*psi*exp_sub.*u; 
      AllSignals(i,:,ClassNumber)= -AFinal+A*sin(2*pi*f*t-phi); 
    end

    %% %%%%%%%%%%%Class 6 - Oscillatory transient%%%%%%%%%%%%%%%%%
    ClassNumber=6; 
    for i=1:1:ns
      phi= phi_min+(phi_max-phi_min)*rand;
      [beta,fn,tau,theta]=selectOscTranParam(beta_min,beta_max,fn_min,fn_max,tau_min,tau_max,theta_min,theta_max);%Select oscillatory transiente parameters
      [u,t1]=selectOscTranInterval(fs,f,PointsPerSignal,periodMinOT,periodMaxOT);%Select oscillatory transient occurrence interval and starting and ending points
      AllSignals(i,:,ClassNumber)= A*sin(2*pi*f*t-phi)+A*beta*exp(-(t-t1)/tau).*sin(2*pi*fn*(t-t1)-theta).*u; %Final distorted sample
    end

    %% %%%%%%%%%%%Class 7 - Harmonics%%%%%%%%%%%%%%%%%
    ClassNumber=7;
    for i=1:1:ns
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);%Select harmonics parameters
      alpha7=alpha7_min+(alpha7_max-alpha7_min)*rand;
      theta7=theta_min+(theta_max-theta_min)*rand;
      AllSignals(i,:,ClassNumber)= A*(alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)+alpha7*sin(7*2*pi*f*t-theta7)); %Final distorted sample
    end

    %% %%%%%%%%%%%Class 8 - Harmonics with Sag%%%%%%%%%%%%%%%%%
    ClassNumber=8; 
    for i=1:1:ns
      %Random amplitude and phase parameters selection
      alpha=alpha_min+(alpha_max-alpha_min)*rand;
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Random sag interval selection
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);
      %Final distorted sample
      AFinal=A*(1-alpha*u); 
      AllSignals(i,:,ClassNumber)= A*AFinal.*(alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)); 
    end

    %% %%%%%%%%%%%Class 9 - Harmonics with Swell%%%%%%%%%%%%%%%%%
    ClassNumber=9;
    for i=1:1:ns
      %Random amplitude and phase parameters selection
      beta=beta_min+(beta_max-beta_min)*rand;
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Random swell interval selection
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);
      %Final distorted sample
      AFinal=A*(1+beta*u);
      AllSignals(i,:,ClassNumber)= A*AFinal.*(alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)); 
    end

    %% %%%%%%%%%%%Class 10 - Flicker%%%%%%%%%%%%%%%%%
    ClassNumber=10;
    for i=1:1:ns
      %Random flicker parameters selection
      [lambda,ff,phi]=selectFlickerParam(lambda_min,lambda_max,ff_min,ff_max,phi_min,phi_max);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*sin(2*pi*f*t-phi).*(1+lambda*sin(2*pi*ff*t));
    end

    %% %%%%%%%%%%%Class 11 - Flicker with Sag%%%%%%%%%%%%%%%%%
    ClassNumber=11;
    for i=1:1:ns
      %Random amplitude, phase and frequency parameters selection
      alpha=alpha_min+(alpha_max-alpha_min)*rand;
      [lambda,ff,phi]=selectFlickerParam(lambda_min,lambda_max,ff_min,ff_max,phi_min,phi_max); 
      %Random sag interval selection
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= (A*sin(2*pi*f*t-phi)).*(lambda*sin(2*pi*ff*t)+(1-alpha*u));
    end

    %% %%%%%%%%%%%Class 12 - Flicker with Swell%%%%%%%%%%%%%%%%%
    ClassNumber=12;
    for i=1:1:ns
      %Random amplitude, phase and frequency parameters selection
      beta= beta_min+(beta_max-beta_min)*rand;
      [lambda,ff,phi]=selectFlickerParam(lambda_min,lambda_max,ff_min,ff_max,phi_min,phi_max);  
      %Random swell interval selection
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= (A*sin(2*pi*f*t-phi)).*(lambda*sin(2*pi*ff*t)+(1+beta*u));
    end

    %% %%%%%%%%%%%Class 13 - Sag with Oscillatory transient%%%%%%%%%%%%%%%%%
    ClassNumber=13;
    for i=1:1:ns
      %Random selection of sag and oscillatory transient parameters
      [alpha,phi]=selectAmpPhase(alpha_min,alpha_max,phi_min,phi_max);
      [beta,fn,tau,theta]=selectOscTranParam(beta_min,beta_max,fn_min,fn_max,tau_min,tau_max,theta_min,theta_max);
      %Random sag/osc transient interval selection
      [u,utran,t1tran]=selectOscTranInSagSwellInterval(fs,f,PointsPerSignal,periodMin,periodMax,pointsfithpartperiod);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*sin(2*pi*f*t-phi).*(1-alpha*u)+A*beta*exp(-(t-t1tran)/tau).*sin(2*pi*fn*(t-t1tran)-theta).*utran;
    end

    %% %%%%%%%%%%%Class 14 - Swell with Oscillatory transient%%%%%%%%%%%%%%%%%
    ClassNumber=14;
    for i=1:1:ns
      %Random selection of swell and oscillatory transient parameters
      [beta1,phi]=selectAmpPhase(beta_min,beta_max,phi_min,phi_max);
      [beta2,fn,tau,theta]=selectOscTranParam(beta_min,beta_max,fn_min,fn_max,tau_min,tau_max,theta_min,theta_max);
      %Random swell/osc. transient interval selection
      [u,utran,t1tran]=selectOscTranInSagSwellInterval(fs,f,PointsPerSignal,periodMin,periodMax,pointsfithpartperiod);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*sin(2*pi*f*t-phi).*(1+beta1*u)+A*beta2*exp(-(t-t1tran)/tau).*sin(2*pi*fn*(t-t1tran)-theta).*utran;
    end

    %% %%%%%%%%%%%Class 15 - Sag with Harmonics%%%%%%%%%%%%%%%%%
    ClassNumber=15;
    for i=1:1:ns
      %Random selection of sag and harmonics parameters
      alpha=alpha_min+(alpha_max-alpha_min)*rand;
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Random sag interval selection
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*(alpha1*sin(2*pi*f*t-theta1)+(alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)).*(-alpha*u)); 
    end

    %% %%%%%%%%%%%Class 16 - Swell with Harmonics%%%%%%%%%%%%%%%%%
    ClassNumber=16;
    for i=1:1:ns
      %Random selection of swell and harmonics parameters
      beta=beta_min+(beta_max-beta_min)*rand;
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Random swell interval selection
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*(alpha1*sin(2*pi*f*t-theta1)+(alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)).*(beta*u)); 
    end

    %% %%%%%%%%%%%Class 17 - Notch%%%%%%%%%%%%%%%%%
    ClassNumber=17;
    for i=1:1:ns
        phi=phi_min+(phi_max-phi_min)*rand;
        %Select number of notchs per period randomly
        pos = randi(length(c));
        notch_number=c(pos);
        %Select notch time occurrence
        td_max=(1/(notch_number*f));
        factor=1/(notch_number*f);
        ut=zeros(1,PointsPerSignal);
        k=k_min+(k_max-k_min)*rand;
        tdminustc= tdminustc_min+(tdminustc_max-tdminustc_min)*rand;
        td=td_min+(td_max-td_min)*rand;
        tc=td-tdminustc;
        while(tc<tc_min)
            td=td_min+(td_max-td_min)*rand;
            tc=td-tdminustc;
        end
        %Iterate for all notches in the whole sample
        for nn=0:1:(n)*notch_number-1 
            pointt1=round((tc+factor*nn)*(fs)); %The period is converted into points
            pointt2=round((td+factor*nn)*(fs));
            u1=[zeros(1,pointt1) ones(1,PointsPerSignal-(pointt1))];
            u2=[zeros(1,pointt2) ones(1,PointsPerSignal-(pointt2))];
            u=k*(u1-u2); %Difference function, 0 all except the interval (tc+T*nn,td+T*nn)
            ut=ut+u; %Final unitary function with all 0 except the notch occurrences
        end
        %Final distorted sample
        AllSignals(i,:,ClassNumber)=A*(sin(2*pi*f*t-phi)-sign(sin(2*pi*f*t-phi)).*ut);
    end

    %% %%%%%%%%%%%Class 18 - Harmonics with Sag with Flicker%%%%%%%%%%%%%%%%%
    ClassNumber=18;
    for i=1:1:ns
      %Selection of sag and harmonics parameters randomly
      alpha=alpha_min+(alpha_max-alpha_min)*rand;
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Selection of flicker parameters randomly
      lambda= lambda_min+(lambda_max-lambda_min)*rand;
      ff= ff_min+(ff_max-ff_min)*rand;  
      %Random sag interval selection
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);
      %Final distorted sample
      AFinal=A*(1-alpha*u);
      AllSignals(i,:,ClassNumber)= A*AFinal.*(alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)).*(1+lambda*sin(2*pi*ff*t)); 
    end

    %% %%%%%%%%%%%Class 19 - Harmonics with Swell with Flicker%%%%%%%%%%%%%%%%%
    ClassNumber=19;
    for i=1:1:ns
      %Selection of sag and harmonics parameters randomly
      beta=beta_min+(beta_max-beta_min)*rand;
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Selection of flicker parameters randomly
      lambda= lambda_min+(lambda_max-lambda_min)*rand;
      ff= ff_min+(ff_max-ff_min)*rand;  
      %Random swell interval selection
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);
      %Final distorted sample
      AFinal=A*(1+beta*u);
      AllSignals(i,:,ClassNumber)= A*AFinal.*(alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)).*(1+lambda*sin(2*pi*ff*t)); 
    end

    %% %%%%%%%%%%%Class 20 - Sag with Harmonics with Flicker%%%%%%%%%%%%%%%%%
    ClassNumber=20;
    for i=1:1:ns
      %Selection of sag and harmonics parameters randomly
      alpha=alpha_min+(alpha_max-alpha_min)*rand;
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Selection of flicker parameters randomly
      lambda= lambda_min+(lambda_max-lambda_min)*rand;
      ff= ff_min+(ff_max-ff_min)*rand;  
      %Random sag interval selection
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin); 
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*(alpha1*sin(2*pi*f*t-theta1)+(alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)).*(-alpha*u).*(1+lambda*sin(2*pi*ff*t))); 
    end

    %% %%%%%%%%%%%Class 21 - Swell with Harmonics with Flicker%%%%%%%%%%%%%%%%%
    ClassNumber=21;
    for i=1:1:ns
      %Selection of swell and harmonics parameters randomly
      beta=beta_min+(beta_max-beta_min)*rand;
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Selection of flicker parameters randomly
      lambda= lambda_min+(lambda_max-lambda_min)*rand;
      ff= ff_min+(ff_max-ff_min)*rand;  
      %Random swell interval selection
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin); 
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*(alpha1*sin(2*pi*f*t-theta1)+(alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)).*(beta*u).*(1+lambda*sin(2*pi*ff*t))); 
    end

    %% %%%%%%%%%%%Class 22 - Sag with Harmonics with Oscillatory transient%%%%%%%%%%%%%%%%%
    ClassNumber=22;
    for i=1:1:ns
      %Selection of sag and harmonics parameters randomly
      alpha= alpha_min+(alpha_max-alpha_min)*rand; 
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Selection of oscillatory transient parametras randomly
      [beta,fn,tau,theta]=selectOscTranParam(beta_min,beta_max,fn_min,fn_max,tau_min,tau_max,theta_min,theta_max);
      %Random sag/osc transient interval selection
      [u,utran,t1tran]=selectOscTranInSagSwellInterval(fs,f,PointsPerSignal,periodMin,periodMax,pointsfithpartperiod);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*(alpha1*sin(2*pi*f*t-theta1)+(alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)).*(-alpha*u)+beta*exp(-(t-t1tran)/tau).*sin(2*pi*fn*(t-t1tran)-theta).*utran); %Final distorted function
    end

    %% %%%%%%%%%%%Class 23 - Swell with Harmonics with Oscillatory transient%%%%%%%%%%%%%%%%%
    ClassNumber=23;
    for i=1:1:ns
      %Selection of swell and harmonics parameters randomly
      beta1= beta_min+(beta_max-beta_min)*rand;
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Selection of oscillatory transient parametras randomly
      [beta2,fn,tau,theta]=selectOscTranParam(beta_min,beta_max,fn_min,fn_max,tau_min,tau_max,theta_min,theta_max);
      %Random swell/osc transient interval selection
      [u,utran,t1tran]=selectOscTranInSagSwellInterval(fs,f,PointsPerSignal,periodMin,periodMax,pointsfithpartperiod);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*(alpha1*sin(2*pi*f*t-theta1)+(alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)).*(beta1*u)+beta2*exp(-(t-t1tran)/tau).*sin(2*pi*fn*(t-t1tran)-theta).*utran); %Final distorted function
    end
    
    %% %%%%%%%%%%%Class 24 - Harmonics with Sag with Oscillatory transient%%%%%%%%%%%%%%%%%
    ClassNumber=24;
    for i=1:1:ns
      %Selection of sag and harmonics parameters randomly
      alpha= alpha_min+(alpha_max-alpha_min)*rand; 
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Random sag interval selection
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);
      %Random oscillatory transient interval selection
      [utran,t1tran]=selectOscTranInterval(fs,f,PointsPerSignal,periodMinOT,periodMaxOT);
      %Selection of oscillatory transient parameters randomly
      [beta,fn,tau,theta]=selectOscTranParam(beta_min,beta_max,fn_min,fn_max,tau_min,tau_max,theta_min,theta_max);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*((alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)).*(1-alpha*u)+beta*exp(-(t-t1tran)/tau).*sin(2*pi*fn*(t-t1tran)-theta).*utran); %Final distorted function
    end

    %% %%%%%%%%%%%Class 25 - Harmonics with Swell with Oscillatory transient%%%%%%%%%%%%%%%%%
    ClassNumber=25;
    for i=1:1:ns
      %Selection of swell and harmonics parameters randomly
      beta1= beta_min+(beta_max-beta_min)*rand;
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Random swell interval selection
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);
      %Random oscillatory transient interval selection
      [utran,t1tran]=selectOscTranInterval(fs,f,PointsPerSignal,periodMinOT,periodMaxOT);
      %Selection of oscillatory transient parameters randomly
      [beta2,fn,tau,theta]=selectOscTranParam(beta_min,beta_max,fn_min,fn_max,tau_min,tau_max,theta_min,theta_max);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*((alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)).*(1+beta1*u)+beta2*exp(-(t-t1tran)/tau).*sin(2*pi*fn*(t-t1tran)-theta).*utran); %Final distorted function
    end

    %% %%%%%%%%%%%Class 26 - Harmonics with Sag with Flicker with Oscillatory transient%%%%%%%%%%%%%%%%%
    ClassNumber=26;
    for i=1:1:ns
      %Selection of sag and harmonics parameters randomly
      alpha= alpha_min+(alpha_max-alpha_min)*rand; 
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Random sag interval selection
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);
      %Flicker
      lambda= lambda_min+(lambda_max-lambda_min)*rand;
      ff= ff_min+(ff_max-ff_min)*rand;  
      %Random oscillatory transient interval selection
      [utran,t1tran]=selectOscTranInterval(fs,f,PointsPerSignal,periodMinOT,periodMaxOT);
      %Selection of oscillatory transient parameters randomly
      [beta,fn,tau,theta]=selectOscTranParam(beta_min,beta_max,fn_min,fn_max,tau_min,tau_max,theta_min,theta_max);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*(((alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)).*(1-alpha*u))+beta*exp(-(t-t1tran)/tau).*sin(2*pi*fn*(t-t1tran)-theta).*utran).*(1+lambda*sin(2*pi*ff*t)); %Final distorted function
    end

    %% %%%%%%%%%%%Class 27 - Harmonics with Swell with with Flicker with Oscillatory transient%%%%%%%%%%%%%%%%%
    ClassNumber=27;
    for i=1:1:ns
      %Selection of swell and harmonics parameters randomly
      beta1= beta_min+(beta_max-beta_min)*rand; 
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Random swell interval selection
      u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin);
      %Flicker
      lambda= lambda_min+(lambda_max-lambda_min)*rand;
      ff= ff_min+(ff_max-ff_min)*rand;  
      %Random oscillatory transient interval selection
      [utran,t1tran]=selectOscTranInterval(fs,f,PointsPerSignal,periodMinOT,periodMaxOT);
      %Selection of oscillatory transient parameters randomly
      [beta2,fn,tau,theta]=selectOscTranParam(beta_min,beta_max,fn_min,fn_max,tau_min,tau_max,theta_min,theta_max);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*(((alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)).*(1+beta1*u))+beta2*exp(-(t-t1tran)/tau).*sin(2*pi*fn*(t-t1tran)-theta).*utran).*(1+lambda*sin(2*pi*ff*t)); %Final distorted function
    end

    %% %%%%%%%%%%%Class 28 - Sag with Harmonics with Flicker with Oscillatory transient%%%%%%%%%%%%%%%%%
    ClassNumber=28;
    for i=1:1:ns
      %Selection of sag and harmonics parameters randomly
      alpha= alpha_min+(alpha_max-alpha_min)*rand; 
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Flicker
      lambda= lambda_min+(lambda_max-lambda_min)*rand;
      ff= ff_min+(ff_max-ff_min)*rand;
      %Random sag/osc. transient interval selection
      [u,utran,t1tran]=selectOscTranInSagSwellInterval(fs,f,PointsPerSignal,periodMin,periodMax,pointsfithpartperiod);
      %Selection of oscillatory transient parameters randomly
      [beta,fn,tau,theta]=selectOscTranParam(beta_min,beta_max,fn_min,fn_max,tau_min,tau_max,theta_min,theta_max);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*(alpha1*sin(2*pi*f*t-theta1)+((alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)).*(-alpha*u)+beta*exp(-(t-t1tran)/tau).*sin(2*pi*fn*(t-t1tran)-theta).*utran).*(1+lambda*sin(2*pi*ff*t))); %Final distorted function
    end

    %% %%%%%%%%%%%Class 29 - Swell with Harmonics with Flicker with Oscillatory transient%%%%%%%%%%%%%%%%%
    ClassNumber=29; %ID of the class
    for i=1:1:ns
      %Selection of swell and harmonics parameters randomly
      beta1= beta_min+(beta_max-beta_min)*rand;
      [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max);
      %Flicker
      lambda= lambda_min+(lambda_max-lambda_min)*rand;
      ff= ff_min+(ff_max-ff_min)*rand;
      %Random sag/osc. transient interval selection
      [u,utran,t1tran]=selectOscTranInSagSwellInterval(fs,f,PointsPerSignal,periodMin,periodMax,pointsfithpartperiod);
      %Selection of oscillatory transient parameters randomly
      [beta2,fn,tau,theta]=selectOscTranParam(beta_min,beta_max,fn_min,fn_max,tau_min,tau_max,theta_min,theta_max);
      %Final distorted sample
      AllSignals(i,:,ClassNumber)= A*(alpha1*sin(2*pi*f*t-theta1)+((alpha1*sin(2*pi*f*t-theta1)+alpha3*sin(3*2*pi*f*t-theta3)+alpha5*sin(5*2*pi*f*t-theta5)).*(beta1*u)+beta2*exp(-(t-t1tran)/tau).*sin(2*pi*fn*(t-t1tran)-theta).*utran).*(1+lambda*sin(2*pi*ff*t))); %Final distorted function
    end
    %% %%%%%%---------------------%%%%%%%%
    out= AllSignals; %Output matrix
end

%%%%%%%%%%%%%%%%%%% AUXILIARY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%
% Auxiliary function to generate a random interval for flicker, sag, swell and interruption distortions
% Inputs
%   fs: sampling frequency
%   f: fundamental frequency
%   PointsPerSignal: number of points per signal
%   periodMax: maximum duration of the distortion
%   periodMin: minimum duration of the distortion
% Output
%   u: unitary function, 0 all except the interval in which the distortion
%   occurs
%%%%%%%%%%%%
function u=selectInterval(fs,f,PointsPerSignal,periodMax, periodMin)
      t1t2= periodMin+(periodMax-periodMin)*rand; % Random number between periodMin-periodMax indicating how many periods the disturbace last
      pointst1t2=round(t1t2*(fs/f)); %The period is converted into points
      pointt1=round(0+(PointsPerSignal-0)*rand); %Initial point of the disturbace selected randomly in the range 0-PointsPerSignal
      while (pointt1+pointst1t2)>PointsPerSignal %We check that the disturbace ends before the signal ends, otherwise a new starting point is generated
        pointt1=round(0+(PointsPerSignal-0)*rand);
      end
      u1=[zeros(1,pointt1) ones(1,PointsPerSignal-pointt1)]; %Step funcion, translated pointt1 points to the right
      u2=[zeros(1,pointt1+pointst1t2) ones(1,PointsPerSignal-(pointt1+pointst1t2))]; %Step funcion, translated pointt1+pointst1t2 points to the right
      u=u1-u2; %Difference function, 0 all except the interval (pointt1,pointt1+pointst1t2)
end

%% %%%%%%%%
% Auxiliary function to generate a random value of amplitude and phase
% Inputs
%   amp_min: minimum possible value of the amplitude
%   amp_max: maximum possible value of the amplitude
%   phase_min: minimum possible value of the phase
%   phase_max: maximum possible value of the phase
% Outputs
%   amp: random amplitude value in the interval (amp_min, amp_max)
%   phase: random phase value in the interval (phase_min, phase_max)
%%%%%%%%%%%%
function [amp,phase]=selectAmpPhase(amp_min,amp_max,phase_min,phase_max)
      amp=amp_min+(amp_max-amp_min)*rand;
      phase= phase_min+(phase_max-phase_min)*rand; 
end

%% %%%%%%%%
% Auxiliary function to generate the oscillatory transient parameters
% randomly
% Inputs
%   beta_min: minimum possible value of the beta parameter
%   beta_max: maximum possible value of the beta parameter
%   fn_min: minimum possible value of the fn parameter
%   fn_max: maximum possible value of the fn parameter
%   tau_min: minimum possible value of the tau parameter
%   tau_max: maximum possible value of the tau parameter
%   theta_min: minimum possible value of the theta parameter
%   theta_max: maximum possible value of the theta parameter
% Outputs
%   beta: random value of beta in the interval (beta_min, beta_max)
%   fn: random value of fn in the interval (fn_min, fn_max)
%   tau: random value of tau in the interval (tau_min, tau_max)
%   theta: random value of theta in the interval (theta_min, theta_max)
%%%%%%%%%%%%
function [beta,fn,tau,theta]=selectOscTranParam(beta_min,beta_max,fn_min,fn_max,tau_min,tau_max,theta_min,theta_max)
      beta= beta_min+(beta_max-beta_min)*rand;      
      fn= fn_min+(fn_max-fn_min)*rand; 
      tau= tau_min+(tau_max-tau_min)*rand;
      theta= theta_min+(theta_max-theta_min)*rand;
end

%% %%%%%%%%
% Auxiliary function to select the oscillatory transient interval (duration
% and starting and ending points) randomly
% Inputs
%   fs: sampling frequency
%   f: fundamental frequency
%   PointsPerSignal: number of points per signal
%   periodMaxOT: maximum duration of the distortion (in periods)
%   periodMinOT: minimum duration of the distortion (in periods)
% Outputs
%    u: unitary function, 0 all except the interval in which the distortion
%   occurs
%   t1: starting time of the oscillatory transient distortion
%%%%%%%%%%%%
function [u,t1]=selectOscTranInterval(fs,f,PointsPerSignal,periodMinOT,periodMaxOT)
      t1t2= periodMinOT+(periodMaxOT-periodMinOT)*rand; % Random number between periodMinOT-periodMaxOT, indicating how many periods the disturbace last
      pointst1t2=round(t1t2*(fs/f)); %The period is converted into points
      pointt1=round(0+(PointsPerSignal-0)*rand); %Initial point of the disturbace: random number between 0 and PointsPerSignal
      while (pointt1+pointst1t2)>PointsPerSignal
        pointt1=round(0+(PointsPerSignal-0)*rand);
      end
      t1=pointt1*(1/fs); %convert pointt1 into time
      u1=[zeros(1,pointt1) ones(1,PointsPerSignal-pointt1)]; 
      u2=[zeros(1,pointt1+pointst1t2) ones(1,PointsPerSignal-(pointt1+pointst1t2))]; 
      u=u2-u1;
end

%% %%%%%%%%
% Auxiliary function to select the harmonics parameters randomly
% Inputs
%   alpha3_min: minimum possible value of the alpha3 parameter
%   alpha3_max: maximum possible value of the alpha3 parameter
%   alpha5_min: minimum possible value of the alpha5 parameter
%   alpha5_max: maximum possible value of the alpha5 parameter
%   theta_min: minimum possible value of the theta parameter
%   theta_max: maximum possible value of the theta parameter
% Outputs
%   alpha3: random value of alpha3 in the interval (alpha3_min, alpha3_max)
%   alpha5: random value of alpha5 in the interval (alpha5_min, alpha5_max)
%   theta1: random value of theta1 in the interval (theta1_min, theta1_max)
%   theta3: random value of theta3 in the interval (theta3_min, theta3_max)
%   theta5: random value of theta5 in the interval (theta5_min, theta5_max)
%%%%%%%%%%%%
function [alpha3,alpha5,theta1,theta3,theta5]=selectHarmParam(alpha3_min,alpha3_max,alpha5_min,alpha5_max,theta_min,theta_max)
      alpha3=alpha3_min+(alpha3_max-alpha3_min)*rand;
      alpha5=alpha5_min+(alpha5_max-alpha5_min)*rand;
      theta1=theta_min+(theta_max-theta_min)*rand;
      theta3=theta_min+(theta_max-theta_min)*rand;
      theta5=theta_min+(theta_max-theta_min)*rand;
end

%% %%%%%%%%
% Auxiliary function to generate the flicker parameters randomly
% Inputs
%   lambda_min: minimum possible value of the lambda parameter
%   lambda_max: maximum possible value of the lambda parameter
%   ff_min: minimum possible value of the ff parameter
%   ff_max: maximum possible value of the ff parameter
%   phi_min: minimum possible value of the phi parameter
%   phi_max: maximum possible value of the phi parameter
% Outputs
%   lambda: random value of lambda in the interval (lambda_min, lambda_max)
%   ff: random value of ff in the interval (ff_min, ff_max)
%   phi: random value of tau in the interval (phi_min, phi_max)
%%%%%%%%%%%%
function [lambda,ff,phi]=selectFlickerParam(lambda_min,lambda_max,ff_min,ff_max,phi_min,phi_max)
      lambda=lambda_min+(lambda_max-lambda_min)*rand;
      ff=ff_min+(ff_max-ff_min)*rand;  
      phi=phi_min+(phi_max-phi_min)*rand;
end

%% %%%%%%%%
% Auxiliary function to randomly select the oscillatory transient interval (duration
% and starting and ending points) within a sag/swell interval selected also
% randomly
% Inputs
%   fs: sampling frequency
%   f: fundamental frequency
%   PointsPerSignal: number of points per signal
%   periodMax: maximum duration of the sag/swell distortion (in periods)
%   periodMin: minimum duration of the sag/swell distortion (in periods)
%   pointsfithpartperiod: number of points equivalent to the fith part of a
%   period
% Outputs
%   u: unitary function, 0 all except the interval in which the sag/swell distortion
%   occurs
%   utran: unitary function, 0 all except the interval in which the oscillatory transient distortion
%   occurs
%   t1tran: starting time of the oscillatory transient distortion
%%%%%%%%%%%%
function [u,utran,t1tran]=selectOscTranInSagSwellInterval(fs,f,PointsPerSignal,periodMin,periodMax,pointsfithpartperiod)
      %Random sag or swell interval selection      
      t1t2= periodMin+(periodMax-periodMin)*rand; 
      pointst1t2=round(t1t2*(fs/f)); 
      pointt1=round(0+(PointsPerSignal-0)*rand); 
      while (pointt1+pointst1t2)>PointsPerSignal
        pointt1=round(0+(PointsPerSignal-0)*rand);
      end
      u1=[zeros(1,pointt1) ones(1,PointsPerSignal-pointt1)]; 
      u2=[zeros(1,pointt1+pointst1t2) ones(1,PointsPerSignal-(pointt1+pointst1t2))];
      u=u1-u2; 
      %Random oscillatory transient interval selection
      pointst1t2tran= round(pointsfithpartperiod+(pointst1t2-pointsfithpartperiod)*rand); % Random number to select the duration of the oscillatory transient (it has to be whithin the sag duration)
      pointt1tran=round(pointt1+((pointt1+pointst1t2-pointst1t2tran)-pointt1)*rand); %Initial point of the oscillatory transient disturbace
      t1tran=pointt1tran*(1/fs); %Convert intial oscillatory transient point into time
      u1=[zeros(1,pointt1tran) ones(1,PointsPerSignal-pointt1tran)]; %Step funcion, translated pointt1tran points to the right
      u2=[zeros(1,pointt1tran+pointst1t2tran) ones(1,PointsPerSignal-(pointt1tran+pointst1t2tran))]; %Step funcion, translated pointt1tran+pointst1t2tran points to the right
      utran=u2-u1; %Difference function, 0 all except the interval (pointt1tran,pointt1tran+pointst1t2tran)
end





