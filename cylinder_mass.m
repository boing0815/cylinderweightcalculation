%MATLAB® example code for calculation of cylinder mass
% Author Tobias Radermacher, TU Dresden, Chair of Fluid-Mechatronic Systems
% License: GNU GPL
% download: github.com/boing0815
clc,clear, close all
%calculation of mass for F*l=const
sig=235E6; %N/m^2,max tension,reversibel
Re=sig;
sf=1.5;%security factor
alpha=2; %ratio
leg="";
l_light=[];
force=700; %niveau of force for f*l
 
f=figure(); hold on; grid on;
% calculation of lines with F*l=const
l=[0.1:0.3:1.6];
skal=0.02:0.02:4; %scaling ofo F*l
for h=1:length(l)
    ls=l(h)./skal;
    fs=force.*skal;
    for k=1:length(skal) %for every variant with f*l constant
        sp(k)=pfl(fs(k)*1000,ls(k),alpha,sf); %neccesary pressure level in bar
        m(k)=zylindermasse_von_r(sig,Re,sp(k)*1E5,fs(k)*1000,ls(k),alpha,sf); %calc cyl. mass; change force to N
    end
    [m_min_f_mal_l(h),ind(h)]=min(m);
    p_min_f_mal_l(h)=sp(ind(h));
    l_light(h)=ls(ind(h));
    f_light(h)=fs(ind(h));
    
    pl=plot(sp,m);
    uistack(pl,"bottom");
    test=convertCharsToStrings(sprintf('%.0f', force*l(h)));
    leg=[leg,test];
end
title(["mass for equivalent lightest cylinders","\alpha=2, sf=1.5, \sigma=235 N/mm^2"]);
xlabel('pressure level [bar]');
ylabel('mass [kg]');
xlim([0,1800])
ylim([0,600])
%marking min.
plot(p_min_f_mal_l,m_min_f_mal_l,'*')
%scale position for paper size
f.Position(3:4) = [400 300];
%make legend
le=legend(leg(1,2:end));
le.String=flip(le.String); %flip legend
le.Position = [0.4 0.5 0.1 0.2];
saveas(f,"lightest cylinders_mass.fig");
saveas(f,"lightest cylinders_mass.emf");

function p=pfl(F,l,alpha, sf)
    %calc of p for defined area ratio
    E=2.1E10; %E-module in N/m^2;
    %go for it
    % min. radius der of rod against buckling
    % buckling length 90%
    r_s=(F*sf*2*(0.9*l)^2/(E*pi^3))^(1/4);
    r_i=sqrt(alpha/(alpha-1))*r_s; %area ratio   
    p=F/(r_i^2*pi);
    p=p/1E5; %in bar
end

function [m] = zylindermasse_von_r(sig,Re,p,F,l,alpha,sf)
%change mass to 1E-12, if useless geometry
% sig, Re in N/m^2, p in N/m^2, F in N, l in m
 
%input:
    %mat. properties
    rho=7850;%kg/m^3, dens. steel
    rho_oel=880; %kg/m^3
    E=2.1E10; %E-module in N/m^2
%calculation
    %rod radius
    r_s=(F*sf*2*(0.9*l)^2/(E*pi^3))^(1/4);
    %piston radius
    r_i=sqrt(F/(p*pi));
    %wall thickness (fat, plasticity=yes)
    r_a=r_i/sqrt(1-p*sf*sqrt(3)/sig);
    s=r_a-r_i;
    %wall thickness (thin, elastic):
    %s=r_i*p/(sig-p)
%cap heigth 
    h=(r_i+s)*sqrt(0.8*p/(sig/sf));
%finally: the result:
    %mass:
    m_platten=rho*pi*(4*h*r_a^2);
    m_rohr=rho*pi*(2*r_a*s-s^2)*l;
    m_stange=rho*pi*l*r_s^2;
    mz=m_platten+m_rohr+m_stange;
    m_oel=(l-h)*r_i^2*pi*rho_oel;
    m=mz+m_oel;
    if s>0.1 || s<=0 %delete useless wall values
        m = NaN;
    end
end

