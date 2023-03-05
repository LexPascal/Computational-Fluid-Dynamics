function Trabajo_Final_GildeTomas_pruebas1New
clear all
close all
clc

%datos
mumax=1/240;
Ks=20;
Ki=100;
Y=0.67;
X0=2000;
S0=200;
tf=50;

%funcionesi
Cinit=[X0 S0];
tr=linspace(0,tf,100);
opciones=odeset('abstol', 1e-5, 'reltol', 1e-5);
[t,var]=ode45(@andrews,tr,Cinit,opciones);              %llamar a la ODE
Xp=var(:,1); Sp=var(:,2);
mu=(mumax.*Sp)./(Sp+Ks+((Sp.^2)/Ki));
dXdtn=mu.*Xp;
OURn=((1-Y)/Y)*dXdtn;

error=0.03;
OURdata=(OURn+randn(size(OURn))*error);                 %Fake del our experimental

%iter number of mesurments
    for j=1:10
    
OURexp2=OURdata(1:j:end);
tr2=tr(1:j:end);
N=length(tr2);
Q=norm(error*OURdata).^2/(N-3);
texp_n=linspace(0,30);

for i=1:length(texp_n)
Qour_real(i)=0.20;
end
for a=1:length(OURdata)
    error_(a)=(0.20-OURdata(a))^2;
end
sse=cumsum(error_)
p=3
for b=1:length(OURdata)
    N_=b
    SSE_=sse(b)
    ssquared(b)=SSE_/(N_-p)
    Q_(b)=1/ssquared(b)
end


%Figure 3
texpn=linspace(0,30);
for i=1:length(texpn)
    Qour_real(i)=0.20;
end
OUR_exp=normrnd(Qour_real,error) 


%Grafica X %figura 3
%Figure 3 - COMPROBAR!!
texpn=linspace(0,30);
for i=1:length(texpn)
    Qour_real(i)=0.20;
end
errorplot=0.01;
OUR_exp=normrnd(Qour_real,errorplot);


%OUR_exp=normrnd(Qour_real,error) 
%figure(99)
%plot(texpn, Qour_real, texpn, OUR_exp, '*')
%alternativa figure 3


%DATOS EXPERIMENTALES
texp=linspace(0,50);

%Parámetros óptimos
[pop]=runnested([1/240 20 100],OURexp2)
mumax=pop(1);
Ks=pop(2);
Ki=pop(3);
%{
popt=fminsearch(@fobj,[1/240 20 100]);
mumax=popt(1);
Ks=popt(2);
Ki=popt(3);
%}

mumax_(j)=mumax;
Ks_(j)=Ks;
Ki_(j)=Ki;


OURmod=((1-Y)/Y)*dXdtn;

%% 
%Sensibilidad
%%
%mumax
tsens=linspace(0,50);
mumaxi=mumax;
mumax=mumaxi*1.01;
[t,yh]=ode45(@andrews,tsens,Cinit,opciones);
OURh=calculateOUR(yh);
mumax=mumaxi*0.99;
[t,yl]=ode45(@andrews,tsens,Cinit,opciones);
OURy=calculateOUR(yl);
mumax=mumaxi;  
SensMU=(OURh(:,1)-OURy(:,1))/(mumax*0.02);

%%
%ks
Ksi=Ks;
Ks=Ksi*1.01;
[t,yh]=ode45(@andrews,tsens,Cinit,opciones);
OURh=calculateOUR(yh);
Ks=Ksi*0.99;
[t,yl]=ode45(@andrews,tsens,Cinit,opciones);
OURy=calculateOUR(yl);
Ks=Ksi;
SensKS=(OURh(:,1)-OURy(:,1))/(Ks*0.02);

%%
%ki
Kii=Ki;
Ki=Kii*1.01;
[t,yh]=ode45(@andrews,tsens,Cinit,opciones);
OURh=calculateOUR(yh);
Ki=Kii*0.99;
[t,yl]=ode45(@andrews,tsens,Cinit,opciones);
OURy=calculateOUR(yl);
Ki=Kii;
SensKi=(OURh(:,1)-OURy(:,1))/(Ki*0.02);

%%
%Calculate the sensitibity matrix
sens=[SensMU SensKS SensKi]; 

%FIM
%Dos columnas,seis filas
% 6. Calculation of the FIM ( Fisher Matrix Information)
%to do --> sens
fim=zeros(3,3);
for i=1:length(sens)
    fim=fim+sens(i,:)'*1./Q*sens(i,:);
end

%mumax,ks,ki
standard_error=abs(sqrt(diag(inv(fim)))) 

frac_mu(j)=(standard_error(1)/mumax)*100
frac_ks(j)=(standard_error(2)/Ks)*100
frac_ki(j)=(standard_error(3)/Ki)*100

error_mu(j)=standard_error(1);
error_ks(j)=standard_error(2);
error_ki(j)=standard_error(3);

end

  %%

  %definición de función objetivo
    function fval=runnested(p,OUR_c)
        fval=fminsearch(@fobj_nest,p);
        
            function fval=fobj_nest(p)
                mumax=p(1); Ks=p(2); Ki=p(3);

                [tm,varm]=ode45(@andrews,tr2,Cinit,opciones);
                Xp_opt=varm(:,1); Sp_opt=varm(:,2);
                mu_opt=(mumax.*Sp_opt)./(Sp_opt+Ks+((Sp_opt.^2)/Ki));
                dXdtn_opt=mu_opt.*Xp_opt;
                OURmod=((1-Y)/Y)*dXdtn_opt;
                fval=norm(OURmod-OUR_c);
            end
    end
function f=fobj(p)
    mumax=p(1); Ks=p(2); Ki=p(3);
    [tm,varm]=ode45(@andrews,texp,Cinit,opciones);
		Xp_opt=varm(:,1); Sp_opt=varm(:,2);
    mu_opt=(mumax.*Sp_opt)./(Sp_opt+Ks+((Sp_opt.^2)/Ki));
		dXdtn_opt=mu_opt.*Xp_opt;
    OURmod=((1-Y)/Y)*dXdtn_opt;
    f=norm(OURmod-OURdata);
end

function [OUR]=calculateOUR(var)
    %Calculate the our for each andrews
    
    Xp=var(:,1); Sp=var(:,2);
    mu=(mumax.*Sp)./(Sp+Ks+((Sp.^2)/Ki));
    dXdtn=mu.*Xp;
    OUR=((1-Y)/Y)*dXdtn;
end

% funcion Andrews
function [dvar OUR]=andrews(t,var)
        X=var(1);
        S=var(2);

        mu=(mumax*S)/(S+Ks+((S^2)/Ki));

        dvar=zeros(2,1);
        dvar(1)=mu*X;
        dvar(2)=(-1/Y)*dvar(1);  
end

%%
%if needed 
ERROR=OURdata-OURmod;
SSE=(sum(ERROR))^2;
N=length(ERROR);
p=3;
s2=(SSE)/(N-p);
S2_IN=1/s2;

medidas=linspace(0,5,10);
one_medidas=ones(length(medidas));

mumax_real=6;
mumax_real=mumax_real*one_medidas;
ks_real=20;
ks_real=ks_real*one_medidas;
ki_real=100;
ki_real=ki_real*one_medidas;
%
mumax_=mumax_*(60*24)
Ks_;
%hacemos vectores tan largos como m=10

ki_realsum=(ki_real(1)+ki_real(2)+ki_real(3)+ki_real(4)+ki_real(5)+ki_real(6)+ki_real(7)+ki_real(8)+ki_real(9)+ki_real(10));

%PLOTS
figure(1)
plot(t,OURn)
figure(2)
plot(t,Sp,'r')
figure(3)
plot(t,Xp,'g')

%FIGURE 2
figure(4)
line(t,OURn)
hold on
plot(t,OURdata,'.b')

%FIGURE 3
figure(5)
plot(texpn, Qour_real, texpn, OUR_exp, '--')
xlim([0 30])
ylim([0.10 0.30])

%FIGURE 4
figure(6)
plot(medidas,mumax_real,'k',medidas,mumax_,'*')
hold on
plot(medidas,ks_real,medidas,Ks_,'*')
plot(medidas,ki_real,medidas,Ki_,'*')

%Figure 5
figure(7)
plot(tsens, SensKi, tsens, SensKS)
figure(8)
plot(tsens, SensMU)

%figure 4 b --> ploted as 8
figure(9)
plot(medidas,frac_mu,'k*')
hold on
plot(medidas,frac_ks,'r*')
plot(medidas,frac_ki,'b*')
legend('mumax','Ks','Ki')





end