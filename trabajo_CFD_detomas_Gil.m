function deTomas_Gil
%Based on 
close all
clear all
clc


%%
%Data for the model
mumax=1/240;
Ks=20 %mg CODsL
Ki=100%mg CODs L
Y=0.67;
%Initital conditions
Init=[2000 200];
%tspan
opciones=odeset('abstol', 1e-7, 'reltol', 1e-7);
tspan=linspace(0,50);
[t,yval]=ode45(@andrews,tspan,Init,opciones)

%%
for a=1:length(t)
    b=t(a);
    c=yval(a,:);
    [dvar,var(a)]=andrews(b,c);
end
v=[((1-Y)*var./Y)];
for i=1:length(v)
    oure_(i)=normrnd(v(i),0.03)
end

figure(6)
plot(t,v,'-k',t,oure_,'-o')
%OUR
figure(1)
plot(t,yval(:,1),'-g')
figure(2)
plot(t,yval(:,2),'-r')
x=yval(:,1);


    function [dvar,b]=andrews(t,y)
        X=y(1);
        S=y(2);
        %dvar
        dvar=zeros(2,1)
        dvar(1)=(mumax*S*X)/((S+Ks)*(1+(S)/Ki));
        dvar(2)=-1*dvar(1)/Y;    
       b=dvar(1)
        our=[((1-Y)*dvar(1)./Y)]
        
        if our<0
            our=0
        end
        oure=normrnd(our,0.03);
        
    end

end
