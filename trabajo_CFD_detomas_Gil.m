function deTomas_Gil
%Based on 
close all
clear all
clc

%%
%Data for the model
mumax=1/240;
Ks=20/1000 %mg CODsL
Ki=100/1000%mg CODs L
Y=0.67;
%Initital conditions
Init=[2000/1000 200/1000];
%tspan
opciones=odeset('abstol', 1e-7, 'reltol', 1e-7);
tspan=linspace(0,50);
[t,yval]=ode45(@andrews,tspan,Init,opciones)
%OUR
figure(1)
plot(t,yval(:,1),'-g')
figure(2)

plot(t,yval(:,2),'-r')


%plot(t,yval(:,3),'-b')

    function dvar=andrews(t,y)
        X=y(1);
        S=y(2);
        %dvar
        dvar=zeros(2,1);
        dvar(1)=(mumax*S*X)/((S+Ks)*(1+(S)/Ki));
        dvar(2)=-1*dvar(1)/Y;

        our=((1-Y)*dvar(1)/Y); 
        if our<0
            our=0;
        end

        figure(3)
        plot(t,our,'*r')
        hold on  
    end


end