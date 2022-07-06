%MECH307 Project 1(Inside a Rectangular) by Kagan Ucak, 0069435
%23.02.2021

clc, clear, close all
%% Needed parameters. 
th = 0:pi/50:2*pi;                          %angle for drawing circiles
txt=["k","b","m","r","y","g","w","c","k"]; %color
dt=0.0001;
fps= linspace(0,1,1/dt +1);                   %the precise of simulation. dt is equal to 0.0001 s.
n=1;                                        %number of collusion
cons=0;                                     %a paramater to save energy values to array "data"
data=[0 0; 1 1];                            %for energy values
data_2=[0 0 0 0 0; 1 1 1 1 1];              %for save time and momentum paramaters before and after any collision, also to prevent the errors
sp=100;                                     %to increase the speed of simulation
N= randi([2 8],1);  
N=input("Number of ball: "); %If you wish you can define a number.
v_N=zeros(1,N);                             %normal velocity initial array
v_T=zeros(1,N);                             %for tangantial velocity initial array

figure(1)
pause(3)
%% To define the positions of center of balls. 
s=[0.5 1;0.5 1];                  %the initial position coordinates
r=randi([1000 3000],[1 N])/10000; %To define radius of balls. 
M=r.^2;                           %To define the mass of balls. 
for i=[3:N]                         %Because I have already assigned two balls. 
    s_1=randi([300 2700],1);        
    s_2=randi([300 1700],1);
    s(1,i)= s_1/1000;               %X directional position of ball "i"
    s(2,i)= s_2/1000;               %Y directional position of ball "i"
    j=0;
    while j~=i                      %If ball "i" overlaps with previous balls, the position is reassigned.
        j=j+1;
        while ((s(1,i)-s(1,j))^2 +(s(2,i)-s(2,j))^2)<=(r(i)+r(j))^2 && i~=j
            s_1=randi([300 2700],1);
            s_2=randi([300 1700],1);
            s(1,i)= s_1/1000;
            s(2,i)= s_2/1000;
            j=0;
            break
        end
    end
end

%% To define X and Y directional velocities. 
v= randi([round(1000*sqrt(2)) round(2000*sqrt(2))],[1 N])/100;
delta=rand(1,N)*2*pi;

v_x=zeros(1,N);
v_y=zeros(1,N);
for i=[1:N]
    v_x(i)=v(i).*sin(delta(i));
    v_y(i)=v(i).*cos(delta(i));
end


%% To analyze every dt second situation
for p =fps
    for i= [1:N]            %step by step every ball
        s(1,i)= s(1,i)+(dt*v_x(i));
        s(2,i)= s(2,i)+(dt*v_y(i));
        if mod(p,dt*sp)==0    %to cpu and time saving
            hold on  
            grid on
            set(gcf, 'WindowState', 'maximized'); %in order to maximize the window
            rectangle('Position',[0 0 3 2],'EdgeColor','r', 'LineWidth',3)
            xunit = r(i) * cos(th) + s(1,i);
            yunit = r(i)* sin(th) + s(2,i);
            plot1(i)=plot(xunit, yunit);
            fill_1(i)=fill(xunit, yunit,txt(i));
            title("Time: " + p + "(s)"+ ",  The Number of Collisions: "+(n-1))
            axis equal
            axis([-1 4 -1 3])

            e_1=s(1,i)+v_x(i);              %to draw arrows
            e_2=s(2,i)+v_y(i);
            q1(i)=quiver(s(1,i),s(2,i), 0.01*v_x(i) ,0.01*v_y(i), 'linewidth',2);
        end
        %% If ball will hit wall after one dt time, balls change one direction of their velocities.  
        if (s(1,i)+(dt*v_x(i))+r(i) >=3 ||s(1,i)+(dt*v_x(i))-r(i) <=0) %to avoid repetition one dt time before
            v_x(i)=-v_x(i);
        end
        if (s(2,i)+(dt*v_y(i))+r(i) >=2 || s(2,i)+(dt*v_y(i))-r(i)<=0) 
            v_y(i)=-v_y(i);
        end
    end
    %% To calculate velocities of balls after collisions
    b_v_x= v_x;             
    b_v_y= v_y;
    for i= [1:N-1]          %to take a ball
        for i_2= [i+1:N]    %to compare it with another ball
            if ((s(1,i)-s(1,i_2))^2 +(s(2,i)-s(2,i_2))^2)<=(r(i)+r(i_2))^2 %if the distance between their radius is small then dum of their radius, they collect
                n=n+1;
                b_v_x(i)= v_x(i); %to bake up the first velocities if the case of repetition
                b_v_y(i)= v_y(i);
                b_v_x(i_2)=v_x(i_2);
                b_v_y(i_2)=v_y(i_2);
                %% To change axis
                teta= atan((s(2,i)-s(2,i_2))/(s(1,i)-s(1,i_2))); %the angle between the normal axis and X axis. 
              
                v_N(i)=v_x(i)*cos(teta)+sin(teta)*v_y(i); 
                v_T(i)=v_y(i)*cos(teta)-v_x(i)*sin(teta);
                v_N(i_2)=v_x(i_2)*cos(teta)+sin(teta)*v_y(i_2);
                v_T(i_2)=v_y(i_2)*cos(teta)-v_x(i_2)*sin(teta);
                %% to calculate new velocities after collections 
                b2_v_N= v_N; %because v_N(i) will change, purpuse f b2_V_N is to calculate v_N(i_2)
                b2_v_T= v_T;
                
                v_N(i)=M(i)/(M(i)+M(i_2))*(v_N(i)+M(i_2)/M(i)*v_N(i_2)-M(i_2)/M(i)*(v_N(i)-v_N(i_2)));
                v_N(i_2)= M(i_2)/(M(i)+M(i_2))*(v_N(i_2)+M(i)/M(i_2)*b2_v_N(i)-M(i)/M(i_2)*(v_N(i_2)-b2_v_N(i)));
                %% rotate N and T axis to x and y axis
                v_x(i)= cos(teta)*v_N(i)-v_T(i)*sin(teta);
                v_y(i)= sin(teta)*v_N(i)+v_T(i)*cos(teta);
                v_x(i_2)= (cos(teta)*v_N(i_2)-v_T(i_2)*sin(teta));
                v_y(i_2)= (sin(teta)*v_N(i_2)+v_T(i_2)*cos(teta));
                %% to prevent repetitions,
                data_2(n,1)=b_v_x(i)*M(i)+b_v_x(i_2)*M(i_2); %x momentum values before the collision
                data_2(n,2)=v_x(i)*M(i)+v_x(i_2)*M(i_2);     %x momentum values after the collision
                data_2(n,3)=b_v_y(i)*M(i)+b_v_y(i_2)*M(i_2); %y momentum values before the collision
                data_2(n,4)=v_y(i)*M(i)+v_y(i_2)*M(i_2);     %y momentum values after the collision
                data_2(n,5)=p;                               %the time of collision
                
                if  round((data_2(n-1,1:4)*10000))== round(10000*data_2(n,1:4)) %if momentum values of consecutive collisions is exactly same, it is probably a repetition
                    if round(10000*(data_2(n-1,5)+dt))<=round(10000*data_2(n,5)) %if the time differences between these collisions is less than or equal to dt. It is most probably a repetition. 
                        if round(10000*(data_2(n-1,5)+dt*10))>round(10000*data_2(n,5)) % If the time difference between these collisions is very much, there may not be a repetition. So I eliminate that. 
                            v_x(i)= b_v_x(i); %velocity values turn back before collision
                            v_y(i)= b_v_y(i) ;
                            v_x(i_2)=b_v_x(i_2);
                            v_y(i_2)=b_v_y(i_2);
                            n=n-1;              % number of collision
                        end
                    end
                end                    
            end
        end
    end
    
    %% to collect energy values at time p
    cons=cons+1;
    data(cons,1)=sum(M.*v_x.^2 + M.*v_y.^2)/2;
    
    if mod(p,dt*sp)==0
         pause(0.00001)
    end
     delete(fill_1) %to refresh figure without closing
     delete(plot1)
     delete(q1)
end

%% to plot energy values
clf
plot(fps, data(:,1),"linewidth",3)
title("Sum of energy of the system")
xlabel("Time(s)")
ylabel("Energy(J)")
axis equal
grid on