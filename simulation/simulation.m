%% Andreas Vezakis
%% AM : 03117186
clear;
close all;
clc;

l = zeros(1,7);
l(4) = 10.0; 
l(6) = 7.0;  
l(7) = 5.0;
dt = 0.1;    
Tf = 20.0; 
t=0:dt:Tf; 

for i=0:3
    if (mod(i,2)==0)
        p_A=[11.00;12.00;20.00]; 
        p_B=[1.00;1.00;20.00];
    else                     
        p_B=[11.00;12.00;20.00];
        p_A=[1.00;1.00;20.00];
    end


    a1 = 1; 
    b1 = 1;
    c1 = 1;

    a2 = (3/Tf^2)*(p_B(1) - p_A(1));
    a3 = -(2/Tf^3)*(p_B(1) - p_A(1));

    b2 = (3/Tf^2)*(p_B(2) - p_A(2));
    b3 = -(2/Tf^3)*(p_B(2) - p_A(2));

    c2 = (3/Tf^2)*(p_B(3) - p_A(3)); 
    c3 = -(2/Tf^3)*(p_B(3) - p_A(3)); 
    xd = zeros(length(t)); 
    yd = zeros(length(t));
    xd(1) = p_A(1); 
    yd(1) = p_A(2); 
    zd(1) = p_A(3); 
    for k=2:length(t)    
       xd(k) = p_A(1) + a2*t(k)^2 + a3*t(k)^3;    
       yd(k) = p_A(2) + b2*t(k)^2 + b3*t(k)^3;
       zd(k) = p_A(3) + c2*t(k)^2 + c3*t(k)^3;
    end  
    disp(sprintf('Kinematic Simulation number %d\n', i)); 
    
    for m=1:length(t)
        pe_x=xd(m);
        pe_y=yd(m);
        pe_z=zd(m);
        
        qd3(m) = acos((pe_x^2 + pe_y^2 + pe_z^2 - l(4)^2 - l(6)^2 - l(7)^2)/(2 * l(6) * l(7)));
        qd2(m) = acos((sin(qd3(m)) * l(7) * pe_y + (l(6) + l(7) * cos(qd3(m))) * sqrt(l(6)^2 + l(7)^2 + 2 * cos(qd3(m)) * l(6) * l(7) - pe_y^2))/(l(6)^2 + l(7)^2 + 2 * cos(qd3(m)) * l(6) * l(7)));
        c23 = cos(qd2(m) + qd3(m));
        qd1(m) = acos(((l(6) * cos(qd2(m)) + l(7) * c23) * (-1) * pe_z + l(6) * sqrt((l(7) * cos(qd2(m)) + l(7) * c23)^2 + l(6)^2 - pe_z^2))/((l(7) * cos(qd2(m)) + l(7) * c23)^2 + l(4)^2));
    end
            
    for it=1:length(t)
        xd1(it) = 0;  
        yd1(it) = 0; 
        zd1(it) = 0;

        xd2(it) = cos(qd2(it)) * sin(qd1(it)) * l(6); 
        yd2(it) = l(6) * sin(qd2(it));
        zd2(it) = -cos(qd1(it)) * cos(qd2(it)) * l(6);
        
        xde(it) = cos(qd1(it)) * l(4) + sin(qd1(it)) * cos(qd2(it) + qd3(it)) * l(7) + sin(qd1(it)) * cos(qd2(it)) * l(6);
        yde(it) = sin(qd2(it) + qd3(it)) * l(7) + sin(qd2(it)) * l(6);
        zde(it) = sin(qd1(it)) * l(4) - cos(qd1(it)) * cos(qd2(it)) * l(6) - cos(qd1(it)) * cos(qd2(it) + qd3(it)) * l(7);
    end
    fig1 = figure; 
    subplot(2,3,1); 
    plot(t,xde(:)); 
    ylabel('Pe,x (cm)'); 
    xlabel('Time (sec)');  

    subplot(2,3,2); 
    plot(t,yde(:)); 
    ylabel('Pe,y (cm)'); 
    xlabel('Time (sec)');
    title('Desired Position of the End-Effector');

    subplot(2,3,3); 
    plot(t,zd(:)); 
    ylabel('Pe,z (cm)'); 
    xlabel('Time (sec)');
    
    Pe=[xde(:) yde(:) zde(:)];
    fig2 = figure;  
    Ve = zeros(length(t),3);
    for q=2:length(t)
        Ve(q,:) = Pe(q,:) - Pe(q-1,:);
    end

    subplot(2,3,1); 
    plot(t,Ve(:,1)*(length(t)-1)*180/(Tf*pi)); 
    ylabel('Vx (cm/sec)'); 
    xlabel('Time (sec)');  

    subplot(2,3,2); 
    plot(t,Ve(:,2)*(length(t)-1)*180/(Tf*pi));  
    ylabel('Vy (cm/sec)'); 
    xlabel('Time (sec)');   
    title('Linear Velocity of the End-Effector');

    subplot(2,3,3); 
    plot(t,Ve(:,3)*(length(t)-1)*180/(Tf*pi)); 
    ylabel('Vz (cm/sec)'); 
    xlabel('Time (sec)');  
    
    fig3=figure;
    subplot(2,3,1); 
    plot(t,qd1(:)*180/pi); 
    ylabel('q1 (degrees)'); 
    xlabel('Time (sec)');  

    subplot(2,3,2); 
    plot(t,qd2(:)*180/pi); 
    ylabel('q2 (degrees)'); 
    xlabel('Time (sec)');   
    title('Angles of the Joints');

    subplot(2,3,3); 
    plot(t,qd3(:)*180/pi); 
    ylabel('q3 (degrees)'); 
    xlabel('Time (sec)'); 


    qd=[qd1(:) qd2(:) qd3(:)];
    fig4 = figure;  
    subplot(1,3,1);
    u = zeros(length(t),3);
    for q=2:length(t)
        u(q,:) = qd(q,:) - qd(q-1,:);
    end
    
    subplot(2,3,1); 
    plot(t,u(:,1)*(length(t)-1)*180/(Tf*pi)); 
    ylabel('ù1 (deg/sec)'); 
    xlabel('Time (sec)');  

    subplot(2,3,2); 
    plot(t,u(:,2)*(length(t)-1)*180/(Tf*pi));  
    ylabel('ù2 (deg/sec)'); 
    xlabel('Time (sec)');   
    title('Angular Velocity of the Joints');

    subplot(2,3,3); 
    plot(t,u(:,3)*(length(t)-1)*180/(Tf*pi)); 
    ylabel('ù3 (deg/sec)'); 
    xlabel('Time (sec)');
    fig2 = figure;  
    axis([-l(4) l(4) -(l(6)+l(7)) (l(6)+l(7)) 0 25]) %setting axis
    axis on 
    hold on 
    xlabel('x (cm)'); 
    ylabel('y (cm)'); 
    zlabel ('z (cm)');
    plot3([0],[0], [0] ,'o');
    dtk=20; 
    for tk=1:dtk:length(t) 
           pause(1);   
           plot3([0,xd1],[0,yd1], [0,zd1]);
           plot3([xd1],[yd1],[zd1],'o');  
           plot3([xd1(tk),xd2(tk)],[xd1(tk),yd2(tk)], [zd1(tk),zd2(tk)]);
           plot3([xd2(tk)],[yd2(tk)], [zd2(tk)],'o');  
           plot3([xd2(tk),xde(tk)],[yd2(tk),yde(tk)],[zd2(tk),zde(tk)]);	
           plot3([xde(tk)],[yde(tk)],[zde(tk)],'r*');   
    end

end
