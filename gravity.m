% 1 %Intialisation of values%  
Re     =600;
Nx     = 32; 
Ny     = 32;
U_inf  = 1;
Ra     =1000;
lx     = 0.1;
ly     = 0.1;
nu     = U_inf*lx/Re;
dx     = lx/Nx; 
dy     = ly/Ny; 
b      = dx/dy;
alpha  = 0.1;
accuracy1=10^-3;
accuracy2=10^-1;
Pr=1.7; %assumed


 %2 %Initial Vaues matrices  
    w   = ones(Ny+1, Nx+1);
    psi = zeros(Ny+1, Nx+1);
    u   = zeros(Ny+1, Nx+1);
    v   = zeros(Ny+1, Nx+1);
    T   = zeros(Ny+1, Nx+1);
    psi_actual = zeros(Ny+1, Nx+1);
    w_actual = ones(Ny+1, Nx+1);
    u_actual = zeros(Ny+1, Nx+1);
    v_actual = zeros(Ny+1, Nx+1);
    T_actual = zeros(Ny+1, Nx+1) ;       
    error1 =100;
    error2 =100;
    a=1.5;   %Over relaxation Factor assumed
    
    while (error2>accuracy2)||(error1>accuracy1) 
        psi_old = psi;
        w_old   = w;
        T_old   = T;
  % 3 % Boundary conditions
        % Left  and right wall
        for i=2:Ny
            psi(i,1) = 0;
            u(i,1) = 0;
            v(i,1) = 0;
            psi(i,Nx+1) = 0;
            u(i,Nx+1)=0;
            v(i,Nx+1)=0;
            T(i,1)=1;
            T(i,Nx+1)=0;
            
        end
        %  Top  and bottom wall
        for j=1:Nx+1
            psi(1,j) = 0;
            psi(Ny+1,j) = 0;
            u(1,j) = 0;
            u(Ny+1,j) = 0;
            v(1,j) = 0;
            v(Ny+1,j) = 0;
        end
  % 4 % Stream function at interior points
        for i= 2:Ny
            for j=2:Nx
                psi(i,j) = ((psi(i,j-1) + psi(i,j+1)) + b^2 * ( psi(i+1,j) + psi(i-1,j) ) + w(i,j) * dx^2 )/(2*(1 + b^2));
            end
        end
  % 5 %Vorticity boundary condition
        % Left wall
        for i=2:Ny
            j=1;
            w(i,j) = -2*( psi(i,2) - psi(i,1) )/dx^2;
        end
        % Right wall
        for i=2:Ny
            j=Nx+1;
            w(i,j) = -2*( psi(i,Nx) - psi(i,Nx+1) )/dx^2;
        end
        % Bottom wall
        for j=1:Nx+1
             i=1;
            w(i,j) = -2*( psi(2,j) - psi(1,j) )/dy^2;
        end
        % Top wall
        for j=1:Nx+1
            i=Ny+1;
            w(i,j) = -2*( psi(Ny,j)- psi(Ny+1,j) )/dy^2;
        end
        
 % 6  % Velocity at the interior points
        for i= 2:Ny
            for j=2:Nx
                u(i,j)  = (psi(i+1,j) - psi(i-1,j))/(2* dy);
                v(i,j)  = -(psi(i,j+1) - psi(i,j-1))/(2* dx);
            end
        end
        
  % 7 %Temperature at the interior points 
        for i=1:Ny+1
            for j=2:Nx
                if i==1
                    LHS     = ( u(i,j) * (T(i,j+1) - T(i,j-1))/(2*dx))/(2*dy);
                    RHS     = ( b^2 * ( T(i+1,j) + T(i+1,j) ) +  T(i,j-1) + T(i,j+1) ) /(dx^2);
                    T(i,j)  = dx^2 * ( RHS - LHS )/(2* (b^2 + 1));
                elseif i==Ny+1
                    LHS     = ( u(i,j) * (T(i,j+1) - T(i,j-1))/(2*dx))/(2*dy);
                    RHS     = ( b^2 * ( T(i-1,j) + T(i-1,j) ) +  T(i,j-1) + T(i,j+1) ) /(dx^2);
                    T(i,j)  = dx^2 * ( RHS - LHS )/(2* (b^2 + 1));
                else
                    LHS     = ( u(i,j) * (T(i,j+1) - T(i,j-1))/(2*dx) + v(i,j) * (T(i+1,j) - T(i-1,j))/(2*dy) );
                    RHS     = ( b^2 * ( T(i+1,j) + T(i-1,j) ) +  T(i,j-1) + T(i,j+1) ) /(dx^2);
                    T(i,j)  = dx^2 * ( RHS - LHS )/(2* (b^2 + 1));
                end
            end
        end
                    
                    
                    
        
  % 8 %Vorticity at the interior points
        for i= 2:Ny
            for j=2:Nx
                LHS     = ( u(i,j) * (w(i,j+1) - w(i,j-1))/(2*dx) + v(i,j) * (w(i+1,j) - w(i-1,j))/(2*dy) );
                RHS     = (( b^2 * ( w(i+1,j) + w(i-1,j) ) +  w(i,j-1) + w(i,j+1) )*Pr /(dx^2))+(Ra*Pr*(T(i,j+1)-T(i,j-1))/(2*dx));
                w(i,j)  = dx^2 * ( RHS - LHS )/(2*Pr*(b^2 + 1));
            end
        end
        
        w     = w_old *(1-a) + w*a;    %vorticity by using Over relaxation
        psi   = psi_old *(1-a) + psi*a; %Psi by using Over relaxation
        %T     = T_old*(1-a) + T*a;
        error1 =0;
        error2 =0;
        error1 = (norm(psi - psi_old)/norm(psi))*100;
        error2 = (norm(T - T_old)/norm(T))*100;
    end
   
           
        
    for i=1:Nx+1
        for j=1:Ny+1
            psi_actual((Ny+2-i),j)=psi(i,j);
            w_actual((Ny+2-i),j)=w(i,j);
            u_actual((Ny+2-i),j)=u(i,j);
            v_actual((Ny+2-i),j)=v(i,j);
            T_actual((Ny+2-i),j)=T(i,j);
        end
    end
            
        x=[0:dx:lx];
        y=[0:dy:ly];
%The figures are plotted 
figure(1)
contour(x,y,(flipud(psi_actual)));
title('stream function for buoyancy driven cavity ')
figure(2)
contour(x,y,(flipud(w_actual)),200);
title('Vorticity for buoyancy driven cavity')
figure(3)
contour(x,y,(flipud(u_actual)),100);
title('Horizontal velocity for buoyancy driven cavity')
figure(4)
contour(x,y,(flipud(v_actual)),100);
title('vertical velocity for bouyancy driven cavity ')
figure(5)
contour(x,y,(flipud(T_actual)),200);
title('Temperature for buoyancy driven cavity ')
    for i=1:Ny+1
        Q(i)=(T_actual(2,i)-T_actual(1,i))/dx;
    end
    figure(6)
    plot(Q,y);
    title('gradient at x=0')
    
    figure(7)
    plot(u(:,16),y);
    title('Horizontal velocity at x=0.05')
    
    figure(8)
    plot(y,v(16,:));
    title('Vertical velocity at x=0.05')
    
    
    
   
