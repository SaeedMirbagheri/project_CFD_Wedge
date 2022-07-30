clc 
clear                               
close all
disp('programmer:Seid Saeed Mirbagheri (400126116)')

deg=1;  % 1 for 15 degree
        % 2 for 35 degree
n_N=5;
Fac=zeros(4,n_N);
U_half=zeros(4,n_N);
H=zeros(1,n_N); Error1=zeros(4,n_N,n_N);
%% SET FREESTREAM CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
free_rho      = 1.2;                            % kg/m^3
free_T        = 300;                         % Kelvin
free_P        = 100000;                         % Pa
free_M        = 1.2;                            % Mach
free_aoa      = 0;                              % deg

free_a        = speedsound(free_P,free_rho);    % m/s
free_u        = free_M*free_a*cosd(free_aoa);   % m/s
free_v        = free_M*free_a*sind(free_aoa);   % m/s
free_vel      = [free_u free_v];

% Freestream Primitive State Vector (PSV)
V_free        = [free_rho free_u free_v free_P];

diary('output.txt')
fprintf('Freestream Conditions:\n');
fprintf('Mach:         %10.2f\n',free_M);
fprintf('Flow AOA:     %10.2f deg\n',free_aoa);
fprintf('u Velocity:   %10.2f m/s\n',free_u);
fprintf('v Velocity:   %10.2f m/s\n',free_v);
fprintf('Pressure:     %10.2e Pa\n',free_P);
fprintf('Temperature:  %10.2f K\n',free_T);
fprintf('Density:      %10.2f kg/m^3\n\n',free_rho);

%% SET ITERATION VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteration variables
iterations   = 2000;             % Number of iterations to run
                               
timestep     = 1e-3;            % Timestep for global timestepping
CFL          = 0.3;             % Courant number

m_stage      = 4;               % m-stage time stepping
                                % e.g. 1 for Euler step
                                %     4 for 4th-order RK

% Output variables
fprintf('Iteration Variables:\n');
fprintf('Iterations:   %5d\n',iterations);
fprintf('CFL:          %5.2f\n',CFL);
fprintf('M-Stage:      %5d\n',m_stage);


for ii=1:n_N    

resid={zeros(1,1)}  ; V={zeros(1,1)} ; U={zeros(1,1)};
e_xlen=zeros(1,1) ; e_ylen=zeros(1,1);
n_xlen=zeros(1,1) ; n_ylen=zeros(1,1);
w_xlen=zeros(1,1) ; w_ylen=zeros(1,1);
s_xlen=zeros(1,1) ; s_ylen=zeros(1,1);
xmid=zeros(1,1)   ; ymid=zeros(1,1);
volume=zeros(1,1) ; mov=zeros(1,1);
sE=zeros(1,1)     ; sN=zeros(1,1) ; sW=zeros(1,1) ; sS=zeros(1,1);
nE={zeros(1,1)}     ; nN={zeros(1,1)} ; nW={zeros(1,1)} ; nS={zeros(1,1)};

%% LOAD/GENERATE GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Open the file and assign it an ID
    % C:\Users\King\Desktop\New folder\name file
if deg==1
    if ii==1
    fid = fopen('C:\Users\King\Desktop\data for run matlab\Wedge15_1.x');
    elseif ii==2
    fid = fopen('C:\Users\King\Desktop\data for run matlab\Wedge15_2.x');
    elseif ii==3
    fid = fopen('C:\Users\King\Desktop\data for run matlab\Wedge15_3.x');
    elseif ii==4
    fid = fopen('C:\Users\King\Desktop\data for run matlab\Wedge15_4.x');
    else
    fid = fopen('C:\Users\King\Desktop\data for run matlab\Wedge15_5.x');
    end

else
    if ii==1
    fid = fopen('C:\Users\King\Desktop\data for run matlab\Wedge35_1.x');
    elseif ii==2
    fid = fopen('C:\Users\King\Desktop\data for run matlab\Wedge35_2.x');
    elseif ii==3
    fid = fopen('C:\Users\King\Desktop\data for run matlab\Wedge35_3.x');
    elseif ii==4
    fid = fopen('C:\Users\King\Desktop\data for run matlab\Wedge35_4.x');
    else
    fid = fopen('C:\Users\King\Desktop\data for run matlab\Wedge35_5.x');
    end
end
    
    if fid >= 1
        % Read in file headers
        zones = fscanf(fid, '%d', 1);
        % Code only handles 1 zone
        % Therefore, check for number of zones
        if (zones == 1)
            % Read in number of i,j,k points
            npi = fscanf(fid, '%d', 1);
            npj = fscanf(fid, '%d', 1);
            npk = fscanf(fid, '%d', 1);
       
            % Retrieve i,j,k coordinates
            x = fscanf(fid, '%f', [npi,npj]);
            y = fscanf(fid, '%f', [npi,npj]);
            z = fscanf(fid, '%f', [npi,npj]);
            disp('Grid read successfully');
          
           
           x=x.';
           y=y.';
           nn=npi;
           npi=npj;
           npj=nn;
          
        end
        fclose(fid);
    end

%% GRID METRICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z_width = [0,0,1];                   % Unit vector in z-dir
nci = npi-1;                         % Number of cells (pts-1) in i dir
ncj = npj-1;                         % Number of cells (pts-1) in j dir

for i = 1:nci
    for j = 1:ncj
        % Assemble the face lengths
        e_xlen(i,j) = x(i+1,j+1)-x(i+1,j);
        e_ylen(i,j) = y(i+1,j+1)-y(i+1,j);
        
        n_xlen(i,j) = x(i,j+1)-x(i+1,j+1);
        n_ylen(i,j) = y(i,j+1)-y(i+1,j+1);

        w_xlen(i,j) = x(i,j)-x(i,j+1);
        w_ylen(i,j) = y(i,j)-y(i,j+1);

        s_xlen(i,j) = x(i+1,j)-x(i,j);
        s_ylen(i,j) = y(i+1,j)-y(i,j);
        
        % Compute midpoint of cell (for plotting)
        xmid(i,j) = (x(i,j) + x(i+1,j))/2;
        ymid(i,j) = (y(i,j) + y(i,j+1))/2;
        
        % Compute volume of cell using A.BxC (volume of parallelepiped)
        volume(i,j) = abs(dot(z_width,cross(-1*[s_xlen(i,j),s_ylen(i,j),0],...
                      [e_xlen(i,j),e_ylen(i,j),0])));
        
        % Compute area of cell
        sE(i,j) = sqrt((e_xlen(i,j))^2 + (e_ylen(i,j))^2);
        sN(i,j) = sqrt((n_xlen(i,j))^2 + (n_ylen(i,j))^2);
        sW(i,j) = sqrt((w_xlen(i,j))^2 + (w_ylen(i,j))^2);
        sS(i,j) = sqrt((s_xlen(i,j))^2 + (s_ylen(i,j))^2);
        
        % Compute outward normal of faces (return 3 component vector)
        temp_nE = cross([e_xlen(i,j),e_ylen(i,j), 0]/sE(i,j), z_width);
        temp_nN = cross([n_xlen(i,j),n_ylen(i,j), 0]/sN(i,j), z_width);
        temp_nW = cross([w_xlen(i,j),w_ylen(i,j), 0]/sW(i,j), z_width);
        temp_nS = cross([s_xlen(i,j),s_ylen(i,j), 0]/sS(i,j), z_width);
        
        % Truncate normal vector to 2 components
        nE{i,j} = [temp_nE(1) temp_nE(2)];
        nN{i,j} = [temp_nN(1) temp_nN(2)];
        nW{i,j} = [temp_nW(1) temp_nW(2)];
        nS{i,j} = [temp_nS(1) temp_nS(2)];
        
        
    end 
end


%% INITIALIZE SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resid_i      = 0;               % Iterative residual
resid_0      = 0;               % Step 0 residual
start_iter   = 0;               % Used for multiple runs
end_iter     = 0;               % Used for multiple runs
residReduced = 0;               % If divergence detected
fm = 1;                         % Used for capturing movies

% Combine normals and areas into big cell array which will be passed
% to the function which computes the residual
normals = {nE nN nW nS};
areas   = {sE sN sW sS};

% Initalize variables which will allow for visualization
% i.e. Plot Contours
con_density = zeros([nci,ncj]);
con_uvel    = zeros([nci,ncj]);
con_vvel    = zeros([nci,ncj]);
con_pres    = zeros([nci,ncj]);

% Loop through all cells and init PSV to freestream conditions
% Convert PSV to conservative state vector
% Init residual to 0
for i = 1:nci
    for j = 1:ncj
        V{i,j}     = V_free;
        U{i,j}     = convV_U(V{i,j});
        resid{i,j} = [0 0 0 0];
    end
end



%% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_iter = start_iter + 1;        % Start at iteration 1

% Main loop in time
for iter = start_iter:(end_iter + iterations)
    % Time variable used to measure time/iteration
    ti1 = cputime;
    
    % Initialize iteration residual to 0
    resid_i = 0;
    
    % Save CSV from this timestep (to be used in m-stage)
    U0 = U;
    
    % M-stage timestepping scheme
    for m = 1:m_stage
        % Calculate residual using function calcResid
        % Passes PSV, normals, areas cell array,
        %        freestream PSV, and nci and ncj
        resid = calcResid(V, V_free, normals, areas, nci, ncj);
        
        % Loop through all cells to update solution
        for i = 1:nci
            for j = 1:ncj
                
                
                    vel = [V{i,j}(2) V{i,j}(3)];
                    cell_a = speedsound(V{i,j}(4),V{i,j}(1));
                    dt(1) = CFL * sE(i,j)/(abs(vel(1)*nE{i,j}(1) +...
                            vel(2)*nE{i,j}(2))+cell_a);
                    dt(2) = CFL * sN(i,j)/(abs(vel(1)*nN{i,j}(1) +...
                            vel(2)*nN{i,j}(2))+cell_a);
                    dt(3) = CFL * sW(i,j)/(abs(vel(1)*nW{i,j}(1) +...
                            vel(2)*nW{i,j}(2))+cell_a);
                    dt(4) = CFL * sS(i,j)/(abs(vel(1)*nS{i,j}(1) +...
                            vel(2)*nS{i,j}(2))+cell_a);
                    timestep = min(dt);
            
                
                % Update solution using the saved CSV
                % Multiply by 'alpha' constant 
                U{i,j} = U0{i,j} -...
                         1/(m_stage-(m-1))*timestep/volume(i,j)*resid{i,j};
                
                % Update cell PSV
                V{i,j} = convU_V(U{i,j});
            
                % Update contour arrays used for plotting
                con_density(i,j) = V{i,j}(1);
                con_uvel(i,j)    = V{i,j}(2);
                con_vvel(i,j)    = V{i,j}(3);
                con_pres(i,j)    = V{i,j}(4);

                % Assemble first part of L2 norm for residual
                resid_i = resid_i + (resid{i,j}.^2);
            end
        end
    end   
    
    % Assemble second part of L2 norm for residual
    resid_i = (resid_i).^.5/(nci*ncj);
    
    % Assign normalization value in first interation
    if iter == 1
        resid_0 = resid_i;
         
    end
    
    % Detects extreme divergence (at the point of no return)
    % and shuts down simulation
    if isnan(resid_i/resid_0)
        disp('Solution corrupt.');
        break;
    end
    
    % Detects divergence happening in x-mom resid and cuts CFL in half
        ee=(resid_i(2)/resid_0(2));
        if ee >= (1e1)
            if residReduced == 0
                CFL = CFL/2; 
                notice = sprintf('Divergence detected.  CFL reduced to %5.2f',CFL);
                disp(notice);
                residReduced = residReduced + 1;
            end
        end

    % Computes time/iteration
    ti2 = cputime-ti1;
    Error1(:,iter,ii)=[resid_i(1)/resid_0(1);resid_i(2)/resid_0(2);resid_i(3)/resid_0(3);resid_i(4)/resid_0(4)];
       
  if(max(abs(Error1)))<=1e-3 
      break;
  end
end

start_iter = iter;
end_iter = iter;

%% Mach matrix calculation
 u_final=con_uvel;
 h=mean(mean(e_ylen));

%% PLOT PRIMITIVE STATE VECTOR CONTOURS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates a figure with 4 contour subplots
% Plot 1: Density (kg/m^3) contour
% Plot 2: U Velocity (m/s) contour
% Plot 3: V Velocity (m/s) contour
% Plot 4: Pressure (Pa)


    figure(1);
    subplot(2,3,ii);
    contourf(xmid',ymid',con_density');
    axis image;
    colorbar('peer',gca,'SouthOutside'); 
    title(sprintf('Density (kg/m^3) h= %d ' ,h))

    figure(2);
    subplot(2,3,ii);
    contourf(xmid',ymid',con_uvel');
    axis image;
    colorbar('peer',gca,'SouthOutside'); 
    title(sprintf('U Velocity (m/s) h= %d ' ,h))

    figure(3);
    subplot(2,3,ii);
    contourf(xmid',ymid',con_vvel');
    axis image;
    colorbar('peer',gca,'SouthOutside'); 
    title(sprintf('V Velocity (m/s) h= %d ' ,h))

    figure(4);
    subplot(2,3,ii);
    contourf(xmid',ymid',con_pres');
    axis image;
    colorbar('peer',gca,'SouthOutside'); 
    title(sprintf('Pressure (Pa) h= %d ' ,h))


    % Error 
    ite=1:iterations;
    figure(10);
    subplot(2,3,ii)
    plot(ite,Error1(1,:,ii),'r')
    hold on
    plot(ite,Error1(2,:,ii),'b')
    plot(ite,Error1(3,:,ii),'g')
    plot(ite,Error1(4,:,ii),'k')
    legend('cont. resid','x-mom resid','y-mom resid','energy resid')
    xlabel('Iterations')
    ylabel('Residual Error')
    title(sprintf('Residual Error(Iterations) for h= %d ' ,h))

    yy=0:h:3;
    n_points=length(yy);
    H(ii)=h;
    Fac(:,ii)=[ ncj/2 ;ncj*(5/8); ncj*(3/4);ncj*(7/8)];
for i=1:4
    if deg==1
        U_half(i,ii)=u_final(nci,Fac(i,ii));
    else
        U_half(i,ii)=u_final((npi+1)/2,Fac(i,ii));
    end

end


end

Er=Error(U_half,n_N);
R_e_delta_y=Error_Slope( Er,H,n_N );

%%%%%Draw the relative error value in terms of H
COL=['s','*','d','o'];
for i=1:4
figure(12)
subplot(1,2,1)
loglog(H(1:end-1),abs(Er(i,:)),sprintf(COL(i)),'linewidth',2);
title('Solution convergence )');
ylabel('successive error');
xlabel('h');
hold on
grid on

%%%%Draw the error slope in terms of h
subplot(1,2,2)
loglog(H(1:end-2),abs(R_e_delta_y(i,:)),'linewidth',2);
title('Error slope');
ylabel('Slope');
xlabel('h');
hold on
grid on
end
legend('y=1.5','y=1.875','y=2.25','y=2.625')
