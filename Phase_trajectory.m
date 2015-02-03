clear all
close all
clc

% 1) y' < constr
% 2) y' > constr
% 3) constr1 < y' < constr2

Y_0 =   [0.1  0.3  -0.2];
Y_end = [2.5  0.5  -0.5];

N = 1000;

% Строится фазовая кривая вида
% Psi = c0 + c1*(y - y0) + c2*(y - y0)^2 + c3*(y - y0)^3 + d*(y - y0)^2*(y - yend)^2,
% соединяющая
% начальное положение Y_0 и конечное положение Y_end
delta = Y_end(1) - Y_0(1);

c0 = Y_0(2);
c1 = Y_0(3)/Y_0(2);

Ma = [delta^2     delta^3;
      2*delta     3*delta^2];
  
Mb = [Y_end(2) - c1*delta - c0;
      Y_end(3)/Y_end(2) - c1];
Mc = inv(Ma)*Mb;

c2 = Mc(1);
c3 = Mc(2);

Expr = c0 + c1*delta + c2*delta^2 + c3*delta^3;

Expr_1 = c1 + 2*c2*delta + 3*c3*delta^2;

d = 0;

figure(1);
hold on; grid on;
title('Psi(y)');
xlabel('y');
ylabel('dy / d\tau');

dy = (Y_end(1) - Y_0(1)) / N;
y=Y_0(1):dy:Y_end(1);

Psi = c0 + c1*(y - Y_0(1)) + c2*(y - Y_0(1)).^2 + c3*(y - Y_0(1)).^3;
dPsi = c1 + 2*c2*(y - Y_0(1)) + 3*c3*(y - Y_0(1)).^2;

plot(y,Psi);


d_min = -100;
d_max = 100;
hd = 0.01;

N_shift = 30;

cond = false;




%% Убираю выход за ограничение dPsi / dy < 0.5

global dy_constr1
dy_constr1 = 0.5;

N_shift = 50;

cond = false;



% while not(cond)
    flag = 0;
    % Ищем границы интервала, на котором нарушается ограничение
    for i=1:(N+1)
        % ------- 1) y' < constr ------------------------------------------
        if ((flag == 0)&&(dPsi(i) >= dy_constr1))
        % ------- 2) y' > constr ------------------------------------------
        % if ((flag == 0)&&(Psi(i) <= constr1))
        % -----------------------------------------------------------------
            if ((i - N_shift) < 1)
                y_left = Y_0;
            else
                y_left(1) = Y_0(1) + (i-1 - N_shift)*dy;  % отступили N_shift шагов назад от критической точки
                y_left(2) = Psi(i - N_shift);
                y_left(3) = (c1 + 2*c2*(y_left(1) - Y_0(1)) + 3*c3*(y_left(1) - Y_0(1))^2)*y_left(2);
            end
            flag = 1;
        end
        % -------1)  y' < constr ------------------------------------------
        if ((flag==1)&&(dPsi(i) < dy_constr1))
        % -------2) y' > constr -------------------------------------------
        % if ((flag==1)&&(Psi(i) > constr1))
        % -----------------------------------------------------------------
        if ((i + N_shift) > N + 1)
            y_right = Y_end;
        else
            y_right(1) = Y_0(1) + (i-1 + N_shift)*dy;
            y_right(2) = Psi(i + N_shift);
            y_right(3) = (c1 + 2*c2*(y_right(1) - Y_0(1)) + 3*c3*(y_right(1) - Y_0(1))^2)*y_right(2);
        end
            break;
        end
    end
    
    y_loc = y_left(1):dy:y_right(1);
    d = -0.1;
    tau = (y_loc - y_left(1)) / (y_right(1) - y_left(1));     % tau = \tilde{y}
    Psi_1 = c0 + c1*(y_loc - Y_0(1)) + c2*(y_loc - Y_0(1)).^2 + c3*(y_loc - Y_0(1)).^3 + d*tau.^2.*(3 - 2*tau);
    dPsi_1 = c1 + 2*c2*(y_loc - Y_0(1)) + 3*c3*(y_loc - Y_0(1)).^2 + d*(6*tau - 6*tau.^2);
    
    plot(y_loc,Psi_1,'g','LineWidth',2);
    
    Y_2_end = [y_right(1) Psi_1(end) dPsi_1(end)*Psi_1(end)];
    
    replace_part = curve_synthesis(Y_2_end, Y_end, dy, 0);
    
    Psi_2 = replace_part(1,:);
    dPsi_2 = replace_part(2,:);
    
    y_2 = y_right(1):dy:Y_end(1);
    plot(y_2, Psi_2, 'm');
   

%%



figure(2);
hold on; grid on;

title('dPsi(y) / dy');
plot(y,dPsi,'b');
plot(y_loc,dPsi_1,'g');
plot(y_2, dPsi_2, 'm');
xlabel('y');
ylabel('dy / d\tau');

%%

global constr2
constr2 = 0.6;

plot_data = upper_constraint_Psi(Y_2_end, Y_end, dy);
y_3 = plot_data(1,:);
Psi_3 = plot_data(2,:);
dPsi_3 = plot_data(3,:);

set(0,'CurrentFigure',1);
plot(y_3, Psi_3, 'r','LineWidth',2);
set(0,'CurrentFigure',2);
plot(y_3, dPsi_3, 'r');