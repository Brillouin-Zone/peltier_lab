% PELTIER-EFFECT
    fname = 'C:\von_Server\ETH\BSc Physics\5\Praktikum 3\Peltier_effect2\plots';
    
% GIVEN VALUES:
R1 = 1200; % Ohm
   err_R1 = 120; % ohm; 10 percent
R2 = 389*10^-3; % ohm
   err_R2 = 38.9*10^-3; % ohm; 10 percent
P1 = 200; % ohm
   err_P1 = 50; %ohm; 25percent
I = [4, 8, 12, 16, -4, -8, -12, -16]; %A
I_positiv = [4, 8, 12, 16]; % A
Tb = [303.15, 323.15, 353.15, 383.15]; % K

    % copper; we do not know the error
    l_1 = 50*10^-3; % m; length
    r_1 = 2.03/2 * 10^-3; % m; radius
    F_1 = (r_1)^2 *pi; % surface in m^2
    rho_1 = 0* 10^-8; % omega m
        % interpolation of rho_1:
        x_rho_1 = [293:1:353]; % temperature in K
        y_rho_1 = linspace(1.68, 2.06, 61);
        x_T = [30, 50, 80, 110]+273.15;
        AF = polyfit(x_rho_1, y_rho_1, 1);
        rho_1 = polyval(AF, x_T) % = [1.7443    1.8710    2.0610    2.2510]
    
    % constantan
    l_2 = 50*10^-3; % m; length
    r_2 = 7.04/2 * 10^-3; % m; radius
    F_2 = (r_2)^2 *pi; % surface in m^2
    rho_2 = 44* 10^-8; % omega m

% INTERPOLATION OF k_1 AND k_2
    % k_1; copper
    % k_2; constantan
    x1 = [250, 400]; % Kelvin
    y1 = [401, 391]; % W m^-1 K^-1
    x2 = [275, 400]; % kelvin
    y2 = [21.9, 26.6]; % W m^-1 K^-1
   %{
    figure
    plot (x1, y1, 'm-', x2, y2, 'r-');
    xlabel('temperature [K]');
    ylabel('\kappa');
    title('\kappa for copper and constantan');
   
    % linear fitlines:   copper: y1 = -0.066667*x1 + 417.67         constantan: y2 = 0.0376*x2 + 11.56
    %}
    kappa_copper = -0.066667*Tb(:) + 417.67; % we assume T_b to be precise, such that these kappas have no error
    kappa_constantan = 0.0376*Tb(:) + 11.56;

% MEASUREMENTS:
   % errors
     err_V_I = 10^-5; % V; last digit
     err_I = 20/(50*10^-3) * err_V_I; %=0.0040 A; I = V_I *20/(50*10^-3)-> err_I = 20/(50*10^-3) * err_V_I
     err_I_T = 8.5*10^-6; % A; ein Strich auf Galvanometer-Skala
     err_I_T_array = 8.5*10^-6*ones(4,1);
     % der Fehler von V_T ist abhängig von den verschiedenen Strömen I_T_??
     err_V_P = 10^-5; % V; last digit 
    
   % T_b = 303.15K = 30°C and I = [4, 8, 12, 16, -4, -8, -12, -16] A
     V_I_30 = [9.9, 20.0, 30.2, 40.5, -9.8, -19.9, -30.1, -40.3].*10^-3; % V;  wegen ungenauigkeit des current-supply dieses V in I umrechnen: *20A / (50*10^-3V)
     I_30 = V_I_30 .* 20 /(50*10^-3); % A
     I_T_30 = [172.6, 296.5, 499.2, 697.0, -34.4, -99.6, -111.4, -93.4].*10^-6; % A; read from A1-galvanometer
     V_T_30 = I_T_30 .* R2; % V
       err_V_T_30 = sqrt((R2 * err_I_T)^2 + (I_T_30(:) .* err_R2).^2); % V; =1.0e-04 *[0.0671, 0.1153, 0.1942, 0.2711, 0.0134, 0.0387, 0.0433, 0.0363]
     V_P_30 = [3.33, 6.69, 10.11, 13.51, -3.25, -6.57, -9.91, -13.2].*10^-3; % V    
                
    % T_b = 323.15K = 50°C and I = [4, 8, 12, 16, -4, -8, -12, -16] A
      V_I_50 = [9.9, 20.1, 30.3, 40.5, -9.98, -20.19, -29.1, -40.56].*10^-3; % V;  wegen ungenauigkeit des current-supply dieses V in I umrechnen: *20A / (50*10^-3V)
      I_50 = V_I_50 .* 20 /(50*10^-3); % A
      I_T_50 = [200.0, 324.1, 492.8, 901.0, -68.5, -83.1, -81.4, -108.5].*10^-6; % A; read from A1-galvanometer
      V_T_50 = I_T_50 .* R2; % V
        err_V_T_50 = sqrt((R2 * err_I_T)^2 + (I_T_50(:) .* err_R2).^2); % V; =[0.0778, 0.1261, 0.1917, 0.3505, 0.0266, 0.0323, 0.0317, 0.0422]*10^-4
      V_P_50 = [3.44, 6.93, 10.54, 13.96, -3.37, -6.81, -9.74, -13.65].*10^-3; % V        
  %{ 
    % T_b = 353.15K = 80°C and I = [4, 8, 12, 16, -4, -8, -12, -16] A
      V_I_80 = [10.01, 20.23, 30.19, 40.55, -10.07, -19.96, -30.32, -40.41].*10^-3; % V;  wegen ungenauigkeit des current-supply dieses V in I umrechnen: *20A / (50*10^-3V)
      I_80 = V_I_80 .* 20 /(50*10^-3); % A
      I_T_80 = [333.7, 396.3, 441.4, 464.5, -217.0, -37.9, -99.8, -380.5].*10^-6; % A; read from A1-galvanometer
      V_T_80 = I_T_80 .* R2; % V
        err_V_T_80 = sqrt((R2 * err_I_T)^2 + (I_T_80(:) .* err_R2).^2); % V; =1.0e-04 *[0.0844, 0.0147, 0.0388, 0.1480, 0.1298, 0.1542, 0.1717, 0.1807]
      V_P_80 = [4.06, 7.79, 11.22, 14.92, -3.14, -6.67, -10.27, -13.82].*10^-3; % V

    % T_b = 383.15K = 110°C and I = [4, 8, 12, 16, -4, -8, -12, -16] A
      V_I_110 = [9.84, 20.27, 30.3, 40.56, -9.93, -20.06, -30.22, -40.62].*10^-3; % V;  wegen ungenauigkeit des current-supply dieses V in I umrechnen: *20A / (50*10^-3V)
      I_110 = V_I_110 .* 20 /(50*10^-3); % A
      I_T_110 = [548.2, 545.5, 514.3, 690.0, -204.3, -308.1, -372.5, -415.4].*10^-6; % A; read from A1-galvanometer
      V_T_110 = I_T_110 .* R2; % V
         err_V_T_110 = sqrt((R2 * err_I_T)^2 + (I_T_110(:) .* err_R2).^2) % V
      V_P_110 = [3.81, 7.61, 11.57, 15.43, -3.28, -7.0, -10.72, -14.37].*10^-3; % V
%}
      
    % T_b = 353.15K = 80°C and I = [4, 8, 12, 16, -4, -8, -12, -16] A
      V_I_80 = [10.29, 20.17, 29.93, 40.15, -9.46, -20.06, -30.14, -40.47].*10^-3; % V;  wegen ungenauigkeit des current-supply dieses V in I umrechnen: *20A / (50*10^-3V)
      I_80 = V_I_80 .* 20 /(50*10^-3); % A
      I_T_80 = [201.8, 319.7, 457.7, 663.0, -115.7, -202.9, -255.9, -276.9].*10^-6; % A; read from A1-galvanometer
      V_T_80 = I_T_80 .* R2; % V
        err_V_T_80 = sqrt((R2 * err_I_T)^2 + (I_T_80(:) .* err_R2).^2); % V; =1.0e-04 *[0.0844, 0.0147, 0.0388, 0.1480, 0.1298, 0.1542, 0.1717, 0.1807]
      V_P_80 = [3.83, 7.40, 10.29, 14.60, -3.11, -6.84, -10.37, -13.92].*10^-3; % V

    % T_b = 383.15K = 110°C and I = [4, 8, 12, 16, -4, -8, -12, -16] A
      V_I_110 = [9.83, 20.07, 30.19, 40.41, -9.83, -20.29, -30.13, -40.75].*10^-3; % V;  wegen ungenauigkeit des current-supply dieses V in I umrechnen: *20A / (50*10^-3V)
      I_110 = V_I_110 .* 20 /(50*10^-3); % A
      I_T_110 = [510.9, 619.0, 803.0, 1046.0, 316.5, 128.0, 52.6, 8.4].*10^-6; % A; read from A1-galvanometer
      V_T_110 = I_T_110 .* R2; % V
         err_V_T_110 = sqrt((R2 * err_I_T)^2 + (I_T_110(:) .* err_R2).^2); % V
      V_P_110 = [3.77, 7.56, 11.34, 15.15, -3.50, -7.39, -10.97, -14.77].*10^-3; % V
      
      
      
% TABLE 2: calculation of DeltaT^{\pm} with reference table:
voltage = [-0.392, -0.353, -0.314, -0.275, -0.236, -0.197, -0.157, -0.118, -0.079, -0.039, 0.000, 0.039, 0.079, 0.119, 0.158, 0.198, 0.238, 0.277, 0.317, 0.357, 0.397, 0.437, 0.477, 0.517, 0.557, 0.597, 0.637, 0.677, 0.718, 0.758, 0.798].*10^-3; % V
temperature = [-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]; % is the temperature difference
   %{ 
    % interpolation of the voltages:
        xq = -10:0.1:20;
        figure;
        vp = interp1(temperature, voltage, xq);
        plot(temperature,voltage,'o',xq,vp,':.');
        xlim([-10 20]);
        xlabel('temperature differences');
        ylabel('voltages in mV');
        title('Linear Interpolation of table 2'); 
     %} 
    coeffs = polyfit(voltage, temperature, 1); % linear interpolation: temperature = m*voltage +q
        lookfor_deltaT_30 = polyval(coeffs, V_T_30);
        lookfor_deltaT_50 = polyval(coeffs, V_T_50);
        lookfor_deltaT_80 = polyval(coeffs, V_T_80);
        lookfor_deltaT_110 = polyval(coeffs, V_T_110);
        
        %calculation of DeltaT^{\pm} with formula (6) and (7):
        DeltaTcopper_formula =0;
        DeltaTconstantan_formula = 0;
        for i = 1:4
            DeltaTcopper_formula(i) = rho_1(i) / (2*kappa_copper(i) * (F_1)^2)* (l_1)^2 * I_positiv(i)^2 + l_1 * 1 %%%%%%%%% hier am schluss muss die ableitung bei 0 stehen
            DeltaTconstantan_formula(i) = rho_2 / (2*kappa_constantan(i) * (F_2)^2)* (l_2)^2 * I_positiv(i)^2 + l_2 * 1 %%%%%%%%% hier am schluss muss die ableitung bei 0 stehen
        end
        
        
% CALCULATION OF THE PELTIER-COEFFICIENTS:
% formula 10: PI_I = ((k_1 * F_1)/(l_1) + (k_2 * F_2)/(l_2)) * (Delta_T_p - Delta_T_m)/(2* I)
    PI_12I_30 = zeros(4,1); % 4 Zeilen, 1 Spalte
    PI_12I_50 = zeros(4,1);
    PI_12I_80 = zeros(4,1);
    PI_12I_110 = zeros(4,1);
    for i = 1:4
        PI_12I_30(i) = ((kappa_copper(i) * F_1)/(l_1) + (kappa_constantan(i) * F_2)/(l_2)) * (lookfor_deltaT_30(i) - lookfor_deltaT_30(i+4))/(2* I_30(i));
        PI_12I_50(i) = ((kappa_copper(i) * F_1)/(l_1) + (kappa_constantan(i) * F_2)/(l_2)) * (lookfor_deltaT_50(i) - lookfor_deltaT_50(i+4))/(2* I_50(i));
        PI_12I_80(i) = ((kappa_copper(i) * F_1)/(l_1) + (kappa_constantan(i) * F_2)/(l_2)) * (lookfor_deltaT_80(i) - lookfor_deltaT_80(i+4))/(2* I_80(i));
        PI_12I_110(i) = ((kappa_copper(i) * F_1)/(l_1) + (kappa_constantan(i) * F_2)/(l_2)) * (lookfor_deltaT_110(i) - lookfor_deltaT_110(i+4))/(2* I_110(i));
    end
    PI_12I_30; 
    PI_12I_50; 
    PI_12I_80; 
    PI_12I_110; 
    
    
% formula 12: PI_V = V_P / 2 * ((Delta_T_p - Delta_T_m)/(Delta_T_p + Delta_T_m))
   % for 30, 50_1, 50_2, 80, 110 °C
    PI_12V_30 = zeros(4,1); % 4 Zeilen, 1 Spalte
    PI_12V_50 = zeros(4,1);
    PI_12V_80 = zeros(4,1);
    PI_12V_110 = zeros(4,1);
    for i = 1:4;
        PI_12V_30(i) = V_P_30(i) / 2 * ((lookfor_deltaT_30(i) - lookfor_deltaT_30(i+4)))/((lookfor_deltaT_30(i) + lookfor_deltaT_30(i+4)));
        PI_12V_50(i) = V_P_50(i) / 2 * ((lookfor_deltaT_50(i) - lookfor_deltaT_50(i+4)))/((lookfor_deltaT_50(i) + lookfor_deltaT_50(i+4)));
        PI_12V_80(i) = V_P_80(i) / 2 * ((lookfor_deltaT_80(i) - lookfor_deltaT_80(i+4)))/((lookfor_deltaT_80(i) + lookfor_deltaT_80(i+4)));
        PI_12V_110(i) = V_P_110(i) / 2 * ((lookfor_deltaT_110(i) - lookfor_deltaT_110(i+4)))/((lookfor_deltaT_110(i) + lookfor_deltaT_110(i+4)));
    end
    PI_12V_30; 
    PI_12V_50; 
    PI_12V_80;        
    PI_12V_110; 
     
   
% ERROR CALCULATION USING GAUSS ERROR PROPAGATION:    
  % Fehler von I, V_I, V_P, V_T und I_T weiter oben in MEASUREMENTS
  % weitere fehlerbehaftete Grössen: lookfor_delta_T_??, PI_12V_??,
  % PI_12I_??
        % Fehler von lookfor_delta_T_??: arrays mit 8 einträgen
            err_lookfor_delta_T_30 = err_V_T_30(:) .* coeffs(1); % nach Gauss: y=mx+q -> erry = m*errx
            err_lookfor_delta_T_50 = err_V_T_50(:) .* coeffs(1);
            err_lookfor_delta_T_80 = err_V_T_80(:) .* coeffs(1);
            err_lookfor_delta_T_110 = err_V_T_110(:) .* coeffs(1);
        % Fehler von lookfor_delta_T_??_plus - lookfor_delta_T_??_minus:
            err_delta_T_30_minus = sqrt(err_lookfor_delta_T_30(1:4).^2 + err_lookfor_delta_T_30(5:8).^2);
            err_delta_T_50_minus = sqrt(err_lookfor_delta_T_50(1:4).^2 + err_lookfor_delta_T_50(5:8).^2);
            err_delta_T_80_minus = sqrt(err_lookfor_delta_T_80(1:4).^2 + err_lookfor_delta_T_80(5:8).^2);
            err_delta_T_110_minus = sqrt(err_lookfor_delta_T_110(1:4).^2 + err_lookfor_delta_T_110(5:8).^2);
            % NOTE: SUMME UND DIFFERENZ DER DELTA_T HABEN GLEICHEN FEHLER
            
        % Fehler von PI_12V_??: arrays mit 4 einträgen
           err_PI_12V_30 = zeros(4,1);
           err_PI_12V_50 = zeros(4,1);
           err_PI_12V_80 = zeros(4,1);
           err_PI_12V_110 = zeros(4,1);
            for i = 1:4
                err_PI_12V_30(i) = sqrt((((lookfor_deltaT_30(i) - lookfor_deltaT_30(i+4)))/(2*(lookfor_deltaT_30(i) + lookfor_deltaT_30(i+4)))* err_V_P)^2 + (err_lookfor_delta_T_30(i) * V_P_30(i) * lookfor_deltaT_30(i+4)/(lookfor_deltaT_30(i) + lookfor_deltaT_30(i+4))^2)^2 + (err_lookfor_delta_T_30(i+4) * V_P_30(i) * lookfor_deltaT_30(i)/(lookfor_deltaT_30(i) + lookfor_deltaT_30(i+4))^2)^2);
                err_PI_12V_50(i) = sqrt((((lookfor_deltaT_50(i) - lookfor_deltaT_50(i+4)))/(2*(lookfor_deltaT_50(i) + lookfor_deltaT_50(i+4)))* err_V_P)^2 + (err_lookfor_delta_T_50(i) * V_P_50(i) * lookfor_deltaT_50(i+4)/(lookfor_deltaT_50(i) + lookfor_deltaT_50(i+4))^2)^2 + (err_lookfor_delta_T_50(i+4) * V_P_50(i) * lookfor_deltaT_50(i)/(lookfor_deltaT_50(i) + lookfor_deltaT_50(i+4))^2)^2);
                err_PI_12V_80(i) = sqrt((((lookfor_deltaT_80(i) - lookfor_deltaT_80(i+4)))/(2*(lookfor_deltaT_80(i) + lookfor_deltaT_80(i+4)))* err_V_P)^2 + (err_lookfor_delta_T_80(i) * V_P_80(i) * lookfor_deltaT_80(i+4)/(lookfor_deltaT_80(i) + lookfor_deltaT_80(i+4))^2)^2 + (err_lookfor_delta_T_80(i+4) * V_P_80(i) * lookfor_deltaT_80(i)/(lookfor_deltaT_80(i) + lookfor_deltaT_80(i+4))^2)^2);
                err_PI_12V_110(i) = sqrt((((lookfor_deltaT_110(i) - lookfor_deltaT_110(i+4)))/(2*(lookfor_deltaT_110(i) + lookfor_deltaT_110(i+4)))* err_V_P)^2 + (err_lookfor_delta_T_110(i) * V_P_110(i) * lookfor_deltaT_110(i+4)/(lookfor_deltaT_110(i) + lookfor_deltaT_110(i+4))^2)^2 + (err_lookfor_delta_T_110(i+4) * V_P_110(i) * lookfor_deltaT_110(i)/(lookfor_deltaT_110(i) + lookfor_deltaT_110(i+4))^2)^2);
            end
            err_PI_12V_30;
            err_PI_12V_50;
            err_PI_12V_80;
            err_PI_12V_110;
           
        % Fehler von PI_12I_?? = ((k_1 * F_1)/(l_1) + (k_2 * F_2)/(l_2)) * (Delta_T_p - Delta_T_m)/(2* I)
           a = ((kappa_copper(:) * F_1)/(l_1) + (kappa_constantan(:) * F_2)/(l_2))/2; % without error
           err_PI_12I_30 = zeros(4,1);
           err_PI_12I_50 = zeros(4,1);
           err_PI_12I_80 = zeros(4,1);
           err_PI_12I_110 = zeros(4,1);
           for i = 1:4
               err_PI_12I_30(i) = sqrt((a(i)/I_30(i) * err_lookfor_delta_T_30(i))^2 + (-a(i)/I_30(i)*err_lookfor_delta_T_30(i+4))^2 + (a(i)*(lookfor_deltaT_30(i+4) - lookfor_deltaT_30(i))/(I_30(i))^2 * err_I)^2);
               err_PI_12I_50(i) = sqrt((a(i)/I_50(i) * err_lookfor_delta_T_50(i))^2 + (-a(i)/I_50(i)*err_lookfor_delta_T_50(i+4))^2 + (a(i)*(lookfor_deltaT_50(i+4) - lookfor_deltaT_50(i))/(I_50(i))^2 * err_I)^2);
               err_PI_12I_80(i) = sqrt((a(i)/I_80(i) * err_lookfor_delta_T_80(i))^2 + (-a(i)/I_80(i)*err_lookfor_delta_T_80(i+4))^2 + (a(i)*(lookfor_deltaT_80(i+4) - lookfor_deltaT_80(i))/(I_80(i))^2 * err_I)^2);
               err_PI_12I_110(i) = sqrt((a(i)/I_110(i) * err_lookfor_delta_T_110(i))^2 + (-a(i)/I_110(i)*err_lookfor_delta_T_110(i+4))^2 + (a(i)*(lookfor_deltaT_110(i+4) - lookfor_deltaT_110(i))/(I_110(i))^2 * err_I)^2);
           end
           err_PI_12I_30;
           err_PI_12I_50;
           err_PI_12I_80;
           err_PI_12I_110;

% FINAL VERSION OF TEMPERATURE PLOTS:
%{
Interval1 = [-20:0.01:20];
Interval2 = [0:0.01:20];
%1. Delta T^{\pm} vs I: quadratisch wegen gleichung 6, 7
    figure
    errorbar(I_30(:), lookfor_deltaT_30(:), err_lookfor_delta_T_30(:), 'bo');
    A1 = polyfit(I_30([8:-1:5 1:4]), lookfor_deltaT_30([8:-1:5 1:4]),2);
    A2 = polyval(A1, Interval1);
    hold on
    plot(Interval1, A2, 'b-');
    hold on
    errorbar(I_50(:), lookfor_deltaT_50(:), err_lookfor_delta_T_50(:), 'mo');
    A3 = polyfit(I_50([8:-1:5 1:4]), lookfor_deltaT_50([8:-1:5 1:4]),2);
    A4 = polyval(A3, Interval1);
    hold on
    plot(Interval1, A4, 'm-');
    hold on
    errorbar(I_80(:), lookfor_deltaT_80(:), err_lookfor_delta_T_80(:), 'go');
    A5 = polyfit(I_80([8:-1:5 1:4]), lookfor_deltaT_80([8:-1:5 1:4]),2);
    A6 = polyval(A5, Interval1);
    hold on
    plot(Interval1, A6, 'g-');
    hold on
    errorbar(I_110(:), lookfor_deltaT_110(:), err_lookfor_delta_T_110(:), 'co');
    A7 = polyfit(I_110([8:-1:5 1:4]), lookfor_deltaT_110([8:-1:5 1:4]),2);
    A8 = polyval(A7, Interval1);
    hold on
    plot(Interval1, A8, 'c-');
    grid on
    xlabel('current I [A]')
    ylabel('\Delta T^{\pm} [K]')
    title('\Delta T^{\pm} vs. current I [A]')
    h = legend('30^\circ C', 'quadratic fit', '50^\circ C', 'quadratic fit', '80^\circ C', 'quadratic fit', '110^\circ C', 'quadratic fit', 'Location','northwest');
    %legend('boxoff');
    saveas(gcf, fullfile(fname, 'peltier_effect_deltaTpm.eps'), 'epsc');      

%2. (Delta T^{+} - Delta T^{-}) vs I : linear wegen  gleichung 12, und weil
%+ qaudratisch ist
    figure
    errorbar(I_30(1:4), lookfor_deltaT_30(1:4)-lookfor_deltaT_30(5:8), err_delta_T_30_minus, 'bo');
    B1 = polyfit(I_30(1:4), lookfor_deltaT_30(1:4)-lookfor_deltaT_30(5:8),1); % y = -0.0020x^2 + 0.5116x + 0.0042
    B2 = polyval(B1,Interval2);
    hold on
    plot(Interval2, B2, 'b-');
    hold on
    errorbar(I_50(1:4), lookfor_deltaT_50(1:4)-lookfor_deltaT_50(5:8), err_delta_T_50_minus, 'mo');
    B3 = polyfit(I_50(1:4), lookfor_deltaT_50(1:4)-lookfor_deltaT_50(5:8), 1); % y= 0.0437x^2 -0.3061x + 3.2771
    B4 = polyval(B3, Interval2);
    hold on
    plot(Interval2, B4, 'm-');
    hold on
    errorbar(I_80(1:4), lookfor_deltaT_80(1:4)-lookfor_deltaT_80(5:8), err_delta_T_80_minus, 'go');
    B5 = polyfit(I_80(1:4), lookfor_deltaT_80(1:4)-lookfor_deltaT_80(5:8), 1); % y = 0.0022x^2 + 0.4633x + 1.1863
    B6 = polyval(B5, Interval2);
    hold on
    plot(Interval2, B6, 'g-');
    hold on
    errorbar(I_110(1:4), lookfor_deltaT_110(1:4)-lookfor_deltaT_110(5:8), err_delta_T_110_minus, 'co');
    B7 = polyfit(I_110(1:4), lookfor_deltaT_110(1:4)-lookfor_deltaT_110(5:8), 1); % y= -0.0013x^2 * 0.6973x - 0.7872
    B8 = polyval(B7, Interval2);
    hold on
    plot(Interval2, B8, 'c-');
    grid on
    xlabel('current I [A]')
    ylabel('\Delta T^{+} - \Delta T^{-} [K]')
    title('\Delta T^{+} - \Delta T^{-} [K] vs. current I [A]')
    h = legend('30^\circ C', 'linear fit', '50^\circ C', 'linear fit', '80^\circ C', 'linear fit', '110^\circ C', 'linear fit', 'Location','northwest');
    %legend('boxoff');
    saveas(gcf, fullfile(fname, 'peltier_effect_deltaTp-m.eps'), 'epsc');
        
%3. (Delta T^{+} + Delta T^{-}) vs I; quadratisch wegen gleichung 11 
    figure
    errorbar(I_30(1:4), lookfor_deltaT_30(1:4)+lookfor_deltaT_30(5:8), err_delta_T_30_minus, 'bo');
    C1 = polyfit(I_30(1:4), lookfor_deltaT_30(1:4)+lookfor_deltaT_30(5:8),2);
    C2 = polyval(C1, Interval2);
    hold on
    plot(Interval2, C2, 'b-');
    hold on
    errorbar(I_50(1:4), lookfor_deltaT_50(1:4)+lookfor_deltaT_50(5:8), err_delta_T_50_minus, 'mo');
    C3 = polyfit(I_50(1:4), lookfor_deltaT_50(1:4)+lookfor_deltaT_50(5:8), 2);
    C4 = polyval(C3, Interval2);
    hold on
    plot(Interval2, C4, 'm-');
    hold on
    errorbar(I_80(1:4), lookfor_deltaT_80(1:4)+lookfor_deltaT_80(5:8), err_delta_T_80_minus, 'go');
    C5 = polyfit(I_80(1:4), lookfor_deltaT_80(1:4)+lookfor_deltaT_80(5:8),2);
    C6 = polyval(C5, Interval2);
    hold on
    plot(Interval2, C6, 'g-');
    hold on
    errorbar(I_110(1:4), lookfor_deltaT_110(1:4)+lookfor_deltaT_110(5:8), err_delta_T_110_minus, 'co');
    C7 = polyfit(I_110(1:4), lookfor_deltaT_110(1:4)+lookfor_deltaT_110(5:8),2);
    C8 = polyval(C7, Interval2);
    hold on
    plot(Interval2, C8, 'c-');
    grid on
    xlabel('current I [A]')
    ylabel('\Delta T^{+} + \Delta T^{-} [K]')
    title('\Delta T^{+} + \Delta T^{-} [K] vs. current I [A]')
    h = legend('30^\circ C', 'quadratic fit', '50^\circ C', 'quadratic fit', '80^\circ C', 'quadratic fit', '110^\circ C', 'quadratic fit', 'Location','BestOutside');
    %legend('boxoff');
    saveas(gcf, fullfile(fname, 'peltier_effect_deltaTp+m.eps'), 'epsc');
%}
           
% FIRST VERSION OF PELTIER-COEFFICIENT-PLOTS:
%{
% PI_I vs I; all graphs in one plot
    % current-peltier-coefficient
    figure
    plot (I_30(1:4),  PI_12I_30(:), 'bo', I_50(1:4),  PI_12I_50(:), 'mo', I_80(1:4),  PI_12I_80(:), 'go', I_110(1:4),  PI_12I_110(:), 'co');
    xlabel('current I [A]');
    ylabel('Peltier coefficient \Pi_{12}^{(I)}');
    title('Peltier coefficient vs. current I [A]');
    h = legend('30^\circ C', '50^\circ C', '80^\circ C', '110^\circ C');
    saveas(gcf, fullfile(fname, 'peltier_effect_PI12_I.eps'), 'epsc');
    
    % voltage-peltier-coefficient
    figure
    plot (I_30(1:4),  PI_12V_30(:), 'bo', I_50(1:4),  PI_12V_50(:), 'mo', I_80(1:4),  PI_12V_80(:), 'go', I_110(1:4),  PI_12V_110(:), 'co');
    xlabel('current I [A]');
    ylabel('Peltier coefficient \Pi_{12}^{(V)}');
    title('Peltier coefficient vs. current I [A]');
    h = legend('30^\circ C', '50^\circ C', '80^\circ C', '110^\circ C', 'Location','northwest');
    saveas(gcf, fullfile(fname, 'peltier_effect_PI12_V.eps'), 'epsc');
   %}  
           
% FINAL VERSION OF PELTIER-COEFFICIENT-PLOTS:

% PI_I vs I; all graphs in one plot
Interval3 = [2:0.01:18];

    % current-peltier-coefficient
    figure
    errorbar(I_30(1:4),  PI_12I_30(:), err_PI_12I_30(:), 'bo')
    A = polyfit(I_30(1:4).', PI_12I_30(:), 1);
    B = polyval(A, Interval3);
    hold on
    plot(Interval3, B, 'b-');
    hold on
    errorbar(I_50(1:4),  PI_12I_50(:), err_PI_12I_50(:), 'mo')
    C = polyfit(I_50(1:4).', PI_12I_50(:), 1);
    D = polyval(C, Interval3);
    hold on
    plot(Interval3, D, 'm-');
    hold on
    errorbar(I_80(1:4),  PI_12I_80(:), err_PI_12I_80(:), 'go')
    E = polyfit(I_80(1:4).', PI_12I_80(:), 1);
    F = polyval(E, Interval3);
    hold on
    plot(Interval3, F, 'g-');
    hold on
    errorbar(I_110(1:4),  PI_12I_110(:), err_PI_12I_110(:), 'co');
    G = polyfit(I_110(1:4).', PI_12I_110(:), 1);
    K = polyval(G, Interval3);
    hold on
    plot(Interval3, K, 'c-');
    grid on
    xlabel('current I [A]');
    ylabel('Peltier coefficient \Pi_{12}^{(I)} [V]');
    title('Peltier coefficient vs. current I [A]');
    h = legend('30^\circ C', 'fitline', '50^\circ C', 'fitline', '80^\circ C', 'fitline', '110^\circ C', 'fitline', 'Location', 'BestOutside');
    %legend('boxoff');
    saveas(gcf, fullfile(fname, 'peltier_effect_PI12_I.eps'), 'epsc');
    
    % voltage-peltier-coefficient
    figure
    errorbar(I_30(1:4),  PI_12V_30(:), err_PI_12V_30(:), 'bo')
    aa = polyfit(I_30(1:4).', PI_12V_30(:), 1);
    b = polyval(aa, Interval3);
    hold on
    plot(Interval3, b, 'b-');
    hold on
    errorbar(I_50(1:4),  PI_12V_50(:), err_PI_12V_50(:), 'mo')
    c = polyfit(I_50(1:4).', PI_12V_50(:), 1);
    d = polyval(c, Interval3);
    hold on
    plot(Interval3, d, 'm-');
    hold on
    errorbar(I_80(1:4),  PI_12V_80(:), err_PI_12V_80(:), 'go')
    e = polyfit(I_80(1:4).', PI_12V_80(:), 1);
    f = polyval(e, Interval3);
    hold on
    plot(Interval3, f, 'g-');
    hold on
    errorbar(I_110(1:4),  PI_12V_110(:), err_PI_12V_110(:), 'co')
    g = polyfit(I_110(1:4).', PI_12V_110(:), 1);
    k = polyval(g, Interval3);
    hold on
    plot(Interval3, k, 'c-');
    grid on
    xlabel('current I [A]');
    ylabel('Peltier coefficient \Pi_{12}^{(V)} [V]');
    title('Peltier coefficient vs. current I [A]');
    h = legend('30^\circ C', 'fitline', '50^\circ C', 'fitline', '80^\circ C', 'fitline', '110^\circ C', 'fitline', 'Location','northwest');
    %legend('boxoff');
    saveas(gcf, fullfile(fname, 'peltier_effect_PI12_V.eps'), 'epsc');

 % Peltier coefficient vs T   
 PI_literature = zeros(1, 4);
 T = [30, 50, 80, 110] + 273.15;
 for i=1:4
   PI_literature(i) = (4.37184*T(i) + 0.1617*T(i)^2 - 1.84371*10^-4 * T(i)^3 + 1.2244*10^-7 * T(i)^4 - 4.47618*10^-11 * T(i)^5)*10^-6; % to get mV
 end
 PI_literature; % = [11.9685   13.2542   15.2486   17.3119]