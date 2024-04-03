    %Ch1603
    clear all;
    close all;
    clc;
    tic;
    format long
    global tdata_cases dna_datas initial_cond tforward
   tinitial=182;
    tforwardlast=385;
    tend=(tforwardlast-tinitial)*10+1;
  tforward = tinitial:0.1:tforwardlast;

%    % experimental data for time (days)
tdata_cases = [182; 189; 196; 203; 210; 217; 224;...
               231; 238; 245; 252; 259; 266; 273; 280; 287; 294; 301;...
               315; 329; 343; 357; 371; 385];


 dna_datas = [1715234000; 2964471000; 1715234000; 1822746000;...
      992429500; 1012746000; 1518861000; 1054635000; 133478900;...
      141845400; 450262100; 217086900; 450262100; 276849800;...
     540347700; 732295700; 2470241000; 992429500; 276849800;...
     77230490; 133478900; 133478900; 150736300; 294202800]; 
     
% Initial value of parameter
% params = lambda,        k,           a,   gamma,    beta,     mu,     delta,     c  
 params = [25829482  1.001e-12        110    0.49     2.0      0.019   0.0380     4.485];
 lb = [25830000      1.000*10^(-12)    105    0.34     2        0.019   0.038      4.480];   
 ub = [25830000      1.001*10^(-12)    110    0.50    2.1      0.030   0.060      5.50];   
 % initial condition
%   initial_cond = [2*10^9.9  0 100 0];
  tinitialp=28;
     tforwardlastp=161;
     tendp=(tforwardlastp-tinitialp)*10+1;
   solution_previous_point=readmatrix('Ch1616_28_182_sol.txt');
  initial_cond = [solution_previous_point(tendp,2),solution_previous_point(tendp,3),...
     solution_previous_point(tendp,4),solution_previous_point(tendp,5)];
  [params,fval] = fminsearchbnd(@err_in_data, params, lb, ub, optimset('Display','iter'));

 [~, Y] = ode45(@(t,y) model_1(t,y, params), tforward, initial_cond);

figure(1)
ylim([0 5000000000])
hold on
semilogy(tforward,Y(:,3),'-');
hold on 
semilogy(tdata_cases, dna_datas, 'r.', 'MarkerSize',20);
title('Chimpanzee-1616');
saveas(gcf,'Ch1616_182_385','fig');
saveas(gcf,'Ch1616_182_385','jpg');
 filename = 'Ch1603_182_385.xlsx';
 D=[tforward'  Y(:,3)];
 writematrix(D,'Ch1616_182_385_sol');

 
%  display('Parameters after data fitting:\n');
% fprintf('lambda = %g\n', params(1));
%  fprintf('k = %g\n',  params(2));   
%  fprintf('a = %g\n', params(3));
%  fprintf('gamma = %g\n', params(4));
%  fprintf('beta = %g\n',  params(5));
% %  fprintf('eta = %g\n', params(6));
%  fprintf('mu = %g\n',  params(6));
%  fprintf('delta = %g\n',  params(7));
%  fprintf('c = %g\n',  params(8));
 estimate_para=[params(1) params(2) params(3) params(4) params(5) params(6) params(7) params(8)];
 filename = 'Ch1616_182_385_para.xlsx';
 writematrix(estimate_para,'Ch1616_182_385_para');
toc

 function dy = model_1(t,y,params)
          dy = zeros(4,1);        
          %  Model Parameters
          lambda = params(1);
          k = params(2);
          a = params(3);
          gamma = params(4);
          beta = params(5);
          mu = params(6);
          delta = params(7);
          c = params(8);
        
          % Model equations 
        dy(1) = lambda-mu*y(1)-k*y(1)*y(4);
        dy(2) = k*y(1)*y(4)-delta*y(2);
        dy(3) = a*y(2)+gamma*(1-0.8)*y(3)-0.8*beta*y(3)-delta*y(3);
        dy(4) = 0.8*beta*y(3)-c*y(4);
 end

function error_in_data = err_in_data(k1)
 tinitialp=28;
     tforwardlastp=182;
     tendp=(tforwardlastp-tinitialp)*10+1;
   solution_previous_point=readmatrix('Ch1616_28_182_sol.txt');
initial_cond = [solution_previous_point(tendp,2),solution_previous_point(tendp,3),...
     solution_previous_point(tendp,4),solution_previous_point(tendp,5)];
 tdata=[182; 189; 196; 203; 210; 217; 224;...
               231; 238; 245; 252; 259; 266; 273; 280; 287; 294; 301;...
               315; 329; 343; 357; 371; 385];

qdata = [1715234000; 2964471000; 1715234000; 1822746000;...
      992429500; 1012746000; 1518861000; 1054635000; 133478900;...
      141845400; 450262100; 217086900; 450262100; 276849800;...
     540347700; 732295700; 2470241000; 992429500; 276849800;...
     77230490; 133478900; 133478900; 150736300; 294202800]; 
 
  tforwardlast=385;
  tforward = 182:0.1:tforwardlast;
  [T Y] = ode23s(@(t,y)(model_1(t,y,k1)),tforward,initial_cond);
  v=1;
  texplast=24;
    q(texplast)=0;
    for i=1:2030 %3931/721
        for j=1:texplast
             if (tforward(i)==tdata(j))
         q(v)=Y(i,3);
         v=v+1;
             end
        end
    end
  error_in_data = sqrt(sum((q' - qdata).^2));
 end
%%
