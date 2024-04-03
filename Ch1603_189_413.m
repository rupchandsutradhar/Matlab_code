    %Ch1603
    clear all;
    close all;
    clc;
    tic;
    format long
    global tdata_cases dna_datas initial_cond tforward
   tinitial=189;
    tforwardlast=413;
    tend=(tforwardlast-tinitial)*10+1;
  tforward = tinitial:0.1:tforwardlast;

%    % experimental data for time (days)
tdata_cases = [189; 196; 203; 238;	245; 252; 266; 280;	287; 294;...
    322; 336; 350; 364; 378; 399; 413];

 dna_datas = [83235710; 34496280; 20319840; 248341.90; 27513.90;...
     7345.10; 1377.29; 3945.27; 166.150; 281.17; 164.91; 362.87; 163.96;...
 	149.71;	149.28;	149.28;	149.28]; 
     
% Initial value of parameter
params=[25167241	  1.365e-12     190     0.450  	     1.55     0.002	     0.20  	5.82];
lb =   [25140000      1.30*10^(-12)   195      0.40        1.5      0.0018     0.19   5.8];    % lower
 ub =  [25170000      1.40*10^(-12)   200      0.45        1.6      0.0021     0.20   5.9];  % upper

  tinitialp=28;
     tforwardlastp=189;
     tendp=(tforwardlastp-tinitialp)*10+1;
   solution_previous_point=readmatrix('Ch1603_28_189_sol.txt');
  initial_cond = [solution_previous_point(tendp,2),solution_previous_point(tendp,3),...
     solution_previous_point(tendp,4),solution_previous_point(tendp,5)];
  [params,fval] = fminsearchbnd(@err_in_data, params, lb, ub, optimset('Display','iter'));

 [~, Y] = ode45(@(t,y) model_1(t,y, params), tforward, initial_cond);

figure(1)
semilogy(tforward,Y(:,3),'-');
hold on 
semilogy(tdata_cases, dna_datas, 'r.', 'MarkerSize',20);
title('Chimpanzee-1603');
% saveas(gcf,'Ch1616_182_413','fig');
% saveas(gcf,'Ch1616_182_413','jpg');
 filename = 'Ch1603_182_413.xlsx';
 D=[tforward'  Y];
 writematrix(D,'Ch1603_182_413_sol');

 
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
 filename = 'Ch1603_182_413_para.xlsx';
 writematrix(estimate_para,'Ch1603_182_413_para');
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
     tforwardlastp=189;
     tendp=(tforwardlastp-tinitialp)*10+1;
   solution_previous_point=readmatrix('Ch1603_28_189_sol.txt');
initial_cond = [solution_previous_point(tendp,2),solution_previous_point(tendp,3),...
     solution_previous_point(tendp,4),solution_previous_point(tendp,5)];
tdata = [189; 196; 203; 238;	245; 252; 266; 280;	287; 294;...
    322; 336; 350; 364; 378; 399; 413];

qdata = [83235710; 34496280; 20319840; 248341.90; 27513.90;...
     7345.10; 1377.29; 3945.27; 166.150; 281.17; 164.91; 362.87; 163.96;...
 	149.71;	149.28;	149.28;	149.28]; 
  tforwardlast=413;
  tforward = 182:0.1:tforwardlast;
  [T Y] = ode23s(@(t,y)(model_1(t,y,k1)),tforward,initial_cond);
  v=1;
  texplast=17;
    q(texplast)=0;
    for i=1:2240 %3931/721
        for j=1:texplast
             if (tforward(i)==tdata(j))
         q(v)=Y(i,3);
         v=v+1;
             end
        end
    end
  error_in_data = sqrt(sum((q' - qdata).^2))
 end
%%
