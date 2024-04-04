    %Ch1603
    clear all;
    close all;
    clc;
    tic;
    format long
    global tdata_cases dna_datas initial_cond tforward
   tinitial=28;
    tforwardlast=189;
    tend=(tforwardlast-tinitial)*10+1;
  tforward = tinitial:0.1:tforwardlast;


tdata_cases = [28; 42; 49; 56; 63; 70; 77; 84; 91; 98; 105;  112; 119;...
    126; 133; 140; 147;	154; 161; 168;	175; 182; 189];


 dna_datas = [100; 851.56; 5390.82;	24005.80; 197851.50; 1252372; 6650129;...
 	24839890; 204725600; 641227000; 910224500; 585574900; 584773700;...
    489680900; 343950200; 221250300; 408930200; 263049400; 262634197;...
    184492800;	154491600; 168462900; 83235710]; 
%   params = lambda,        k,           a,     gamma,     beta,     mu,      delta,     c]    
 params=[25104026	2.259e-12	      198	 0.214	   0.864	  0.045     0.105 	   3.93];
 lb = [25140000      2.25*10^(-12)    195      0.2       0.8      0.040     0.100       3.9];   
 ub = [2516000      2.28*10^(-12)     200      0.25      0.9      0.050     0.110       3.93];   
 % initial condition
 initial_cond = [2*10^9.9  0 100 0];

  [params,fval] = fminsearchbnd(@err_in_data, params, lb, ub, optimset('Display','iter'));


 [~, Y] = ode45(@(t,y) model_1(t,y, params), tforward, initial_cond);

figure(1)
semilogy(tforward,Y(:,3),'-');
hold on 
semilogy(tdata_cases, dna_datas, 'r.', 'MarkerSize',20);
title('Chimpanzee-1603');
% saveas(gcf,'Ch1616_28_189','fig');
% saveas(gcf,'Ch1616_28_189','jpg');
 filename = 'Ch1603_28_189.xlsx';
 D=[tforward'  Y];
 writematrix(D,'Ch1603_28_189_sol');

 
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
 filename = 'Ch1603_28_189_para.xlsx';
 writematrix(estimate_para,'Ch1603_28_189_para');
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

tdata = [28; 42; 49; 56; 63; 70; 77; 84; 91; 98; 105;  112; 119;...
    126; 133; 140; 147;	154; 161; 168;	175; 182; 189];

qdata = [100; 851.56; 5390.82;	24005.80; 197851.50; 1252372; 6650129;...
 	24839890; 204725600; 641227000; 910224500; 585574900; 584773700;...
    489680900; 343950200; 221250300; 408930200; 263049400; 262634197;...
    184492800;	154491600; 168462900; 83235710]; 
  tforwardlast=189;
  tforward = 28:0.1:tforwardlast;
  [T Y] = ode23s(@(t,y)(model_1(t,y,k1)),tforward,[2*10^9.9  0 100 0]);
  v=1;
  texplast=23;
    q(texplast)=0;
    for i=1:1610 %3931/721
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
