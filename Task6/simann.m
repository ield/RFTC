function [xopt,fopt,exitflag,output]=...
                       simann(func,x,LB,UB,sa_t,sa_rt,sa_nt,sa_ns,dispiter)
% SIMANN Simulated Annealing programmed for minimization problem
%      [xopt,fopt]=SIMANN(func,x,LB,UB,sa_t,sa_rt,sa_nt,sa_ns,dispiter)
%
% INPUTS
% func     string variable containing name of function file to be optimized 
% x        starting values (if x=NaN uses random starting values)
% LB       lower bound on optimization parameters
% UB       upper bound on optimization parameters
% sa_t     initial temperature
% sa_rt    temperature reduction factor, 0 < sa_rt < 1, try .85
% sa_nt    number of times through ns loop before temperature reduction 
%                                               (recommended value: 5)
% sa_ns    number of times through function before stepsize adjustment 
%                                              (recommended value: 20)
% dispiter true to show iteration results
%
% OUTPUTS
% xopt     the optimal solution
% fopt     function value at the optimal solution
% exitflag reasons for termination: -1 MaxFunEvals exceeded
%                                    1 Minimum step size reached
%                                    2 Convergence criteria reached
% output   information about process
%
% Translated from original Fortran code in 1999 by Knoek van Soest 
%                                                 <a_j_van_soest@fbw.vu.nl>
% Bug fix in 2000 by Anne Su <sua@bme.ri.ccf.org>
% Minor changes in 2016,2018 by Jaime Esteban
%
% things to do: documentation, input and output argument passing, 
%               printing, error checking
% 
% See original Fortran code by Bill Goffe for documentation and literature 
% references. This code is at 
% http://netlib.org/cgi-bin/netlibfiles.pl?filename=/opt/simann.f
% and at
% http://ideas.repec.org/c/cod/fortra/simanneal.html

x=x(:)';
LB=LB(:)';                                   
UB=UB(:)';
sa_nargs=length(LB);              %number of parameters

sa_neps=4;    %number of times eps tolerance is achieved before termination. Number of failures to stop the program.
sa_TolFun=eps;                   %convergence criteria. The comparison done f'-fopt
sa_MaxFunEvals=40000*sa_nargs;   %maximum number of function evaluations. Other stopping options

sa_nacc=0;                        %number of acceptions
sa_nevals=0;                      %number of evaluations

fstar=Inf*ones(sa_neps,1);        %optima in last sa_neps iterations

VM=(UB-LB);%/2;                      %maximum step size
sa_TolX=sqrt(eps);                       %minimum step size

if any(isnan(x))
  x=LB+(UB-LB).*rand(1,sa_nargs);   %starting values for model parameters
end
f=feval(func,x);                    %function evaluation with parameters x
if dispiter
  disp(['initial loss function value: ' num2str(f)])
end
sa_nevals=sa_nevals+1;
xopt=x;
fopt=f;
fstar(1)=f;

%LOOP
while true
  
  nup=0;                                   %number of uphill movements. Number of acceptances that worsen the fucntion
  nrej=0;                                  %number of rejections.
  ndown=0;                                 %number of downhill movements. Number of acceptances that improve the function.
  nacp=zeros(sa_nargs,1);
  
  for m=1:sa_nt
    for j=1:sa_ns
      for h=1:sa_nargs
        if sa_nevals>=sa_MaxFunEvals
          disp('Too many function evaluations')
          exitflag=-1;
          output.funccount=sa_nevals;
          output.currenttemperature=sa_t;
          output.currentstepsize=VM;
          return
        end
        % generate xp, trial value of x
        xp=x;
        xp(h)=x(h)+VM(h)*(2*rand(1,1)-1.0); %calculate new value for x (xp)
        if (xp(h)<LB(h)) || (xp(h)>UB(h))
          xp(h)=LB(h)+(UB(h)-LB(h))*rand(1,1);
        end       
        % evaluate at xp and return as fp
        fp=feval(func,xp);          %function evaluation with parameters xp
        sa_nevals=sa_nevals+1;

        % we minimize! accept if the function value decreases
        if fp<=f
          x=xp;
          f=fp;
          sa_nacc=sa_nacc+1;
          nacp(h)=nacp(h)+1;
          ndown=ndown+1;
          % if smaller than any previous point, record as new optimum
          if fp<fopt
            xopt=xp;
            fopt=fp;
          end
        else % function value increases
          p=exp((f-fp)/sa_t);                              %random number
          pp=rand(1,1);
          if pp<p
            x=xp;
            f=fp;
            sa_nacc=sa_nacc+1;
            nacp(h)=nacp(h)+1;
            nup=nup+1;
          else
            nrej=nrej+1;
          end
        end
      end  % end of parameters loop
    end  % end of repetitions for a given step size
    
    % adjust maximal step size vm
    for i=1:sa_nargs
      ratio=nacp(i)/sa_ns;
      if ratio>0.6
        VM(i)=VM(i) * (1+5*(ratio-0.6));
      elseif ratio <0.4
        VM(i)=VM(i) / (1+5*(0.4-ratio));
      end
      if VM(i)>(UB(i)-LB(i))
        VM(i)=UB(i)-LB(i);
      end
    end
    if all(VM./xopt<sa_TolX)
       disp('Minimum step size reached')
       exitflag=1;
       output.funccount=sa_nevals;
       output.currenttemperature=sa_t;
       output.currentstepsize=VM;
       return
    end
    nacp=zeros(sa_nargs,1);

    % provide statistics about current state of optimization
    
    if dispiter
      disp(' ');
      disp(['No. of evaluations: ' num2str(sa_nevals)]);
      disp(['  current temperature: ' num2str(sa_t)]);
      disp(['  current optimum function value: ' num2str(fopt)]);
      disp(['  No. of downhill steps: ' num2str(ndown)]);
      disp(['  No. of accepted uphill steps: ' num2str(nup)]);
                           % we minimize, thus downhill is always accepted!
      disp(['  No. of rejections: ' num2str(nrej)]);
      disp(['  current parameter values: ' num2str(xp)]);
      disp(['  current optimum vector: ' num2str(xopt)]);
      disp(['  current step size: ' num2str(VM)]);
    end

  end  % end of repetitions for a given temperature
  
  
  % check termination criteria
  fstar(1)=f;
  quit = ((fstar(1)-fopt) <= sa_TolFun);
  if any(abs(fstar-f)>sa_TolFun)
    quit=0;
  end
  
  if quit
    if dispiter
      disp(['Sim. ann. termination after ', num2str(sa_nevals),' evals']);
    end
    exitflag=2;
    output.funccount=sa_nevals;
    output.currenttemperature=sa_t;
    output.currentstepsize=VM;
    return
  end
  
  % reduce temperature  
  sa_t=sa_t*sa_rt;
  fstar(2:end)=fstar(1:end-1);
  % continue from current optimum
  x=xopt;
  f=fopt;
  
end