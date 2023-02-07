function [ind_winner,time_winner] = function_read_out(t, iFR, r_w, r_tau, r_th, t_interval)
  dt = t(2) - t(1);
  iFRsize = size(iFR);
  nPop = iFRsize(1:end-1); % we discard the temporal dimension

  if nargin < 3
    r_w = 0.05*nPop;
  end

  if nargin < 4
    r_tau = 100; % ms
    % r_tau = 5; % ms
  end

  if nargin < 5
    r_th = 50; % sp/s
  end

  if nargin < 6
    t_interval = [t(1) t(end)];
  end

  r0 = zeros(nPop(:,1), 1)
  for iPop = 1:numel(r0)
    [time readout(iPop,:)] = ode45(@(t,r_prev) r_odefun(r_prev,r_tau,r_w(iPop,:),iFR(iPop,:,1+floor(t/dt))),t,r0(iPop));
  end

  ind1 = find(time >= t_interval(1) - 0.5*dt,1,'first')
  ind2 = find(time >= t_interval(2) - 0.5*dt,1,'first')

  fontSize = 16;

  figure
  plot(time,readout)
  hold on
  plot(minmax(t),[r_th r_th],'k--')
%   xlim(t_interval)
  title(sprintf('r_{tau} = %.3g, r_{w1} = %.3g, r_{w2} = %.3g, r_{th} = %.3g',r_tau,unique(r_w(:,1)),unique(r_w(:,2)),r_th))
  legend('r_A','r_B','Thr','Location','Best')
  xlabel('Time (ms)','fontSize',fontSize)
  ylabel('Readout','fontSize',fontSize)

  ind_winner = inf;
  time_winner = inf;
  for iPop = randperm(unique(nPop)) % when distinct populations reach the threshold simultaneously (same dt), the winner among them is chosen randomly
    ind_time_xthreshold = ind1 + min([inf, find(readout(iPop, ind1:ind2) >= r_th, 1,'first')]);
    if ~isinf(ind_time_xthreshold)
      time_xthreshold = time(ind_time_xthreshold);
      if time_xthreshold < time_winner
        ind_winner = iPop;
        time_winner = time_xthreshold;
      end
    end
  end
end

function rdot = r_odefun(r,r_tau,r_w,iFR)
  r_input = r_w*iFR';
  rdot= -r/r_tau + r_input.*(r_input > 0);
end
