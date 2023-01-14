function xi = recover(r,s)
% Func: get xi from r and s
  xi=[r(:,1).*s(:,1),(1-r(:,1)).*s(:,1), r(:,2).*s(:,2),(1-r(:,2)).*s(:,2)];
end

