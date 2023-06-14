function y = Psi(x, psi)

% 计算函数 y = \psi(x) 的值
% 注意：这里 \psi 的参数 nv, p, alpha, theta 需要与 prox_tnn_psi.m 中的保持一致
% 注意：这里我们仅仅考虑 x>=0 的情况

nv = 0.01; % Cap
p = 0.618; % LP    2/3
alpha = 0.01; % MCP 1
theta = 0.5; % Log

if strcmp(psi, 'L1')
    y = Lp(x, 1);
elseif strcmp(psi, 'Lp')
    y = Lp(x, p);
elseif strcmp(psi, 'MCP')
    y = MCP(x, alpha);
elseif strcmp(psi, 'Log')
    y = Log(x, theta);
elseif strcmp(psi, 'CapL1')
    y = CapLp(x, 1, nv);
elseif strcmp(psi, 'CapLp')
    y = CapLp(x, p, nv);
elseif strcmp(psi, 'CapMCP')
    y = CapMCP(x, alpha, nv);
elseif strcmp(psi, 'CapLog')
    y = CapLog(x, theta, nv);
end

end
