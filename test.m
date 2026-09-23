clear;
clc;

p = params();
g = build_grid(p);
s = state(g);

fprintf('Number of grid cells : %d\n',g.N);
fprintf('Number of states     : %d\n\n',s.Nx);

fprintf('T   : %3d states, x(%d:%d)\n', ...
    s.n.T, s.block.T);

fprintf('H2  : %3d states, x(%d:%d)\n', ...
    s.n.H2, s.block.H2);

fprintf('O2  : %3d states, x(%d:%d)\n', ...
    s.n.O2, s.block.O2);

fprintf('mw  : %3d states, x(%d:%d)\n', ...
    s.n.mw, s.block.mw);

fprintf('mi  : %3d states, x(%d:%d)\n', ...
    s.n.mi, s.block.mi);
