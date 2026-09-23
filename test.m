clear;
clc;

p = params();
g = build_grid(p);
s = state(g);

%% Current-distribution test

j = 0.3e4;     % 0.3 A/cm^2 = 3000 A/m^2

[is,im,iv,currentInfo] = ...
    calc_current_distribution(j,g,p);

fprintf('\nCurrent distribution:\n');
fprintf('  Applied j       = %.3f A/m^2\n',j);
fprintf('  max |is+im-j|   = %.3e A/m^2\n', ...
    currentInfo.currentError);
fprintf('  face current err = %.3e A/m^2\n', ...
    currentInfo.faceCurrentError);

assert(currentInfo.currentConserved);


%% Check pure conduction regions

assert(all(abs(is(g.idx_aGDL)-j) < 1e-10));
assert(all(abs(im(g.idx_aGDL))   < 1e-10));

assert(all(abs(is(g.idx_PEM))    < 1e-10));
assert(all(abs(im(g.idx_PEM)-j)  < 1e-10));

assert(all(abs(is(g.idx_cGDL)-j) < 1e-10));
assert(all(abs(im(g.idx_cGDL))   < 1e-10));


%% Check catalyst-layer conversion and face values

tolCurrent = max(1e-10,1e-12*max(abs(j),1));

% iv is the positive reaction-rate magnitude in both catalyst layers.
assert(abs(sum(iv(g.idx_aCL).*g.dx(g.idx_aCL))-j) < tolCurrent);
assert(abs(sum(iv(g.idx_cCL).*g.dx(g.idx_cCL))-j) < tolCurrent);

% Proton current rises through aCL and falls through cCL.
assert(all(diff(im(g.idx_aCL)) > 0));
assert(all(diff(im(g.idx_cCL)) < 0));

% Signed proton-current gradient is separate from positive iv.
assert(all(currentInfo.dimdx(g.idx_aCL) > 0));
assert(all(currentInfo.dimdx(g.idx_cCL) < 0));

% Face-current conservation and exact interface values.
assert(max(abs(currentInfo.isFace+currentInfo.imFace-j)) < tolCurrent);
assert(abs(currentInfo.imFace(g.face_aGDL_aCL)) < tolCurrent);
assert(abs(currentInfo.imFace(g.face_aCL_PEM)-j) < tolCurrent);
assert(abs(currentInfo.imFace(g.face_PEM_cCL)-j) < tolCurrent);
assert(abs(currentInfo.imFace(g.face_cCL_cGDL)) < tolCurrent);


%% Plot

figure;

plot(g.x_um,is/1e4,'o-','LineWidth',1.5);
hold on;

plot(g.x_um,im/1e4,'s-','LineWidth',1.5);

xlabel('x / \mum');
ylabel('Current density / A cm^{-2}');

legend('i_s','i_m','Location','best');

grid on;

title('Electronic and protonic current distribution');
