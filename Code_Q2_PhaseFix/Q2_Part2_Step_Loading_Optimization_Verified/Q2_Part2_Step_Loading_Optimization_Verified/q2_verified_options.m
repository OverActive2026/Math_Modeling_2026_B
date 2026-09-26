function opt=q2_verified_options()
opt=q2_optimization_options();
opt.RelTol=5e-5; opt.AbsTol=5e-9;
opt.coldTMaxS=600; opt.progressIntervalS=10;
assert(opt.gammaIce==3.5 && opt.kFreeze==.4);
assert(isequal(opt.nCellLayer,[2 3 4 4 2]));
assert(opt.minimumVoltageV==.30 && opt.qMaxCcm2==20 && opt.jMaxAcm2==.5);
end
