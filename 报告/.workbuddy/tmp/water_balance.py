import math
F=96485.0; R=8.314
A=25e-4            # m^2  (25 cm^2)
j=0.5*1e4          # A/m^2 (0.5 A/cm^2)
I=j*A
n_prod=I/(2*F)     # mol/s 产水 (2e- -> 1 H2O)
n_O2=I/(4*F)
xi=2.0             # 阴极化学计量比
n_air=n_O2*xi/0.21 # mol/s 干空气供给
p=101325.0

def psat_ice(Tc):  # Magnus over ice, Pa
    return 611.15*math.exp(22.46*Tc/(272.62+Tc))

print(f"I = {I:.1f} A,  产水速率 = {n_prod:.3e} mol/s = {n_prod*18*1e3:.3f} mg/s = {n_prod*18*1e3/25:.4f} mg/(cm2*s)")
print(f"干空气供给 (stoich 2) = {n_air:.3e} mol/s\n")
print(f"{'T/℃':>6} {'p_sat/Pa':>9} {'携水/mol/s':>12} {'产水/mol/s':>12} {'排出占比%':>10}")
for Tc in [-30,-25,-20,-10,0,20,40,60,80]:
    ps=psat_ice(Tc)
    n_out=n_air*ps/(p-ps)
    print(f"{Tc:>6} {ps:>9.1f} {n_out:>12.3e} {n_prod:>12.3e} {100*n_out/n_prod:>10.3f}")

# 关键孔容 -> 填满所需水量 (mg/cm^2)
rho_ice=917.0
def fill(li,eps):  # um, porosity
    V=li*1e-6*A*eps            # m^3 孔体积
    return V*rho_ice*1e6/25    # mg/cm^2
print("\n填满各层孔隙所需冰量 (mg/cm^2):")
for name,li,eps in [("cCL (11.3um, eps=0.4)",11.3,0.4),("cGDL (150um, eps=0.6)",150,0.6),("aGDL (150um, eps=0.6)",150,0.6)]:
    m=fill(li,eps); print(f"  {name:24s} {m:8.2f}   -> 0.5A/cm2 下 {m/(n_prod*18*1e3/25):6.1f} s 填满")
