import math
F=96485.0; A=25e-4; j=0.5*1e4; I=j*A
n_prod=I/(2*F); n_air=(I/(4*F))*2/0.21; p=101325.0
def psat(Tc):
    return 611.15*math.exp(22.46*Tc/(272.62+Tc)) if Tc<0 else 611.15*math.exp(17.62*Tc/(243.12+Tc))
print(f"产水 {n_prod*18e3:.3f} mg/s = {n_prod*18e3/25:.4f} mg/(cm2*s)；干空气 {n_air:.3e} mol/s")
print(f"{'T/℃':>6}{'p_sat/Pa':>10}{'排气携水/产水 %':>18}")
for Tc in [-30,-20,-10,-5,0,20,40,60,80]:
    ps=psat(Tc); n_out=n_air*ps/(p-ps)
    print(f"{Tc:>6}{ps:>10.1f}{100*n_out/n_prod:>18.2f}")
