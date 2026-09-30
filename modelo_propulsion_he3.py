import numpy as np, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
g0=9.80665
BLUE,ORANGE,AQUA,INK,MUTED,GRID="#2a78d6","#eb6834","#1baf7a","#0b0b0b","#52514e","#e4e3df"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.edgecolor":MUTED,"axes.labelcolor":INK,
 "xtick.color":MUTED,"ytick.color":MUTED,"axes.spines.top":False,"axes.spines.right":False,"axes.grid":True,
 "grid.color":GRID,"grid.linewidth":0.8,"lines.linewidth":2,"figure.dpi":200})
R={}
# ---- Validation 1: Dawn (NSTAR) ----
m0,mxe=1217.7,425.0; mf=m0-mxe
for isp in (3100,2650,1900):
    R[f"dawn_dv_isp{isp}"]=isp*g0*np.log(m0/mf)/1e3
F_nstar=0.092; ve=3100*g0; mdot=F_nstar/ve
R["dawn_mdot"]=mdot; R["dawn_tburn_yr"]=mxe/mdot/3.156e7
R["dawn_accel"]=F_nstar/m0
# ---- Validation 2: Aime et al. 2021 DFD Haumea ----
m0d,mpd,Fd,ispd=7488.0,3595.0,8.0,10000.0
ved=ispd*g0
R["dfd_dv"]=ved*np.log(m0d/(m0d-mpd))/1e3
R["dfd_mdot"]=Fd/ved; R["dfd_tburn_d"]=mpd/(Fd/ved)/86400
R["dfd_accel"]=Fd/m0d
R["accel_ratio"]=R["dfd_accel"]/R["dawn_accel"]
# jet powers
R["next_jet_W"]=0.5*0.237*4190*g0; R["next_eta"]=R["next_jet_W"]/6900
R["dfd_jet_W"]=0.5*Fd*ved
R["nstar_jet_W"]=0.5*F_nstar*ve
# ---- He-3 energetics ----
MeV=1.602176634e-13; u=1.66053907e-27
E=18.35*MeV; mHe3=3.01603*u
R["J_per_kg_He3"]=E/mHe3; R["GJ_per_kg_He3"]=E/mHe3/1e9
yr=3.156e7
R["kgHe3_per_yr_2MW"]=2e6*yr/(E/mHe3)
R["kgHe3_510d_2MW"]=2e6*510*86400/(E/mHe3)
R["regolith_t_per_kgHe3_20ppb"]=1/20e-9/1000
R["dawn_tburn_ratio"]=None
# ---- Lawson (simplified, NRL rates) ----
T=np.array([1,2,5,10,20,50,100,200,500,1000.])
DT=np.array([5.5e-21,2.6e-19,1.3e-17,1.1e-16,4.2e-16,8.7e-16,8.5e-16,6.3e-16,3.7e-16,2.7e-16])*1e-6
DHe=np.array([1.0e-26,1.4e-23,6.7e-21,2.3e-19,3.8e-18,5.4e-17,1.6e-16,2.4e-16,2.3e-16,1.8e-16])*1e-6
Tf=np.logspace(np.log10(3),np.log10(1000),400)
from scipy.interpolate import PchipInterpolator
sv=lambda tab:np.exp(PchipInterpolator(np.log(T),np.log(tab))(np.log(Tf)))
keV=1.602176634e-16; CB=5.34e-37
svDT,svDHe=sv(DT),sv(DHe)
# D-T: charged 3.5 MeV, thermal 3nT; D-He3: charged 18.35 MeV, thermal 3.75nT; n=total ion density, equimolar
tp_DT=12*Tf**2/(svDT*3500)                     # keV s m^-3, no brem
tp_DHe=15*Tf**2/(svDHe*18350)
# with bremsstrahlung: nτE = W/(Pc/n^2 - Pb/n^2)
PcDT=svDT*3.5*MeV/4; PbDT=CB*1*np.sqrt(Tf)   # n_e=n, sum n_i Z^2 = n
PcDHe=svDHe*E/4; PbDHe=CB*1.5*2.5*np.sqrt(Tf)
def withb(Pc,Pb,wcoef):
    d=Pc-Pb; out=np.where(d>0,wcoef*Tf*keV/np.where(d>0,d,1)*Tf,np.nan); return out
tpb_DT=withb(PcDT,PbDT,3.0); tpb_DHe=withb(PcDHe,PbDHe,3.75)
for name,arr in [("DT",tp_DT),("DHe",tp_DHe),("DTb",tpb_DT),("DHeb",tpb_DHe)]:
    i=np.nanargmin(arr); R[f"min_{name}"]=float(arr[i]); R[f"Tmin_{name}"]=float(Tf[i])
R["ratio_nobrem"]=R["min_DHe"]/R["min_DT"]; R["ratio_brem"]=R["min_DHeb"]/R["min_DTb"]
# ideal ignition temp (Pc=Pb)
for n,(a,b) in {"DT":(PcDT,PbDT),"DHe":(PcDHe,PbDHe)}.items():
    k=np.where(a>b)[0]; R[f"Tign_{n}"]=float(Tf[k[0]]) if len(k) else None
json.dump({k:(float(v) if v is not None else None) for k,v in R.items()},open("/home/claude/paper/results.json","w"),indent=1)
for k,v in R.items(): print(k,v)

# ===== FIG 1: Tsiolkovsky mass ratio =====
dv=np.linspace(0,80,400)
fig,ax=plt.subplots(figsize=(6.4,3.8))
for isp,c,lab in [(3100,BLUE,"Iónico NSTAR (3 100 s, vuela desde 1998)"),(4190,AQUA,"Iónico NEXT (4 190 s, calificado)"),(10000,ORANGE,"Fusión DFD (≈10 000 s, conceptual)")]:
    ax.plot(dv,np.exp(dv*1e3/(isp*g0)),color=c,label=lab)
ax.scatter([R["dawn_dv_isp3100"]],[m0/mf],s=45,color=BLUE,edgecolor="white",zorder=5)
ax.annotate("Dawn (real):\n1 217,7 → 792,7 kg",(R["dawn_dv_isp3100"],m0/mf),xytext=(19,1.2),fontsize=8.5,color=INK,arrowprops=dict(arrowstyle="-",color=MUTED))
ax.scatter([R["dfd_dv"]],[m0d/(m0d-mpd)],s=45,color=ORANGE,edgecolor="white",zorder=5)
ax.annotate("Aime et al. (2021),\nmisión a Haumea",(R["dfd_dv"],m0d/(m0d-mpd)),xytext=(48,3.2),fontsize=8.5,color=INK,arrowprops=dict(arrowstyle="-",color=MUTED))
ax.set_yscale("log"); ax.set_ylim(1,12); ax.set_xlim(0,80)
ax.set_yticks([1,2,3,5,10]); ax.set_yticklabels(["1","2","3","5","10"])
ax.set_xlabel("Δv requerido (km/s)"); ax.set_ylabel("Razón de masas m₀/m_f")
ax.legend(frameon=False,fontsize=8.5,loc="upper left")
fig.tight_layout(); fig.savefig("/home/claude/paper/fig/fig1_tsiolkovsky.png"); plt.close()

# ===== FIG 2: thrust vs Isp at fixed power =====
isp=np.linspace(1000,20000,400)
fig,ax=plt.subplots(figsize=(6.4,3.8))
for Pj,c,lab in [(4.9e3,BLUE,"Potencia de chorro 4,9 kW (un NEXT)"),(49e3,AQUA,"49 kW (10 NEXT en paralelo)"),(392e3,ORANGE,"392 kW (DFD 2 MW, Aime et al.)")]:
    ax.plot(isp,2*Pj/(isp*g0),color=c,label=lab)
ax.scatter([4190],[0.237],s=45,color=BLUE,edgecolor="white",zorder=5); ax.annotate("NEXT PM1\n0,237 N",(4190,0.237),xytext=(6200,0.55),fontsize=8.5,arrowprops=dict(arrowstyle="-",color=MUTED))
ax.scatter([10000],[8],s=45,color=ORANGE,edgecolor="white",zorder=5); ax.annotate("DFD 2 MW\n8 N",(10000,8),xytext=(5200,40),fontsize=8.5,arrowprops=dict(arrowstyle="-",color=MUTED))
ax.set_yscale("log"); ax.set_xlabel("Impulso específico Isp (s)"); ax.set_ylabel("Empuje F (N)")
ax.legend(frameon=False,fontsize=8.5,loc="upper right")
fig.tight_layout(); fig.savefig("/home/claude/paper/fig/fig2_empuje_isp.png"); plt.close()

# ===== FIG 3: burn time vs dv =====
dvs=np.linspace(0.1,80,400)
def tb(m0,mdot,ve,dvv): return m0/mdot*(1-np.exp(-dvv*1e3/ve))
fig,ax=plt.subplots(figsize=(6.4,3.8))
tA=tb(m0,mdot,ve,dvs)/86400; limA=ve*np.log(m0/mf)/1e3
mask=dvs<=limA
ax.plot(dvs[mask],tA[mask],color=BLUE,label="Nave clase Dawn (1 NSTAR, 92 mN)")
ax.plot(dvs,tb(m0d,Fd/ved,ved,dvs)/86400,color=ORANGE,label="Nave DFD 2 MW (8 N)")
ax.axvline(limA,color=BLUE,lw=1,ls=(0,(4,3))); ax.text(limA+1,3,"límite de xenón\nde Dawn",fontsize=8.5,color=MUTED)
ax.scatter([R["dfd_dv"]],[R["dfd_tburn_d"]],s=45,color=ORANGE,edgecolor="white",zorder=5)
ax.annotate("510 días (publicado)",(R["dfd_dv"],R["dfd_tburn_d"]),xytext=(42,1500),fontsize=8.5,arrowprops=dict(arrowstyle="-",color=MUTED))
ax.set_yscale("log"); ax.set_xlabel("Δv (km/s)"); ax.set_ylabel("Tiempo de encendido (días)")
ax.set_xlim(0,80); ax.legend(frameon=False,fontsize=8.5,loc="lower right")
fig.tight_layout(); fig.savefig("/home/claude/paper/fig/fig3_tiempo.png"); plt.close()

# ===== FIG 4: Lawson =====
fig,ax=plt.subplots(figsize=(6.4,3.8))
ax.plot(Tf,tpb_DT,color=BLUE,label="D–T (con bremsstrahlung)")
ax.plot(Tf,tpb_DHe,color=ORANGE,label="D–³He (con bremsstrahlung)")
ax.plot(Tf,tp_DT,color=BLUE,lw=1.2,ls=(0,(4,3)),label="D–T (sin radiación)")
ax.plot(Tf,tp_DHe,color=ORANGE,lw=1.2,ls=(0,(4,3)),label="D–³He (sin radiación)")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(3,1000); ax.set_ylim(1e21,1e24)
ax.set_xlabel("Temperatura iónica T (keV)"); ax.set_ylabel("n·T·τ_E para ignición (keV·s·m⁻³)")
ax.legend(frameon=False,fontsize=8.5,loc="lower right")
fig.tight_layout(); fig.savefig("/home/claude/paper/fig/fig4_lawson.png"); plt.close()
