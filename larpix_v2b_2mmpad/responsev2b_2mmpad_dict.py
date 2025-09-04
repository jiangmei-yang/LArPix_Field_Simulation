import ROOT as r
import Garfield as g
import os
import ctypes
import math
import numpy as np
import time

start = time.time()
def electron_mobility(efield, temperature):
    """
    Calculation of the electron mobility w.r.t temperature and electric
    field.
    References:
     - https://lar.bnl.gov/properties/trans.html (summary)
     - https://doi.org/10.1016/j.nima.2016.01.073 (parameterization)
     
    Args:
        efield (float): electric field in kV/cm
        temperature (float): temperature
        
    Returns:
        float: electron mobility in cm^2/kV/us

    """
    a0, a1, a2, a3, a4, a5 = ELECTRON_MOBILITY_PARAMS

    num = a0 + a1 * efield + a2 * pow(efield, 1.5) + a3 * pow(efield, 2.5)
    denom = 1 + (a1 / a0) * efield + a4 * pow(efield, 2) + a5 * pow(efield, 3)
    temp_corr = pow(temperature / 89, -1.5)

    mu = num / denom * temp_corr * V / kV

    return mu

ext_x = 2;
ext_y = 2;
ext_z = 31;

lar = r.Garfield.Medium()
lar.SetTemperature(87.17) # Set the Temperature [K]
lar.EnableDrift() # Allow for drifting in this medium

# Electric charge [Q]

e = 1. # electron charge
e_SI = -1.60217733e-19 # electron charge in coulomb
coulomb = e/e_SI # coulomb = 6.24150 e+18 * e

# Energy [E]

megaelectronvolt = 1.

# Electric potential [E][Q^-1]

megavolt = megaelectronvolt / e
kilovolt = 1.e-3 * megavolt
volt = 1.e-6 * megavolt
millivolt = 1.e-3 * volt

V = volt
mV = millivolt
kV = kilovolt

ELECTRON_MOBILITY_PARAMS = 551.6, 7158.3, 4440.43, 4.29, 43.63, 0.2053

ef = range(100,30000,10)
bf = [0]
an = [0]
lar.SetFieldGrid(ef,bf,an)

for i in range(len(ef)):
    lar.SetElectronVelocityE(i,0,0,electron_mobility(ef[i]/1000,87.17)*(ef[i]/1000)/1000)
    
elm=r.Garfield.ComponentElmer("mesh.header","mesh.elements","mesh.nodes","dielectrics.dat","9x9pixel_43um_31cm_response_v2b_3.87975pitch_2pad_small_FR4_LAr.result","cm")

elm.SetMedium(81,lar)
elm.SetWeightingField("9x9pixel_43um_31cm_response_v2b_3.87975pitch_2pad_small_FR4_LAr_weight.result","readout")

# Set up a sensor object
sensor=r.Garfield.Sensor()
sensor.AddComponent(elm)
sensor.SetArea(-ext_x,-ext_y,-ext_z,ext_x,ext_y,ext_z)
sensor.AddElectrode(elm,"readout")

tmin=0
tmax = 195000
tstep= 50
nTimeBins = int((tmax-tmin)/tstep)
sensor.SetTimeWindow(tmin,tstep,nTimeBins)
drift_length = 304.31/10 #cm

drifte = r.Garfield.AvalancheMC();
drifte.SetSensor(sensor)
drifte.DisableDiffusion()
drifte.DisableAttachment()
drifte.SetTimeSteps(tstep)

cD=r.TCanvas('cD','',600,600)
vFE =r.Garfield.ViewFEMesh();
vFE.SetCanvas(cD);
vFE.SetComponent(elm);
vFE.SetArea(-ext_x,-ext_y,-ext_z,ext_x,ext_y,ext_z)
vFE.SetPlane(0,-1,0,0,0,0);
vFE.SetFillMesh(True);
vFE.SetColor(0,r.kBlue);
vFE.SetColor(82,r.kGreen+3);

viewDrift = r.Garfield.ViewDrift()
viewDrift.SetArea(-ext_x,-ext_y,-0.05,ext_x,ext_y,1)
viewDrift.SetCanvas(cD);
viewDrift.SetPlane(0,-1,0,0,0,0);
drifte.EnablePlotting(viewDrift)

response_bin_size = 3.87975/100
step = np.arange(response_bin_size/2,response_bin_size/2+response_bin_size*45,response_bin_size)
resv2b=np.zeros((45,45,3900))

for i in range(len(step)):
    for j in range(len(step)):
        print(i,j)
        sensor.ClearSignal()
        drifte.DriftElectron(step[i],step[j],drift_length,0)
        induce_s = []
        for k in range(nTimeBins):
            induce_s.append(sensor.GetSignal("readout",k))
        induce_s = np.array(induce_s)
        induce_c = induce_s*tstep*1e-15/(-1.60217733e-19)/0.05
        resv2b[i][j] = induce_c
    t = time.time()
    print("Time since start:",t-start)
response_dict = {
    "response": resv2b,
    "drift_length": drift_length,
    "time_tick": tstep
}
np.savez_compressed(f'response_38_v2b_full.npz', response = resv2b, 
                       drift_length = np.array([drift_length]), time_tick = np.array([tstep/1000]), 
                       bin_size =  np.array([response_bin_size]))
print("Done")
end = time.time()
print("Total Time:",end-start)
viewDrift.Plot(True)
vFE.Plot(True)
cD.Draw()
cD.SaveAs("v2b_2mmpad_drift_dict.png")