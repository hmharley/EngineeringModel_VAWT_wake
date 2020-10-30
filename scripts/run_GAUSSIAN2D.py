import sys, os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'src'))

import math
import numpy as np
import EngineeringModels as EM
import matplotlib.pyplot as plt

x = np.linspace(1,6,50)
y = 0; z = 0; zh = 0
k_star = 0.04; Ct = 0.94
D = 1; H=3; k_starD=k_star; k_starH =k_star;
DU_U0,sigmaD,sigmaH = EM.GaussianVAWT(Ct, x,y,z,zh,D,H, k_starD, k_starH);
Ur_over_U0 = 1 - DU_U0


plt.plot(x,Ur_over_U0,'--r')
plt.show()

# %%
x = np.linspace(0,6,50)
Dr = 1;
Ct = 0.9;
H = 3;

k = 0.066;
kD = k;kH = k;
Ur_over_U0,Hw, Dw = EM.TophatVAWT(Ct,H,Dr,x,kD,kH)
plt.plot(x,Ur_over_U0,'--r')
plt.show()
# %%
