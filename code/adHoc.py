#%%
import sys
sys.path.append('/home/MntStrBrk/BSPT')
from BSPT import tools as bspt

# %%
bspt.hello_world("haixi")
# %%
com1 = bspt.Computer('i5', 36)
com2 = bspt.Computer('i7', 16)
# %%
com1.config()
com2.config()
# %%
