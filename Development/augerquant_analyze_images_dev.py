# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 13:23:25 2017

@author: tkc
"""
import os
from skimage.feature import register_translation
if 'C:\\Users\\tkc\\Documents\\Python_Scripts\\Montage' not in sys.path:
    sys.path.append('C:\\Users\\tkc\\Documents\\Python_Scripts\\Montage')
import Image_montage_functions as mont

#%%
os.chdir('H:\\Research_data\\Thin_sections\\Auger\\Acfer094\\6Sep17')
chargemap1=np.load('Acfer094_Omap_area8_chargemap.npy')

os.chdir('H:\\Research_data\\Thin_sections\\Auger\\Acfer094\\27Sep17')
chargemap2=np.load('Acfer094_Omap_area8_chargemap.npy')

# register translation (pretty accurate)
shift, error, diffphase = register_translation(chargemap1, chargemap2)

plt.imshow(chargemap1)
plt.imshow(chargemap2)

# correlation coeff between chargemaps
np.corrcoef(chargemap1, chargemap2)

# other method of determining relative shifts between images (ORB detector)
ORBshift, scaling, rotation = mont.findORBshift(chargemap1, chargemap2)
