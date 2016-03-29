import models as m
import pickle

name='WL16-Ulrich-Env_Mass=0.1-Rc=2000-Lsun=5.6'
folder='WL16/'
WL16 = m.modelLoad(folder,name)

WL16.plotModel(show=True)

