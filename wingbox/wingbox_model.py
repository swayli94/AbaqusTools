'''

'''
import os
import time
import json
import numpy as np

from AbaqusTools import Model, IS_ABAQUS
from wingbox_part import WingboxPart

if IS_ABAQUS:
    from abaqus import *
    from abaqusConstants import *
    from sketch import *
    import mesh


class WingboxModel(Model):
    '''
    Class for wingbox model.
    '''
    def __init__(self, name_job, pGeo, pMesh, pRun):
        super(WingboxModel,self).__init__(pGeo=pGeo, pMesh=pMesh, pRun=pRun)
        self.name_job = name_job

    def initialization(self):
        
        self.model = mdb.models[str(self.name_model)]
    
    def setup_parts(self):
        
        self.wingbox = WingboxPart(name_part='wingbox',
                        model=self.model, pGeo=self.pGeo, pMesh=self.pMesh)
        self.wingbox.build()
    
    def setup_assembly(self):
        
        a = self.rootAssembly
        
        p = self.model.parts['lofted_part']
        a.Instance(name='lofted_part', part=p, dependent=ON)
        
        for i in range(len(self.wingbox.sections)):
            p = self.model.parts['rib_%d'%(i)]
            a.Instance(name='section_%d'%(i), part=p, dependent=ON)
    
    def setup_steps(self):
        '''
        Define a static analysis step
        '''
        self.create_static_step(
            timePeriod= self.pRun['timePeriod'],
            maxNumInc=  self.pRun['maxNumInc'],
            initialInc= self.pRun['initialInc'],
            minInc=     self.pRun['minInc'],
            maxInc=     self.pRun['maxInc'],
            nlgeom=     self.pRun['nlgeom'],
        )
    
    def setup_loads(self):
        pass
        
    def setup_outputs(self):
        pass

    def setup_jobs(self):
        '''
        Define jobs
        '''
        mdb.Job(name=self.name_job, model=str(self.name_model), description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, 
            memory=self.pRun['memory_max_percentage'], 
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', resultsFormat=ODB, 
            numThreadsPerMpiProcess=self.pRun['numThreadsPerMpiProcess'], 
            multiprocessingMode=DEFAULT, 
            numCpus=self.pRun['numCpus'], numDomains=self.pRun['numCpus'], numGPUs=0)
        

if __name__ == '__main__':

    fname = 'default-parameters.json'
    with open(fname, 'r') as f:
        parameters = json.load(f)

    pGeo = parameters['pGeo']
    pMesh = parameters['pMesh']
    pRun = parameters['pRun']

    #* Build model

    name_job = 'Job_WB'

    model = WingboxModel(name_job, pGeo, pMesh, pRun)
    model.build()
    model.set_view()
    model.save_cae('WingBox.cae')
    