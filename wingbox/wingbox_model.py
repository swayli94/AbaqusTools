'''

'''
import os
import time
import json
import numpy as np

from AbaqusTools import Model, IS_ABAQUS
from lofting_part import LoftingPart
from rib_part import RibPart

if IS_ABAQUS:
    from abaqus import mdb
    from abaqusConstants import *


class WingboxModel(Model):
    '''
    Class for wingbox model.
    '''
    def __init__(self, name_job, pGeo, pMesh, pRun):
        super(WingboxModel,self).__init__(pGeo=pGeo, pMesh=pMesh, pRun=pRun)
        self.name_job = name_job

    def initialization(self):
        
        self.model = mdb.models[str(self.name_model)]
        
        self.create_material_IM785517(elastic_type='ENGINEERING_CONSTANTS')

    def setup_parts(self):
        
        self.lofting = LoftingPart(name_part='lofting',
                        model=self.model, pGeo=self.pGeo, pMesh=self.pMesh)
        self.lofting.build()
        
        self.ribs = []
        for i in range(len(self.pGeo['sections'])):
            rib = RibPart(model=self.model, pGeo=self.pGeo, pMesh=self.pMesh, index_rib=i)
            rib.build()
            self.ribs.append(rib)

    def setup_assembly(self):
        
        a = self.rootAssembly
        p = self.model.parts['lofting']
        a.Instance(name='lofting', part=p, dependent=ON)
        
        for i, rib in enumerate(self.ribs):
            p = self.model.parts[rib.name_part]
            a.Instance(name=rib.name_part, part=p, dependent=ON)
    
    def setup_steps(self):

        if self.pRun['nlgeom']:
            nlgeom = ON
        else:
            nlgeom = OFF
        
        self.model.StaticStep(
            name= 'Loading', 
            previous= 'Initial', 
            description= 'Static simulation',
            timePeriod= self.pRun['timePeriod'],
            maxNumInc= self.pRun['maxNumInc'],
            initialInc= self.pRun['initialInc'],
            minInc=self.pRun['minInc'],
            maxInc=self.pRun['maxInc'],
            nlgeom=nlgeom,
            )
            
    def setup_interactions(self):
        
        n_sections = len(self.lofting.sections)

        #* Tie between covers and ribs
        pairs = []
        for i_section in range(n_sections-1):
            if i_section < n_sections-2:
                index_ribs = [i_section]
            else:
                index_ribs = [i_section, i_section+1]
            for side in ['upper', 'lower']:
                nsi = 'lofting'
                nss = 'face_wingbox%d_cover_%s' % (i_section, side)
                for k in index_ribs:
                    nmi = self.ribs[k].name_part
                    nms = 'edge_rib_cover_%s' % (side)
                    pairs.append((nmi, nms, nsi, nss))

        for i_pair, (nmi, nms, nsi, nss) in enumerate(pairs):
            name_constraint = 'tie_cover_rib_%d' % i_pair
            self._create_tie_constraint(name_constraint,
                    name_master_instance=nmi, name_master_set=nms,
                    name_slave_instance=nsi, name_slave_set=nss)
        
        #* Tie between spars and ribs
        pairs = []
        for i_section in range(n_sections-1):
            section = self.lofting.sections[i_section]
            if i_section < n_sections-2:
                index_ribs = [i_section]
            else:
                index_ribs = [i_section, i_section+1]
            for j in range(section.n_spars):
                nmi = 'lofting'
                nms = 'face_wingbox%d_spar_%d' % (i_section, j)
                for k in index_ribs:
                    nsi = self.ribs[k].name_part
                    nss = 'edge_rib_spar_%d' % (j)
                    pairs.append((nmi, nms, nsi, nss))
        
        for i_pair, (nmi, nms, nsi, nss) in enumerate(pairs):
            name_constraint = 'tie_spar_rib_%d' % i_pair
            self._create_tie_constraint(name_constraint,
                    name_master_instance=nmi, name_master_set=nms,
                    name_slave_instance=nsi, name_slave_set=nss)
        
        #* Tie between stringers and ribs
        pairs = []
        for i_section, section in enumerate(self.lofting.sections):
            nmi = self.ribs[i_section].name_part
            nms = 'face_rib'
            for j in range(section.n_stringers):
                for side in ['upper', 'lower']:
                    for feature in ['web', 'flange']:
                        nsi = 'lofting'
                        nss = 'edge_sec%d_stringer_%d_%s_%s' % (i_section, j, side, feature)
                        pairs.append((nmi, nms, nsi, nss))
                    
        for i_pair, (nmi, nms, nsi, nss) in enumerate(pairs):
            name_constraint = 'tie_stringer_rib_%d' % i_pair
            self._create_tie_constraint(name_constraint,
                    name_master_instance=nmi, name_master_set=nms,
                    name_slave_instance=nsi, name_slave_set=nss)
        
    def _create_tie_constraint(self, name_constraint,
                    name_master_instance, name_master_set,
                    name_slave_instance, name_slave_set):
        '''
        Create tie constraint.
        '''
        a = self.rootAssembly
        region1=a.instances[name_master_instance].sets[name_master_set]
        region2=a.instances[name_slave_instance].sets[name_slave_set]
        mdb.models['Model-1'].Tie(name=name_constraint, main=region1, secondary=region2, 
            positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
    
    def setup_loads(self):
        
        a = self.rootAssembly
        n_sections = len(self.lofting.sections)
        
        region = a.instances['rib_0'].sets['all']
        self.model.EncastreBC(name='sec0', createStepName='Initial', 
            region=region, localCsys=None)
        
        for side in ['upper', 'lower']:
            for i_section in range(n_sections-1):
                region = a.instances['lofting'].surfaces['face_wingbox%d_cover_%s' % (i_section, side)]
                self.model.Pressure(name='pressure_%d_%s' % (i_section, side),
                    createStepName='Loading', region=region, 
                    distributionType=UNIFORM, field='',
                    magnitude=self.pMesh['pressure_%s_cover' % side], 
                    amplitude=UNSET)
        
    def setup_outputs(self):
        
        self.model.FieldOutputRequest(name='F-Output-1', 
            createStepName='Loading', variables=('S', 'E', 'U', 'RF', 'P'), 
            frequency=LAST_INCREMENT)

        for i, name_layup in enumerate(self.lofting.name_layups):
            self.model.FieldOutputRequest(name='Layup-Output-lofting-%d' % i, 
                createStepName='Loading', variables=('S', 'E', 'SDV'),
                frequency=LAST_INCREMENT, 
                layupNames=('lofting.%s' % name_layup, ), 
                layupLocationMethod=SPECIFIED, outputAtPlyTop=False, outputAtPlyMid=True, 
                outputAtPlyBottom=False, rebar=EXCLUDE)
        
        for i, rib in enumerate(self.ribs):
            for j, name_layup in enumerate(rib.name_layups):
                self.model.FieldOutputRequest(name='Layup-Output-rib%d-%d' % (i, j), 
                    createStepName='Loading', variables=('S', 'E', 'SDV'),
                    frequency=LAST_INCREMENT, 
                    layupNames=('%s.%s' % (rib.name_part, name_layup), ),
                    layupLocationMethod=SPECIFIED, outputAtPlyTop=False, outputAtPlyMid=True, 
                    outputAtPlyBottom=False, rebar=EXCLUDE)

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
    
    if pMesh['failure_model'] == 'LaRC05':
        model.write_job_inp()
        model.write_IM785517_property_table_inp(
            method=pMesh['user_subroutine'], fname_input=model.name_job+'.inp')
    
    if not parameters['not_run_job']:
        model.submit_job(name_job)
    