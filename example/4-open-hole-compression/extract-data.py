'''
Extract data in the hole face from the *.odb file.
'''

import os
import time
import numpy as np

from AbaqusTools import OdbOperation

try:

    from abaqus import *
    from abaqusConstants import *
    from caeModules import *

except:
    
    pass


SET_NAME = [ ['PLATE',     ['PARTITION_CIRCLE']] ]


if __name__ == '__main__':

    t0 = time.time()
    
    name_job = 'Job_OHT'

    odb = OdbOperation(name_job)
    
    def get_element_value_on_set(name_instance, name_set):
    
        element_labels, indices_fieldOutput = odb.get_element_labels_and_indices(name_instance, name_set)

        coordinates = odb.probe_element_center_coordinate(name_instance=name_instance, element_label=element_labels)
        values_S    = odb.probe_element_values(variable='S', index_fieldOutput=indices_fieldOutput)

        return indices_fieldOutput, coordinates, values_S

    f0 = open('specimen-stress-field.dat', 'w')
    f0.write('Variables= X Y Z S11 S22 S33 S12 S13 S23 index\n')

    for name_instance, name_sets in SET_NAME:
        for name_set in name_sets:
            
            t1 = time.time()
            name_zone = 'zone T=" %s %s "'%(name_instance, name_set)
            
            print('>>> ==========================================')
            print('>>> '+name_zone)
            
            indices_fieldOutput, coordinates, values_S = \
                get_element_value_on_set(name_instance, name_set)

            n_element = len(indices_fieldOutput)

            f0.write(name_zone+' \n')
            for i in range(n_element):

                for j in range(3):
                    f0.write(' %14.6E'%(coordinates[i][j]))
                
                for j in range(6):
                    f0.write(' %14.6E'%(values_S[i][j]))
                    
                f0.write(' %d \n'%(indices_fieldOutput[i]))
            f0.write('\n')
            
            t2 = time.time()
            print('>>> Number of element: %d, Time = %.2f min'%(n_element, (t2-t1)/60.0))
            print(' ')
        
    f0.close()



