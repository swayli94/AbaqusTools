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


SET_NAME = [ ['PLATE',     ['ALL']] ]

if __name__ == '__main__':

    t0 = time.time()
    
    name_job = 'Job_OHT'

    odb = OdbOperation(name_job)
    
    def get_element_value_on_set(name_instance, name_set):
    
        node_labels, indices_fieldOutput = odb.get_node_labels_and_indices(name_instance, name_set=name_set)
        
        coordinates = odb.probe_node_coordinate(name_instance=name_instance, node_label=node_labels)
        values_U    = odb.probe_node_values(variable='U', index_fieldOutput=indices_fieldOutput)

        return indices_fieldOutput, coordinates, values_U

    f0 = open('specimen-U.dat', 'w')
    f0.write('Variables= X Y Z U1 U2 U3 index \n')

    for name_instance, name_sets in SET_NAME:
        for name_set in name_sets:
            
            t1 = time.time()
            name_zone = 'zone T=" %s %s "'%(name_instance, name_set)
            
            print('>>> ==========================================')
            print('>>> '+name_zone)
            
            
            indices_fieldOutput, coordinates, values_U = \
                get_element_value_on_set(name_instance, name_set)

            n_element = len(indices_fieldOutput)

            f0.write(name_zone+' \n')
            for i in range(n_element):

                for j in range(3):
                    f0.write(' %14.6E'%(coordinates[i][j]))
                
                for j in range(3):
                    f0.write(' %14.6E'%(values_U[i][j]))
                    
                f0.write(' %d \n'%(indices_fieldOutput[i]))
            f0.write('\n')
            
            t2 = time.time()
            print('>>> Number of element: %d, Time = %.2f s'%(n_element, (t2-t1)))
            print(' ')
        
    f0.close()



