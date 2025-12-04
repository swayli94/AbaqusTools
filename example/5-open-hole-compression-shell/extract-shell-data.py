'''
Extract data from the *.odb file.
'''
import time

from AbaqusTools import OdbOperation

try:

    from abaqus import *
    from abaqusConstants import *
    from caeModules import *

except:
    
    pass


# SET_NAME = [ ['PLATE',     ['ALL']] ]
SET_NAME = [ ['PLATE',     ['PARTITION_CIRCLE']] ]


if __name__ == '__main__':

    t0 = time.time()
    
    name_job = 'Job_OHT'
    fname_mean = 'specimen-stress-field-mean.dat'
    fname_3D = 'specimen-stress-field-3D.dat'

    odb = OdbOperation(name_job)
    
    def get_element_values_on_set(name_instance, name_set):
        '''
        Get the variable value of elements in a set of an instance.
        
        Parameters
        --------------
        name_instance: str
            name of the instance in CAPITAL letters, e.g., 'PLATE'.
            
        name_set: str
            name of the set in CAPITAL letters, e.g., 'FACE_HOLE'.
            
        Returns
        --------------
        indices_fieldOutput: list[int]
            indices of the elements in the fieldOutput.
            
        coordinates: ndarray [n_element, n_dim]
            the coordinates of the elements.
            
        data: Dict[str, ndarray]
            the data of the elements.
            The key is the variable name, e.g., 'S11', 'S22', 'S12'.
            The value is the value of the variable, i.e., ndarray [n_element] or [n_element, n_thickness].
            The 'thickness_distribution' is the thickness coordinates of the thickness-direction-distributed data,
            ndarray [n_thickness].
        '''
        element_labels, indices_fieldOutput = odb.get_element_labels_and_indices(name_instance, name_set)
        
        # Element-wise coordinate: ndarray [n_element, n_dim], i.e., [(x, y, z), ...].
        coordinates = odb.probe_element_center_coordinate(name_instance=name_instance, element_label=element_labels)
        
        # Element-wise data: ndarray [n_element, 1], i.e., [(comp1,), ...].
        values_S11 = odb.probe_element_values(variable='S', component='S11', index_fieldOutput=indices_fieldOutput)
        values_S22 = odb.probe_element_values(variable='S', component='S22', index_fieldOutput=indices_fieldOutput)
        values_S12 = odb.probe_element_values(variable='S', component='S12', index_fieldOutput=indices_fieldOutput)
        
        # Thickness-direction-distributed data: ndarray [n_element, n_thickness, 2], i.e., [[(coordinate, value), ...], ...].
        values_thickness_S11 = odb.probe_shell_element_thickness_values(variable='S', component='S11', name_instance=name_instance, element_label=element_labels)
        values_thickness_S22 = odb.probe_shell_element_thickness_values(variable='S', component='S22', name_instance=name_instance, element_label=element_labels)
        values_thickness_S12 = odb.probe_shell_element_thickness_values(variable='S', component='S12', name_instance=name_instance, element_label=element_labels)
        
        data = {
            'S11': values_S11[:,0],
            'S22': values_S22[:,0],
            'S12': values_S12[:,0],
            'thickness_S11': values_thickness_S11[:,:,1],
            'thickness_S22': values_thickness_S22[:,:,1],
            'thickness_S12': values_thickness_S12[:,:,1],
            'thickness_distribution': values_thickness_S11[0,:,0]
        }
        
        return indices_fieldOutput, coordinates, data


    f = open(fname_mean, 'w')
    f.write('Variables= X Y Z S11 S22 S12 index\n')
    f.close()
    
    f = open(fname_3D, 'w')
    f.write('Variables= X Y Z thickness S11 S22 S12 index index_thickness\n')
    f.close()

    for name_instance, name_sets in SET_NAME:
        for name_set in name_sets:
            
            t1 = time.time()
            name_zone = 'zone T=" %s %s "'%(name_instance, name_set)
            
            print('>>> ==========================================')
            print('>>> '+name_zone)
            
            indices_fieldOutput, coordinates, data = \
                get_element_values_on_set(name_instance, name_set)

            n_element = len(indices_fieldOutput)
            
            thickness_distribution = data['thickness_distribution']
            n_thickness = len(thickness_distribution)
            
            #* Write mean stress field
            with open(fname_mean, 'a') as f:
                f.write(name_zone+' \n')
                for i in range(n_element):
                    for j in range(3):
                        f.write(' %14.6E'%(coordinates[i][j]))
                    f.write(' %14.6E'%(data['S11'][i]))
                    f.write(' %14.6E'%(data['S22'][i]))
                    f.write(' %14.6E'%(data['S12'][i]))
                    f.write(' %d'%(indices_fieldOutput[i]))
                    f.write('\n')
                f.write('\n')

            #* Write 3D stress field
            with open(fname_3D, 'a') as f:
                f.write(name_zone+' \n')
                for i in range(n_element):
                    for k in range(n_thickness):
                        for j in range(3):
                            f.write(' %14.6E'%(coordinates[i][j]))
                        f.write(' %14.6E'%(thickness_distribution[k]))
                        f.write(' %14.6E'%(data['thickness_S11'][i,k]))
                        f.write(' %14.6E'%(data['thickness_S22'][i,k]))
                        f.write(' %14.6E'%(data['thickness_S12'][i,k]))
                        f.write(' %d'%(indices_fieldOutput[i]))
                        f.write(' %d'%(k))
                        f.write('\n')
            
            t2 = time.time()
            print('>>> Number of element: %d, Time = %.2f s'%(n_element, (t2-t1)))
            print(' ')
        
    



