
Create load
===================================

In the Abaqus Load module, loads and boundary conditions are defined.

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.setup_loads


Load
-----------------------------------

.. seealso:: 

    `Load commands <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_:
    
    - A specific type of load object and a specific type of load state object are 
      designed for each type of load. 
    
    - A load object stores the non-propagating data of a load 
      as well as a number of instances of the corresponding load state object, 
      each of which stores the propagating data of the load in a single step. 
    
    - Instances of the load state object are created and deleted internally by 
      its corresponding load object.

    Examples:

    - The `Pressure object 
      <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_
      defines a pressure load.

    - The `SurfaceTraction object 
      <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_
      defines surface traction on a region.


A simple example for :py:meth:`setup_loads <AbaqusTools.model.Model.setup_loads>`:

.. code-block:: python
    :linenos:

    def setup_loads(self):
    
        region = self.instance('plate_0').sets['z0']
        self.model.DisplacementBC(name='BC-wall-z0', createStepName='Initial', 
            region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
            amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
        
        region = self.instance('plate_0').surfaces['z1']
        self.model.Pressure(name='Load-pressure-z1', createStepName='Loading', 
            region=region, distributionType=UNIFORM, field='', magnitude=1E8, 
            amplitude=UNSET)


Boundary condition
-----------------------------------

.. seealso:: 

    `Boundary Condition commands 
    <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_:

    A specific type of boundary condition object and a specific type of boundary condition 
    state object are designed for each type of boundary condition. 

    A BoundaryCondition object stores the non-propagating data of a boundary condition 
    as well as a number of instances of the corresponding BoundaryConditionState object, 
    each of which stores the propagating data of the boundary condition in a single step. 
    
    Instances of the BoundaryConditionState object are created and deleted internally by 
    its corresponding BoundaryCondition object.


Reference point
-----------------------------------

Reference point is useful in Abaqus/CAE.

The :py:class:`Model <AbaqusTools.model.Model>` class provides functions to create and access reference points.

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.create_reference_point


.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.get_reference_point


.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.create_reference_point_set



