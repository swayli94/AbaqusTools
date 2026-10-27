
Create assembly
===================================

The :py:class:`Model <AbaqusTools.model.Model>` class contains a 
:py:meth:`setup_assembly <AbaqusTools.model.Model.setup_assembly>` function that assembles Abaqus Parts.


.. seealso:: 

    An `Assembly object <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_
    is a container for instances of Abaqus Parts. 
    The Assembly object has no constructor command. 
    Abaqus creates the rootAssembly member when an Abaqus Model object is created.

    There are `Abaqus commands <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_ 
    operate on Abaqus Model objects.
    For example, the Abaqus `Instance()` function copies an Abaqus PartInstance object from 
    the specified Abaqus Model and creates a new Abaqus PartInstance object.


.. code-block:: python
    :linenos:

    a = myMdl.rootAssembly
    p = myMdl.parts['name_part']
    a.Instance(name='name_instance', part=p, dependent=ON)

    a.translate(instanceList=('name_instance',), vector=tuple(0,0,1))


.. seealso:: 

    Merge multiple Abaqus Parts into a new Abaqus Part:

    The `InstanceFromBooleanMerge <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_
    function of `PartInstance` object creates an Abaqus PartInstance in the instances repository 
    after merging two or more Abaqus Part instances.


.. code-block:: python
    :linenos:

    a = myMdl.rootAssembly

    a.InstanceFromBooleanMerge(name='name_new_part', 
    instances=(a.instances['part-1'], a.instances['part-2']), 
    keepIntersections=ON, originalInstances=DELETE, domain=GEOMETRY)

    del a.features['complete_beam-1']

