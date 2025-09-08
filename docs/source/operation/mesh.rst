
.. _operation_create_mesh:

Create mesh
===================================

In the :py:meth:`build <AbaqusTools.model.Model.build>` function, 
three functions are called to generate mesh for a part, i.e., 
:py:meth:`set_seeding <AbaqusTools.part.Part.set_seeding>`, 
:py:meth:`create_mesh <AbaqusTools.part.Part.create_mesh>`, and 
:py:meth:`set_element_type <AbaqusTools.part.Part.set_element_type>`.


Set seeding
-----------------------------------

The :py:class:`Part <AbaqusTools.part.Part>` class contains a 
:py:meth:`set_seeding <AbaqusTools.part.Part.set_seeding>` 
function that seeds edges, regions and parts.

.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.set_seeding


The following commands are the
`Abaqus Mesh commands <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_
of Abaqus Part objects for seeding:

.. code-block:: python
    :linenos:

    # Seed edge
    myPrt.seedEdgeByBias(biasMethod=SINGLE, end1Edges=edges, 
                    minSize=minSize, maxSize=maxSize, constraint=FINER)
                
    myPrt.seedEdgeByBias(biasMethod=DOUBLE, end1Edges=edges, 
                    minSize=minSize, maxSize=maxSize, constraint=FIXED)

    myPrt.seedEdgeBySize(edges=edges, size=size, 
                    deviationFactor=0.1, minSizeFactor=0.1, constraint=FIXED)

    myPrt.seedEdgeByNumber(edges=edges, number=number, constraint=FIXED)

    # Seed part
    myPrt.seedPart(size=size, deviationFactor=0.1, minSizeFactor=0.1)


Generate mesh
-----------------------------------

The :py:class:`Part <AbaqusTools.part.Part>` class contains a 
:py:meth:`create_mesh <AbaqusTools.part.Part.create_mesh>`
function that sets mesh control and generates mesh.

.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.create_mesh


The following commands are the
`Abaqus Mesh commands <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_
of Abaqus Part objects for mesh control and generation:

.. code-block:: python
    :linenos:

    # Mesh control
    myPrt.setMeshControls(regions=myPrt.cells, technique=STRUCTURED)
    myPrt.setMeshControls(regions=myPrt.cells, technique=SWEEP, algorithm=ADVANCING_FRONT)

    # Set sweep path
    myPrt.setSweepPath(region=cell, edge=edge, sense=FORWARD)

    # Assign stack direction
    myPrt.assignStackDirection(referenceRegion=face, cells=cells)

    # Generate mesh
    myPrt.generateMesh()


Set element type
-----------------------------------

The :py:class:`Part <AbaqusTools.part.Part>` class contains a 
:py:meth:`set_element_type <AbaqusTools.part.Part.set_element_type>` 
function that assigns element types to the specified regions.

.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.set_element_type


The following command is the
`Abaqus Mesh commands <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_
of Abaqus Part objects for assigning element types:


.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.set_element_type_of_part

