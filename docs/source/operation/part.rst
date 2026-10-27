
Create part
===================================

The :py:class:`Part <AbaqusTools.part.Part>` class contains a 
:py:meth:`create_part <AbaqusTools.part.Part.create_part>` function that creates Abaqus Part objects.

.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.create_part


.. seealso:: 

    The `Abaqus Part <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch48pyo01.html>`_
    object defines the physical attributes of a structure. 
    Parts are instanced into the assembly and positioned before an analysis.


Create datum
-----------------------------------

The Abaqus Datum Planes, Axes, Coordinate Systems (CSYS) are defined in the 
:py:meth:`create_part <AbaqusTools.part.Part.create_part>` function. The following commands are the frequently-used 
`Abaqus Datum commands <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_:

.. code-block:: python
    :linenos:

    # Datum axis
    myPrt.DatumAxisByPrincipalAxis(principalAxis=XAXIS)

    # Datum plane
    myPrt.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)


The :py:class:`Part <AbaqusTools.part.Part>` class provides a 
:py:meth:`create_datum_csys_3p <AbaqusTools.part.Part.create_datum_csys_3p>`
function to create a datum coordinate system by three points 
(, i.e., the origin and two direction vectors along the x and y axis).
The datum CSYS is also given a name for easy access.

.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.create_datum_csys_3p


The :py:class:`Part <AbaqusTools.part.Part>` class provides two functions 
to get access of a datum, either by the creation index 
of a datum (:py:meth:`get_datum_by_index <AbaqusTools.part.Part.get_datum_by_index>`), 
or by the name of datum feature (:py:meth:`get_datum_by_name <AbaqusTools.part.Part.get_datum_by_name>`). 

.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.get_datum_by_index


.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.get_datum_by_name


Create partition
-----------------------------------

The :py:class:`Part <AbaqusTools.part.Part>` class contains a 
:py:meth:`create_partition <AbaqusTools.part.Part.create_partition>` 
function that creates partition of faces and cells.

.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.create_partition


.. code-block:: python
    :linenos:

    myPrt.PartitionCellByExtrudeEdge(line=myPrt.datums[1], cells=myPrt.cells, 
            edges=edges, sense=FORWARD)

    myPrt.PartitionCellByPlaneThreePoints(cells=myPrt.cells, 
            point1=tuple([0, 0, 0]), point2=tuple([0, 1, 0]), point3=tuple([1, 0, 0]))

    myPrt.PartitionFaceBySketch(sketchUpEdge=myPrt.datums[1], faces=faces, 
            sketchOrientation=BOTTOM, sketch=mySkt)


.. seealso:: 

    The `Partition commands <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_
    are used to partition edges, faces, and cells into new regions. 
    A partition command can be invoked for a Part object or for an Assembly object. 
    The partition commands create Feature objects.


Create surface
-----------------------------------

The :py:class:`Part <AbaqusTools.part.Part>` class contains a 
:py:meth:`create_surface <AbaqusTools.part.Part.create_surface>` 
function that creates surfaces in parts.

.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.create_surface


A surface is created by the following Abaqus command:

.. code-block:: 
    :linenos:

    faces = self.get_faces(myPrt, (x, y, z))

    myPrt.Surface(side1Faces=faces, name='face')


Create set
-----------------------------------

The :py:class:`Part <AbaqusTools.part.Part>` class contains a
:py:meth:`create_set <AbaqusTools.part.Part.create_set>`
function that creates sets in parts.

.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.create_set


The :py:class:`Part <AbaqusTools.part.Part>` class proves a 
:py:meth:`create_geometry_set <AbaqusTools.part.Part.create_geometry_set>`  
function to create sets.

.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.create_geometry_set

