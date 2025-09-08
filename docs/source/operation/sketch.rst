
Create sketch
===================================

The :py:class:`Part <AbaqusTools.part.Part>` class contains a 
:py:meth:`create_sketch <AbaqusTools.part.Part.create_sketch>` function that creates stand-alone sketches.

.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.create_sketch


.. seealso:: 

    `Sketches <http://130.149.89.49:2080/v6.13/books/usi/default.htm?startat=pt03ch20.html>`_ 
    are two-dimensional profiles that are used to help form the geometry defining 
    an Abaqus/CAE native part. You use the Sketch module to create a sketch that defines 
    a planar part, a beam, or a partition or to create a sketch that might be extruded, 
    swept, or revolved to form a three-dimensional part. 

    The sketch is usually used to:

    - Create or edit a feature while defining a part in the Part module.
    - Sketch a partition on a face or cell while working in the Part or Assembly module.

    The `ConstrainedSketch <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch48pyo01.html>`_
    object contains the entities that are used to create a sketch.

    The `ConstrainedSketch` objects include 
    `ConstrainedSketchGeometry <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch48pyo01.html>`_ 
    objects contained in the Geometry Repository, 
    such as Line, Arc, and Spline. Vertex, Dimension, Constraint, and Parameter objects 
    are contained in their respective repositories.


The function :py:meth:`create_sketch <AbaqusTools.part.Part.create_sketch>` 
is used to create stand-alone sketches for creating parts and partition faces/cells. 
There are some frequently-used Abaqus Scripting functions to create a sketch:

.. code-block:: python
    :linenos:

    myMdl = Model.model

    mySkt = myMdl.ConstrainedSketch(name='sketch', sheetSize=200.0)

    mySkt.Line(point1=tuple([0.0, 0.0]), point2=tuple([0.0, 1.0]))

    mySkt.CircleByCenterPerimeter(
            center=tuple([0.0, 0.0]), point1=tuple([0.0, 1.0]))
