
.. _operations_create_part:

Create a part with class `Part`
===================================

The :py:class:`Part <AbaqusTools.part.Part>` class is used to create an Abaqus Part for assembly, 
which includes creating sketches, partitions, surfaces, sets and mesh, as well as assigning sections
and editing composite layups.

An Abaqus Part is built by calling the 
:py:meth:`build <AbaqusTools.part.Part.build>` function, which contains a series of functions:

- :py:meth:`create_sketch <AbaqusTools.part.Part.create_sketch>` 
  creates stand-alone sketches via Abaqus Module: Sketch;

- :py:meth:`create_part <AbaqusTools.part.Part.create_part>` 
  creates Abaqus Part via Abaqus Module: Part, including:

  - Create Datum Plane/Axis/CSYS;
  - Create Part;
  - Create Solid: Extrude;
  - Create Shell: Extrude;
  - Create Cut: Extrude;
  - Create Round and Fillet;
  - Create Mirror;
  
- :py:meth:`create_partition <AbaqusTools.part.Part.create_partition>` 
  creates partition via Abaqus Module: Part, including:

  - Create Partition;
  - Partition Face: Sketch;
  - Partition Face: Use Datum Plane;
  - Partition Cell: Define Cutting Plane;
  - Partition Cell: Use Datum Plane;
  - Partition Cell: Extend Face;
  - Partition Cell: Extrude/Sweep Edges;
  
- :py:meth:`create_surface <AbaqusTools.part.Part.create_surface>`;
- :py:meth:`create_set <AbaqusTools.part.Part.create_set>`;

- :py:meth:`set_seeding <AbaqusTools.part.Part.set_seeding>`
  sets seeding via Abaqus Module: Mesh, including:

  - Seed Edges;
  - Seed Part;
  
- :py:meth:`create_mesh <AbaqusTools.part.Part.create_mesh>`
  sets mesh control and creates mesh via Abaqus Module: Mesh, including:

  - Assign Mesh Control;
  - Assign Stack Direction;
  - Mesh Region;
  - Mesh Part;

- :py:meth:`set_element_type <AbaqusTools.part.Part.set_element_type>`
  sets element type via Abaqus Module: Mesh;

- :py:meth:`set_section_assignment <AbaqusTools.part.Part.set_section_assignment>`
  sets section assignment via Abaqus Module: Property;

- :py:meth:`set_composite_layups <AbaqusTools.part.Part.set_composite_layups>`
  sets composite layup via Abaqus Module: Property;


.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.build


