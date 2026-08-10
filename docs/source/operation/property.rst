
.. _Create_property:

Create property
===================================

Usually, we need to first create *Material* and *Section* for an Abaqus model.
We also need to assign *Section* or create *Composite Layup* for different parts or cells,
after they are created.

The :py:class:`Model <AbaqusTools.model.Model>` class contains a 
:py:meth:`setup_property <AbaqusTools.model.Model.setup_property>` 
function that creates material and sections, 
these functions are called during the initialization of an Abaqus model, 
i.e., the creation of an instance of :py:class:`Model <AbaqusTools.model.Model>` object.

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.setup_property

Material data is kept separately from the Abaqus API calls that consume it.
The material cards live in :py:mod:`AbaqusTools.materials` as plain dictionaries,
stored in the N-mm-tonne unit system, and
:py:meth:`create_material <AbaqusTools.model.Model.create_material>` turns an
entry of that library into an Abaqus material:

.. code-block:: python
    :linenos:

    def setup_property(self):

        self.create_material('IM7/8551-7', elastic_type='ENGINEERING_CONSTANTS')
        self.create_section('IM7/8551-7')

        self.create_material('Ti-6Al-4V')
        self.create_section('Ti-6Al-4V')

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos:
    :pyobject: Model.create_material

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos:
    :pyobject: Model.create_section

Adding a material is a matter of adding an entry to
``AbaqusTools.materials.MATERIAL_LIBRARY``, not of adding a method. Each entry
declares its tables together with the indices of the columns that carry a stress
unit, so that ``unit_length='m'`` can convert them:

.. literalinclude:: ../../../AbaqusTools/materials.py
    :language: python
    :linenos:
    :start-after: MATERIAL_LIBRARY = {
    :end-before: 'Ti-6Al-4V': {

Because :py:mod:`AbaqusTools.materials` imports nothing from Abaqus, it can also
be read by post-processing scripts running in a normal Python interpreter.


The section assignment and composite layup creation are defined in the :py:class:`Part <AbaqusTools.part.Part>` class,
these operations are carried out during the creation of parts.

.. code-block:: python
    :linenos: 

    class Part(object):

        def set_section_assignment(self):
            '''
            Set section assignment via Abaqus Module: Property
            
            - Assign Section
            '''
            myPrt = self.model.parts[self.name_part]
            
            myPrt.SectionAssignment(region=myPrt.sets['all'], sectionName='Steel', offset=0.0, 
                offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

        def set_composite_layups(self):
            '''
            Set composite layup via Abaqus Module: Property
            
            - Create Composite Layup
            '''
            myPrt = self.model.parts[self.name_part]

            self.set_CompositeLayup_of_set(myPrt, 
                    name_set=           'all', 
                    total_thickness=    1.0, 
                    ply_angle=          [0, 45, -45, 90, 0],
                    symmetric=          True, 
                    numIntPoints=       3)

The built-in function 
:py:meth:`set_CompositeLayup_of_set <AbaqusTools.part.Part.set_CompositeLayup_of_set>` 
of class :py:class:`Part <AbaqusTools.part.Part>` creates the composite layup for a set.

.. literalinclude:: ../../../AbaqusTools/part.py
    :language: python
    :linenos: 
    :pyobject: Part.set_CompositeLayup_of_set


.. seealso:: 

    A `section <https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/usi/default.htm?startat=pt03ch12s02s03.html>`_
    contains information about the properties of a part or a region of a part.
    You can use the Property module to create solid sections, shell sections, 
    beam sections, and other sections.

    `Composite layups <http://dsk-016-1.fsid.cvut.cz:2080/v6.12/books/usi/default.htm?startat=pt04ch23s01.html>`_ 
    in Abaqus/CAE are designed to help you manage a large number of plies 
    in a typical composite model. In contrast, composite sections are a product of 
    finite element analysis and may be difficult to apply to a real-world application.



