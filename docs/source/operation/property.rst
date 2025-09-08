
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

The following code blocks show examples of the creation of materials and sections,
these examples are built-in functions of the :py:class:`Model <AbaqusTools.model.Model>` class.

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.create_material_IM785517

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.create_material_steel

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.create_section_steel


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



