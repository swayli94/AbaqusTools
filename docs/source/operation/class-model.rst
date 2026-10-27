
.. _operations_create_model:

Create a model with class `Model`
===================================

The :py:class:`Model <AbaqusTools.model.Model>` class is used to create an Abaqus model, 
an empty Abaqus model is created during the creation of an instance of 
:py:class:`Model <AbaqusTools.model.Model>` object.
The :py:meth:`__init__ <AbaqusTools.model.Model.__init__>` function calls the 
:py:meth:`initialization <AbaqusTools.model.Model.initialization>`
function to prepare class attributes, create model materials and sections.

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.initialization


An Abaqus Model is built by calling the 
:py:meth:`build <AbaqusTools.model.Model.build>` function, which contains a series of functions:

- :py:meth:`setup_parts <AbaqusTools.model.Model.setup_parts>` 
  builds different parts with their :py:class:`Part <AbaqusTools.part.Part>` objects;

- :py:meth:`setup_assembly <AbaqusTools.model.Model.setup_assembly>` 
  assembles the parts;
  
- :py:meth:`setup_steps <AbaqusTools.model.Model.setup_steps>` 
  creates a static/dynamic analysis step;
  
- :py:meth:`setup_interactions <AbaqusTools.model.Model.setup_interactions>`
  defines interactions between different parts;

- :py:meth:`setup_loads <AbaqusTools.model.Model.setup_loads>`
  defines loads and boundary conditions.
  
- :py:meth:`setup_outputs <AbaqusTools.model.Model.setup_outputs>`
  define outputs;

- :py:meth:`setup_jobs <AbaqusTools.model.Model.setup_jobs>`
  define jobs;


.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.build



