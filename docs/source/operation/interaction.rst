
Create interaction
===================================

The :py:class:`Model <AbaqusTools.model.Model>` class contains a 
:py:meth:`setup_interactions <AbaqusTools.model.Model.setup_interactions>` 
function that contains functions to define interactions.

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.setup_interactions


The :py:class:`Model <AbaqusTools.model.Model>` class provides functions to define contacts.

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.create_interaction_property_contact


.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.create_contact_constraints


.. seealso:: 

    `Interaction commands 
    <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_ :

    - A specific type of interaction object and a specific type of interaction state object 
      are designed for each type of interaction. 
    
    - An interaction object stores the non-propagating data of an interaction as well as 
      a number of instances of the corresponding interaction state object, 
      each of which stores the propagating data of the interaction in a single step. 
    
    - Instances of the interaction state object are created and deleted internally 
      by its corresponding interaction object.

    
