
Create step
===================================

Setup step
-----------------------------------

.. seealso:: 

    Chapter 49, `“Step commands (step),” 
    <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_
    described in this chapter are used to create and configure analysis steps. 

    Chapter 50, `“Step commands (miscellaneous),” 
    <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_ 
    describes the commands used to configure controls, damping, and frequency tables.


The :py:class:`Model <AbaqusTools.model.Model>` class provides a function
to setup a step for static analysis. It pins the conventions used throughout
this package: the step is named ``'Loading'`` and follows ``'Initial'``.

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos:
    :pyobject: Model.create_static_step

Other step types are created directly with the Abaqus API, since they carry no
convention worth wrapping. For an explicit dynamic analysis, note that the
element types must then be selected from the ``EXPLICIT`` element library:

.. code-block:: python
    :linenos:

    self.model.ExplicitDynamicsStep(
            name=               'Loading',
            previous=           'Initial',
            description=        'Dynamic (explicit) simulation',
            nlgeom=             ON,
            improvedDtMethod=   ON)

A linear perturbation buckling step is built in
``wingbox/wingbox_model.py``, which assembles the ``BuckleStep`` options from
the run parameters and supports a preloaded base state.


Setup output
-----------------------------------

.. seealso:: 

    Chapter 51, `“Step commands (output),” 
    <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_ 
    describes the commands used to create and configure output requests and 
    integrated output sections and the commands to configure diagnostic printing, 
    monitoring, and restart. 


A simple example for field output:

.. code-block:: python
    :linenos:

    myMdl.fieldOutputRequests['F-Output-1'].setValues(
        variables=('S', 'E', 'U', 'RF', 'CF'), 
        frequency=1)
