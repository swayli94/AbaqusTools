
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


The :py:class:`Model <AbaqusTools.model.Model>` class provides functions 
to setup a step for static analysis or dynamic analysis.

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.create_static_step

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.create_dynamic_step


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
