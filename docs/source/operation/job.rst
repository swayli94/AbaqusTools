
Create job
===================================

.. seealso:: 

    The `Job commands” 
    <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_ 
    provide methods to create, modify, submit, and control jobs.

A simple example for job:

.. code-block:: python
    :linenos:

    mdb.Job(name='Job-1', model=str(self.name_model), description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
        multiprocessingMode=DEFAULT, numCpus=4, numDomains=4, numGPUs=0)


The :py:class:`Model <AbaqusTools.model.Model>` class provides a 
:py:meth:`submit_job <AbaqusTools.model.Model.submit_job>` function to start the simulation.

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.submit_job


There is also a function to write the `*.inp` file.

.. literalinclude:: ../../../AbaqusTools/model.py
    :language: python
    :linenos: 
    :pyobject: Model.write_job_inp

