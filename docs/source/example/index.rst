
.. _examples:

Examples
=================================

This chapter contains examples of building Abaqus models using the 
:py:class:`Part <AbaqusTools.part.Part>` and 
:py:class:`Model <AbaqusTools.model.Model>` classes.

.. note:: 

    When running the example python scripts, 
    put the scripts in the same directory as the 'AbaqusTools' folder.
    Run the following command in the Command Prompt (CMD) of Windows or the terminal of Linux.
    It is suggested to clean up the `*.pyc` files before calling Abaqus, because the `*.pyc` files
    are the compiled Python files that store the bytecode of the source code. 
    They are most likely compiled by Abaqus in previous jobs.
    Abaqus will not update them once they existed, even if the source code is modified.
    The `clean.bat` gives an example of cleaning Abaqus jobs in Windows.

.. code-block:: bash
    :linenos:

    abaqus cae script=script.py

    clean.bat & abaqus cae script=script.py
