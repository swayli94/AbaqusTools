
Basic operations
=================================

In this `AbaqusTools` Abaqus Python Script, an Abaqus model is built by carrying out operations 
in the following **Modules** of Abaqus CAE. 
The operations in the **Sketch**, **Part**, and **Mesh** modules are carried out within 
the :py:class:`Part <AbaqusTools.part.Part>` class, the operations in other **Modules** 
are carried out within the 
:py:class:`Model <AbaqusTools.model.Model>` class.

1)  **Property**;
2)  **Sketch**;
3)  **Part**;
4)  **Mesh**;
5)  **Assembly**;
6)  **Step**;
7)  **Interaction**;
8)  **Load**;
9)  **Job**;


.. toctree::
    :maxdepth: 2
    :caption: Basic operations:

    abaqus
    class-model
    class-part
    property
    sketch
    part
    mesh
    assembly
    step
    interaction
    load
    job
    pbc
    user-subroutines


.. seealso:: 

    `Abaqus/CAE User's Guide (2016)
    <http://130.149.89.49:2080/v6.13/books/usi/default.htm?startat=pt03ch20s01.html>`_

    `Abaqus Scripting Reference Guide (2016)
    <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch48pyo01.html>`_
