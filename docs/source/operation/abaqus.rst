
Abaqus Scripting
===================================

.. seealso:: 

    `Abaqus/CAE User's Guide (2016)
    <http://130.149.89.49:2080/v6.13/books/usi/default.htm?startat=pt03ch20s01.html>`_

    `Abaqus Scripting Reference Guide (2016)
    <http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch48pyo01.html>`_


Model database (Mdb)
-----------------------------------

The `Mdb object 
<http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_ 
is the high-level Abaqus model database. 
A model database stores models and analysis controls.


Features
-----------------------------------

Abaqus/CAE is a feature-based modeling system, and features are stored in the Feature object. 

Features in Abaqus/CAE include **Parts**, **Datums**, **Partitions**, and **Assembly operations**.

- The Part object defines the physical attributes of a structure. 
  Parts are instanced into the assembly and positioned before an analysis.

- The Datum commands return Feature objects and inherit the methods of Feature objects. 
  Datums can be created using methods on a Part or Assembly object.
  Each command also creates a Datum object in the corresponding datum repository. 
  The Datum object is used as an argument to other commands, such as Part and Partition commands.

- The Partition commands are used to partition edges, faces, and cells into new regions. 
  A partition command can be invoked for a Part object or for an Assembly object. 
  The partition commands create Feature objects.

- An Assembly object is a container for instances of parts. 
  The Assembly object has no constructor command. 
  Abaqus creates the rootAssembly member when a Model object is created.

Feature objects contain both the parameters and the resulting model modification.
Abaqus/CAE modifies the model based on the value of the parameters of features. 
This evaluation of the parameters is called regeneration of the feature. 


Regions
-----------------------------------

`Region commands 
<http://130.149.89.49:2080/v2016/books/ker/default.htm?startat=pt01ch20pyo01.html>`_ 
are used to create part and assembly sets or surfaces from elements, nodes, and geometry. 
For more information, see `“Specifying a region” 
<http://130.149.89.49:2080/v2016/books/cmd/default.htm?startat=pt02ch06s06.html#cmd-int-acl-regions>`_. 

Part and assembly objects have the following member, a repository of Set objects: `sets`.
In turn, a Set object can contain any one of the following types:
`elements`, `nodes`, and `geometry`.

A Set object can contain a number of entities of a single type (nodes, elements, or geometry) 
or a combination of node and element types. 
However, except for nodes and elements, a Set object cannot contain a combination of types.
The following are members of the Set object:
`elements`, `nodes`, `cells`, `faces`, `edges`, `vertices`, and `referencePoints`.

Region commands are also used to create surfaces on the assembly. 
Surfaces are sets with additional “sidedness” information.
Part sets contain regions of a part. 

You can assign section definitions to a set created by selecting a region of a part. 
The part sets can be accessed from the instance; however, the section definition you assigned 
to a region is copied automatically to all instances of that part in the assembly.

Assembly sets contain regions of an assembly and are used by many commands that 
operate on the assembly. 
For example, you can apply a load or boundary condition to a set created by 
selecting a region of the assembly. 
Sets can include regions from multiple part instances.

