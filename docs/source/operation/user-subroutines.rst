
User subroutines
=================================

User subroutines allow advanced users to customize a wide variety of Abaqus capabilities. 
Information on writing user subroutines and detailed descriptions of each subroutine appear online 
in the `Abaqus User Subroutines Guide <https://abaqus-docs.mit.edu/2017/English/SIMACAEANLRefMap/simaanl-c-subroutineover.htm>`_ . 
A listing and explanations of associated utility routines also appear in that guide.

User subroutines:

- are provided to increase the functionality of several Abaqus capabilities for which the usual data input methods alone may be too restrictive;
- provide an extremely powerful and flexible tool for analysis;
- are written as C, C++, or Fortran code and must be included in a model when you execute the analysis;
- must be included and, if desired, can be revised in a restarted run, since they are not saved to the restart files;
- cannot be called one from another; and
- can in some cases call utility routines that are also available in Abaqus.

In this project, the :ref:`umat_introduction`, :ref:`uvarm_introduction`, and :ref:`usdfld_introduction` are used.


.. _umat_introduction :

UMAT
----------------------------------

The `user-defined mechanical material behavior <https://abaqus-docs.mit.edu/2017/English/SIMACAEMATRefMap/simamat-c-usermat.htm>`_ (UMAT/VUMAT) 

- is provided by means of an interface whereby any mechanical constitutive model can be added to the library;
- requires that a constitutive model (or a library of models) is programmed in user subroutine UMAT (Abaqus/Standard) or VUMAT (Abaqus/Explicit); and
- requires considerable effort and expertise: the feature is very general and powerful, but its use is not a routine exercise.

.. note:: 

    The use of this subroutine generally requires considerable expertise. 
    You are cautioned that the implementation of any realistic constitutive model requires extensive development and testing. 
    Initial testing on a single-element model with prescribed traction loading is strongly recommended.


The user subroutine `UMAT <https://abaqus-docs.mit.edu/2017/English/SIMACAESUBRefMap/simasub-c-umat.htm#simasub-c-umat>`_:

- can be used to define the mechanical constitutive behavior of a material;
- will be called at all material calculation points of elements for which the material definition includes a user-defined material behavior;
- can be used with any procedure that includes mechanical behavior;
- can use solution-dependent state variables (SDV);
- must update the stresses and solution-dependent state variables to their values at the end of the increment for which it is called;
- must provide the material Jacobian matrix, :math:`\partial \Delta \sigma / \partial \Delta \varepsilon` , for the mechanical constitutive model;
- can be used in conjunction with user subroutine :ref:`usdfld_introduction` to redefine any field variables before they are passed in; and
- is described further in `user-defined mechanical material behavior <https://abaqus-docs.mit.edu/2017/English/SIMACAEMATRefMap/simamat-c-usermat.htm>`_.


.. _uvarm_introduction :

UVARM
----------------------------------

The `user subroutine UVARM <https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/sub/default.htm?startat=ch01s01asb44.html>`_
generates element output, which is used to define some special variables that Aabaqus itself does not have, such as: 
damage factor, risk coefficient, safety margin, etc. This subroutine only works with the Abaqus/Standard solver.

The user subroutine `UVARM <https://abaqus-docs.mit.edu/2017/English/SIMACAESUBRefMap/simasub-c-uvarm.htm#simasub-c-uvarm>`_:

- will be called at all material calculation points of elements for which the material definition 
  includes the specification of user-defined output variables;
- may be called multiple times for each material point in an increment, as Abaqus/Standard iterates to a converged solution;
- will be called for each increment in a step;
- allows you to define output quantities that are functions of any of the available integration point quantities listed 
  in the Output Variable Identifiers table (Abaqus/Standard output variable identifiers);
- allows you to define the material directions as output variables;
- can be used for gasket elements;
- can call utility routine GETVRM to access material point data;
- cannot be used with linear perturbation procedures; and
- cannot be updated in the zero increment.

.. note:: 

    The UVARM subroutine can obtain the stress, strain and other information on the material integration point in the current analysis step. 
    Users can use this information to define some output variables that abaqus itself does not have. 

    Taking three-dimensional composite laminates as an example, the two-dimensional Hashin's failure model that comes with Abaqus is not applicable. 
    The initial failure criteria such as Cai Wu, Cai Xier, and maximum stress and maximum strain can only be used for plane stress 
    and plane strain problems, and there are not applicable as well for the three-dimensional solid composite laminates. 

    In this case, users can define a three-dimensional failure criterion by themselves to determine which locations in the structure are safe, 
    which locations are dangerous, what the specific risk coefficient is, etc.

.. note:: 

    It should be stated that UVARM can only output some custom variables and cannot change the material constitutive relationship 
    and the original constitutive parameters. 

    If you want to change the parameter values of Abaqus's own constitutive, you can use the :ref:`usdfld_introduction`/VUSDFLD subroutine. 

    If you want to change the constitutive relationship, you need to write :ref:`umat_introduction`/VUMAT subroutines, and the difficulty increases in sequence.

The UVARM subroutine is more suitable for structural engineering analysis, or early strength prediction of structural plans, 
and its scope of application is still very wide.


.. _usdfld_introduction :

USDFLD
----------------------------------

The `user subroutine USDFLD <https://abaqus-docs.mit.edu/2017/English/SIMACAESUBRefMap/simasub-c-usdfld.htm#simasub-c-usdfld>`_
redefines field variables at a material point. It

- allows you to define field variables at a material point as functions of time or of any of the available material point quantities 
  listed in the Output Variable Identifiers table (Abaqus/Standard output variable identifiers) except the user-defined output variables UVARM and UVARMn;
- can be used to introduce solution-dependent material properties since such properties can easily be defined as functions of field variables;
- will be called at all material points of elements for which the material definition includes user-defined field variables;
- must call utility routine GETVRM to access material point data;
- can use and update state variables; and
- can be used in conjunction with user subroutine UFIELD to prescribe predefined field variables.

USDFLD enables user to defined custom field variables. All functions of USDFLD can be achieved using :ref:`umat_introduction`, 
but its application is relatively simpler than UMAT. 
Users can use the constitutive model that comes with Abaqus, instead of re-developing  a new material constitutive model. 
During the simulation, USDFLD reads the field variables at the integration point,
then, calculates and uploads the new custom field variables to Abaqus.
In this way, users effectively modify the built-in constitutive model.



