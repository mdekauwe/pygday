G'DAY (Generic Decomposition And Yield) Model
=============================================

GDAY simulates carbon, nitrogen and water cycling between the plant and the soil. 

The model is coded entirely in `Python <http://www.python.org/>`_ with only a single dependancy package `configobj <http://www.voidspace.org.uk/python/configobj.html>`_. 


Key Reference
=============
1. Comins, H. N. and McMurtrie, R. E. (1993) Long-Term Response of Nutrient-Limited Forests to CO2 Enrichment; Equilibrium Behavior of Plant-Soil Models. Ecological Applications, 3, 666-681.

.. contents:: :local:

Installation
=============

The G'DAY model depends on very few non-standard python packages, in fact only one I think??!

* `configobj <http://www.voidspace.org.uk/python/configobj.html>`_. (Used for reading cfg/ini files)

On the off chance you have a well behaved setup you should only need to do ::
    
    easy_install configobj

Alternatively follow the simple instructions on the above link. Note, replacing this with a small bit of code to read user input files would be
trivial. I might write that...

Once you have downloaded the source code, or clone the repository (go on...) there is a simple makefile, e.g. ::

    sudo make install

or the standard python approach ::

    sudo python setup.py install

Running the model
=================

I need to add some simple scripts (on the todo list!). But within a python script you would require the following lines... ::
    
    from gday import gday as model
    G = model.Gday(cfg_fname, spin_up=True)
    G.spin_up_pools()

will spin the model up. Spin up expects a met forcing file with a 1000 yrs of data, how you recycle this is up to you. The model automatically stops once the soil, plant and litter pools have reached equilibrium (check code for finer details).

Changing the model default parameters for user defined ones is trivial and utilises a python dictionary, e.g. ::

    from gday import adjust_gday_param_file as ad
    replace_dict = { "albedo": "0.123" }
    ad.adjust_param_file(cfg_fname, replace_dict)

And finally running the model... ::

    from gday import gday as model
    G = model.Gday(cfg_fname)
    G.run_sim()

In all cases ``cfg_fname`` is simple a string with a link to the users parameter file, it can be named anything you please.
    
Contacts
========
* `Martin De Kauwe <http://mdekauwe.github.io/>`_  (mdekauwe at gmail.com)
* `Belinda Medlyn <http://bio.mq.edu.au/people/person.php?user=bmedlyn>`_ (bmedlyn at bio.mq.edu.au).
