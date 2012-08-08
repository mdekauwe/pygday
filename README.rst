====================
G'DAY
====================

.. contents:: :local:

Installation
=============

G'DAY model depends on very few other packages, in fact only one I think??

http://www.voidspace.org.uk/python/configobj.html

And replacing this with a small bit of code to read user input files would be
trivial.

Once you have downloaded the source code, or cloned the repository there is 
a simple makefile, e.g.

    sudo make install

or

    sudo python setup.py install

I need to add some simple scripts (to do list!). But...

    from gday import gday as model
    G = model.Gday(cfg_fname)
    G.run_sim()

(Also spin up option, adjusting parameter file...add!).

