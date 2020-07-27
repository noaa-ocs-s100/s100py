Docs need a python 3 environment with sphinx (2.4.4 or 3.0.4), sphinx-automodapi (0.12), sphinx-autodoc-typehints, sphinxcontrib-fulltoc and python-graphviz (0.13.2) which will automatically get graphviz too (really you only need graphviz).
Also there is a modification for sphinx-automod 0.12 and sphinx itself
In addition to the file modifications there is a class template in the docs directory.

The purpose of the setup is to make the docs come at with a full class per html page.  Otherwise the docs either have one function per page or one file per page.  
I prefer docs with each class fully described so I can see the other methods available yet not the other classes to confuse me.
Automodapi effects the one class per page.  The source modifications fix two annoying behaviors of sphinx/automodapi with the diagrams.  
First the inheritance diagrams would go back too far, so setting top classes can fix this.  
Second the diagrams would always be in one direction (horizontal or vertical) but they are more effective when the file/module diagrams are horizontal but the 
class page is shown vertical.

Most of the magic is in the conf.py in this directory to set up automodapi

conda install sphinx sphinx-automodapi sphinx-autodoc-typehints python-graphviz sphinxcontrib-fulltoc 
  or (the current install at time of writing)
conda install sphinx=2.4.4 sphinx-automodapi=0.12 python-graphviz=0.13.2 sphinx-autodoc-typehints=1.10.3 graphviz=2.38 sphinxcontrib-fulltoc=1.2.0
  or
conda install sphinx=3.0.4 sphinx-automodapi=0.12 python-graphviz=0.13.2 sphinx-autodoc-typehints=1.10.3 graphviz=2.38 sphinxcontrib-fulltoc=1.2.0
  then
copy "PydroXL_19\NOAA\site-packages\Python3\s100py\docs\automodsumm_0.12.py" "PydroXL_19\envs\Pydro367\site-packages\sphinx_automodapi\automodsumm.py" 
  (replace x.x.x with 2.4.4 or 3.0.4 as needed)
copy "PydroXL_19\NOAA\site-packages\Python3\s100py\docs\inheritance_diagram_x.x.x_revised.py" "PydroXL_19\envs\Pydro367\site-packages\sphinx\ext\inheritance_diagram.py" 


Then in a console in this directory run sphinx-build -M html . . -E -a
  or 
run make_docs.bat 
  whichs calls activate_pydro.bat   ==  call "%~dp0..\..\..\..\..\..\..\scripts\activate" Pydro367
  then calls compile_html.bat   ==  sphinx-build -M html "%~dp0." "%~dp0." -E -a