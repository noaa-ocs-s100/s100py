Docs need a python 3 environment with sphinx (2.4.4), sphinx-automodapi (0.12), sphinx-autodoc-typehints and python-graphviz (0.13.2) which will automatically get graphviz too (really you only need graphviz).
Also there is a modification for sphinx-automod 0.12
Most of the magic is in the conf.py in this directory to set up automodpai

conda install sphinx sphinx-automodapi sphinx-autodoc-typehints python-graphviz
  or (the current install at time of writing)
conda install sphinx=2.4.4 sphinx-automodapi=0.12 python-graphviz=0.13.2 sphinx-autodoc-typehints=1.10.3 graphviz=2.38
  then
copy "PydroXL_19\NOAA\site-packages\Python3\s100py\docs\automodsumm.py" "PydroXL_19\envs\Pydro367\site-packages\sphinx_automodapi\automodsumm.py" 


Then in a console in this directory run sphinx-build -M html . . -E -a
  or 
run make_docs.bat 
  whichs calls activate_pydro.bat   ==  call "%~dp0..\..\..\..\..\..\..\scripts\activate" Pydro367
  then calls compile_html.bat   ==  sphinx-build -M html "%~dp0." "%~dp0." -E -a