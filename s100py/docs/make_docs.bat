call activate_pydro.bat
rem call "%~dp0..\..\..\..\..\..\..\scripts\activate" Pydro367

rem This is now done automatically by automodapi extension, which is set up in conf.py
rem call generate_rst.bat
rem rem call sphinx-apidoc -o "%~dp0." "%~dp0.."

call compile_html.bat
rem sphinx-build -M html "%~dp0." "%~dp0." -E -a