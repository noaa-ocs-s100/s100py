call del /Q /S "%~dp0.\api\*.*"
sphinx-build -M html "%~dp0." "%~dp0." -E -a