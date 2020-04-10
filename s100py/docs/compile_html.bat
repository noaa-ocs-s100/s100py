call del /Q /S "%~dp0.\api\*.*"
sphinx-build -M html "%~dp0." "%~dp0." -E -a
PUSHD "%~dp0.\html"
python "%~dp0make_img_maps_resizable.py"
POPD