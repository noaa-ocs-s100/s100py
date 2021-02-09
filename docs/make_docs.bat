rd /s /q %~dp0_build
rd /s /q %~dp0api
sphinx-build -T -E -b html -d _build/doctrees -D language=en %~dp0 %~dp0_build/html