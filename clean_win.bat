@ECHO OFF
SET Files=datasets\Windows\*.asn;temp\*.*

ECHO Deleting the following files: %Files%

FOR %%Y IN (%Files%) DO (
ECHO   Deleting "%%Y"
DEL /Q "%%Y"
)

ECHO Done.
REM PAUSE

