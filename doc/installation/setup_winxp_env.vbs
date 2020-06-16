'************************************************* ***********
'* Modifying The System Path With New Entries *
'************************************************* ***********
Dim ExistingPath, NewPath
Set oShell = WScript.CreateObject("WScript.Shell")
Set oEnv = oShell.Environment("SYSTEM")

'************************************************* ***********
'* Add your Path Entry Here *
'************************************************* ***********
ExistingPath = oEnv("PATH")
NewPath = ExistingPath & ";" & "C:\Python25" & ";" & "C:\Program Files\swigwin-1.3.33"
oEnv("PATH") = NewPath

oEnv("PYTHON_PARAM") = "C:\Python25\libs\python25.lib"
oEnv("PYTHON_INCLUDE") = "C:\Python25\include"