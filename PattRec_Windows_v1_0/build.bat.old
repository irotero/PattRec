@echo on
setlocal
set path=C:\Program Files\Java\jdk1.8.0_171\bin;%path%
if not exist build md build
if not exist build\classes md build\classes
javac -source 1.8 -target 1.8 -classpath "lib\commons-net-3.6.jar;lib\AbsoluteLayout.jar;lib\mysql-connector-java-5.1.23-bin.jar;lib\opencsv-3.0.jar;lib\picard.jar;lib\JRI.jar" -implicit:class -d build\classes -sourcepath src src\org\gradiant\main.java
copy src\org\gradiant\UI\*.form build\classes\org\gradiant\UI\ > NUL
cd build\classes
jar cmf ..\..\manifest\MANIFEST.MF ..\..\gridd.jar org
cd ..\..
