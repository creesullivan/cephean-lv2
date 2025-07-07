xcopy "C:\cephean-lv2\shared" "C:\cephean-lv2\dwarf\mod-plugin-builder\plugins\package\cephean-build\shared" /i /y
del "C:\cephean-lv2\dwarf\mod-plugin-builder\plugins\package\cephean-build\shared\cephean-test.h" /q
del "C:\cephean-lv2\dwarf\mod-plugin-builder\plugins\package\cephean-build\shared\cephean-test.cpp" /q
rmdir "C:\cephean-lv2\dwarf\mod-plugin-builder\plugins\package\cephean-build\source" /s /q
mkdir "C:\cephean-lv2\dwarf\mod-plugin-builder\plugins\package\cephean-build\source"
copy "C:\cephean-lv2\projects\%1\deploy\Makefile" "C:\cephean-lv2\dwarf\mod-plugin-builder\plugins\package\cephean-build\source\Makefile" /y
copy "C:\cephean-lv2\projects\%1\deploy\Makefile.mk" "C:\cephean-lv2\dwarf\mod-plugin-builder\plugins\package\cephean-build\source\Makefile.mk" /y
xcopy "C:\cephean-lv2\projects\%1\src" "C:\cephean-lv2\dwarf\mod-plugin-builder\plugins\package\cephean-build\source" /y
del "C:\cephean-lv2\dwarf\mod-plugin-builder\plugins\package\cephean-build\source\%1_test.cpp" /q
docker start mpbc
docker exec mpbc bash -c "source local.env moddwarf-new; cd plugins/package/cephean-build/source; make"
copy "C:\cephean-lv2\dwarf\mod-plugin-builder\plugins\package\cephean-build\source\%1.so" "C:\cephean-lv2\plugins\package\%1.lv2\%1.so" /y
docker stop mpbc