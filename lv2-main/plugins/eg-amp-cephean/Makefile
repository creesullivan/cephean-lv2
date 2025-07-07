#!/usr/bin/make -f

hellow: hellow.o

hellow.o : hellow.cpp cc -c 

# --------------------------------------------------------------
# eg-amp build rules

eg-amp: eg-amp.lv2/amp$(LIB_EXT) eg-amp.lv2/manifest.ttl

eg-amp.lv2/amp$(LIB_EXT): amp.c
	$(CC) $^ $(BUILD_C_FLAGS) $(LINK_FLAGS) -lm $(SHARED) -o $@

eg-amp.lv2/manifest.ttl: eg-amp.lv2/manifest.ttl.in
	sed -e "s|@LIB_EXT@|$(LIB_EXT)|" $< > $@

# --------------------------------------------------------------

clean:
	rm -f eg-amp.lv2/amp$(LIB_EXT) eg-amp.lv2/manifest.ttl

# --------------------------------------------------------------
