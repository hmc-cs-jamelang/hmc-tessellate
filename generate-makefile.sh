#!/usr/bin/env bash

stlib="$1"

if [ -z "$stlib" ]
then
	echo "Enter path for stlib (https://bitbucket.org/seanmauch/stlib): "
	read stlib
fi

cat > Makefile <<END
# Manually configured
STLIB_INCLUDE_PATH = "$stlib"
END

cat >> Makefile <<'END'
# Generic
CC = g++
CC_INCLUDE = -I$(STLIB_INCLUDE_PATH)
CC_FLAGS = -std=c++11 -O3

TARGETS = test

.PHONY: all clean again

all: $(TARGETS)

clean:
	rm -f $(TARGETS)

again: clean $(TARGETS)

test: test.cpp
	$(CC) $< -o $@ $(CC_INCLUDE) $(CC_FLAGS)
END

exit
