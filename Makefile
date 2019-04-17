all: src

src:
	$(MAKE) -C src
test:
	$(MAKE) -C src test
format:
	astyle --options=astyle.conf src/*.cpp,*.h
.PHONY: all src test clean
clean:
	$(MAKE) -C src clean
	
